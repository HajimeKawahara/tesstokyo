#/usr/bin/python
import pylab
import numpy as np
import matplotlib.pyplot as plt
import argparse
from astropy.io import fits
from photutils import MMMBackground
from astropy.stats import SigmaClip
import cv2
#from skimage.restoration import inpaint_biharmonic
import numpy.ma as ma
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='QL for TESS (FFIC)')
    parser.add_argument('-f', nargs=1, default=[""], help='file', type=str)
    parser.add_argument('-b', nargs=1, default=["lowpass"], help='bkgmode', type=str)

    
    args = parser.parse_args()
    bkgmode=args.b[0]

    fitspath = args.f[0]
    with fits.open(fitspath) as hdu:
        # time                                                              
        time_arr = hdu[1].header["TSTART"]
        # flux                                                              
        flux_arr = np.array(hdu[1].data)

    x=447
    y=1971
    ed=70
    miniarr=flux_arr[y-ed:y+ed,x-ed:x+ed]

    flux_arrraw=np.copy(flux_arr)

    #edge
    flux_arr[flux_arr<1]=np.median(flux_arr)
    if bkgmode=="elelike":
        #elenor-like bkg
        sigma_clip = SigmaClip(sigma=10.0)
        bkgc = MMMBackground(sigma_clip=sigma_clip)
        bkg = bkgc.calc_background(miniarr)
    elif bkgmode=="lowpass":
        flux_arr[flux_arr!=flux_arr]=np.median(flux_arr)
        flux_arr[flux_arr>3.0*np.median(flux_arr)]=np.median(flux_arr)
#        plt.imshow(np.log10(flux_arr))
        
        src = np.fft.fft2(flux_arr)
        h, w = src.shape
        cy, cx =  int(h/2), int(w/2)
        a=1.0/100
        rh, rw = int(a*cy), int(a*cx)
        fsrc =  np.fft.fftshift(src)
        fdst = np.zeros(src.shape, dtype=complex)
        fdst[cy-rh:cy+rh, cx-rw:cx+rw] = fsrc[cy-rh:cy+rh, cx-rw:cx+rw]
        fdst =  np.fft.fftshift(fdst)
        dst = np.fft.ifft2(fdst)
        bkg = (np.real(dst))
        
    corf=flux_arrraw-bkg
    corf=corf-np.min(corf)
    fig=plt.figure()
    ax=fig.add_subplot(121)
    plt.imshow(corf,vmin=200,vmax=1000)
    plt.axvline(x,color="white")
    plt.axhline(y,color="white")

    ax=fig.add_subplot(122)
    plt.imshow(bkg,vmin=200,vmax=700)
    plt.axvline(x,color="white")
    plt.axhline(y,color="white")
    plt.show()
