#/usr/bin/python
import pylab
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
from scipy.signal import medfilt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
import sys

import modules.explanation as exp
import modules.misc as misc
import modules.maskmod as maskmod
import modules.diffcen as diffcen
import modules.onpaint as onp

import builtins

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='QL for TESS (h5)')
    parser.add_argument('-f', nargs=1, default=[""], help='h5 file', type=str)
    parser.add_argument('-n', help='do not connect mysql server', action='store_true')
    swicb=[False]
    builtins.expl=exp.explanation_qlh5()
    
    for line in expl:
        print(line)
    args = parser.parse_args()    
    h5path = args.f[0]
    with h5py.File(h5path, "r") as f:
        try:
            builtins.TID=f["header"]["TID"].value
        except:
            builtins.TID=-1 #guest TID
#        for k in (f["header"]):
#            print(k)
        builtins.sector=f["header"]["sector"].value
        builtins.camera=f["header"]["camera"].value
        try:
            builtins.chip=f["header"]["chip"].value
        except:
            #OLD format
            builtins.chip=f["header"]["CCD"].value
        
        builtins.ffix0 = f["header"]["x"].value
        builtins.ffiy0 = f["header"]["y"].value
        builtins.tagn=str(TID)+"_"+str(sector)+"_"+str(camera)+"_"+str(chip)+" FFI: x="+str(ffix0)+" y="+str(ffiy0)

        builtins.time = f["LC"]["TIME"].value
        builtins.flux_orig = f["LC"]["SAP_FLUX"].value
        builtins.quality = f["LC"]["QUALITY"].value

        builtins.tpft=f["TPF"]["TIME"].value
        builtins.cnts=f["TPF"]["ROW_CNTS"].value
        builtins.medcnts=(np.median(cnts,axis=0))
        builtins.ap=f["APERTURE_MASK"]["FLUX"].value
        builtins.apbkg=f["APERTURE_MASK"]["FLUX_BKG"].value

        #### ap position
        

        
    #######################################
    builtins.nt=cnts.shape[0]
    builtins.nx=cnts.shape[1]
    builtins.ny=cnts.shape[2]
    builtins.cntsf=cnts.reshape(nt,nx*ny)
    builtins.medcntsf=medcnts.reshape(nx*ny)
    builtins.cbar = [1]


    ###########################################
    #
    # make lc from the original aperture
    #
    ###########################################
    builtins.masklcc=ap.reshape(nx*ny).astype(np.bool)
    builtins.maskbkgc=apbkg.reshape(nx*ny).astype(np.bool)
    builtins.nlcc=np.sum(masklcc)
    builtins.nbkgc=np.sum(maskbkgc)
    builtins.mask_qzero=(quality==0) + (quality==2)

    builtins.lcc=np.mean(cntsf[:,masklcc],axis=1)
    builtins.bkgc=np.mean(cntsf[:,maskbkgc],axis=1)
    builtins.flux = lcc - bkgc
    
    ######## median filtered BKG (not good)#######
    #    bkgcmed=medfilt(bkgc,kernel_size=21)
    #    flux = lcc - bkgcmed
    ######### ADAPTIVE BKG (not good) ######################
    #    flux=[]
    #    for i in range(0,nt):
    #        diff = cnts[i,:,:]-medcnts
    #        maskW=maskmod.adaptive_bkgmask(diff)
    #        apbkgx=np.copy(apbkg)
    #        apbkgx[maskW]=False
    #        maskbkgx=apbkgx.reshape(nx*ny).astype(np.bool)
    #       flux.append(lcc[i] - np.median(cntsf[i,maskbkgx]))
    #   flux=np.array(flux)

    ##############################################
    #
    #　Asteroid Indicator
    #
    ##############################################
    builtins.ndiff = cnts[1:,:,:]-cnts[:-1,:,:]
    ndiff = np.concatenate([np.zeros((1,nx,ny)),ndiff])
    ndiff = ndiff.reshape(nt,nx*ny)
    builtins.ndmax = np.max(np.abs(ndiff[:,maskbkgc]),axis=1)
    ndmax = ndmax/medfilt(ndmax,kernel_size=51)

    Q2 = stats.scoreatpercentile(ndmax, 50+34.1)
    Q1 = stats.scoreatpercentile(ndmax, 50-34.1)
    Qsigma = (Q2-Q1)/2.0
    crit = 4.0
    builtins.mask_asteroid = (ndmax - np.nanmedian(ndmax)) < crit*Qsigma
    ##############################################
    #
    #  EB in Background Aperture
    #
    ##############################################


    ##############################################
    #
    # set masks
    #
    #############################################
    builtins.mask = mask_qzero#*mask_asteroid


    print("FLUX ORIGINAL STD",np.std(flux_orig[mask]/np.median(flux_orig[mask])))
    print("FLUX STD",np.std(flux[mask]/np.median(flux[mask])))
    ###############################################
    builtins.masklc=np.zeros(nx*ny,dtype=bool)
    builtins.maskbkg=np.zeros(nx*ny,dtype=bool)
    builtins.logsw=np.ones(1,dtype=bool)
    builtins.pind=np.zeros(1,dtype=int)
    builtins.diffind=[-1,-1]
    builtins.pdiffind=[-1]
    builtins.diffrangeind=[-1,-1,-1,-1,-1,-1]

    ######################################

    builtins.fig=plt.figure(figsize=(25,10))
    builtins.ax=fig.add_subplot(211)
    ax.title.set_text("qlh5 (*_*)/ "+" File: "+tagn)
    ax.plot(time[mask],flux_orig[mask]/np.abs(np.median(flux_orig[mask])),".",label="LC (original)",color="C1",alpha=0.2)
    ax.plot(time[mask_qzero],flux[mask_qzero]/np.abs(np.median(flux[mask])),".",color="gray",alpha=0.5)
    ax.plot(time[mask],flux[mask]/np.abs(np.median(flux[mask])),".",label=tagn,color="C0")
    ax.title.set_text(tagn)
    ax.legend()

    builtins.ax2=fig.add_subplot(245)
    ax2.title.set_text("PRESS h for key binds.")
    ax2.imshow(np.log10(medcnts))

    builtins.ax3 = fig.add_subplot(426)
    ax3.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),".",color="black")
#    ax3.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),color="black")

#    ax3.plot(time[mask],flux[mask]/np.median(flux[mask]),".",label="flux")
    ax3.legend()

    builtins.ax4 = fig.add_subplot(428)
#    ax4.plot(time,ndsum,".",label="asteroid indicator")
    ax4.plot(time[mask_qzero],ndmax[mask_qzero],".",label="asteroid indicator (masked)",color="red")
    ax4.plot(time[mask],ndmax[mask],".",label="asteroid indicator",color="C0")

    ax4.legend()

    #explanation
    builtins.ax5=fig.add_subplot(4,4,14)
    ax5.text(0.1,0.8,tagn,color="C0",fontsize=12)
#    if True:
    try:
        print("Try to connect mysql server...")
        print("Press CTRL-C if you are outside the tess tokyo network.")

        if TID>0 and not args.n:
            import mysql.connector
            from urllib.parse import urlparse
            url = urlparse('mysql://fisher:atlantic@133.11.229.168:3306/TESS')
            conn = mysql.connector.connect(
                host = url.hostname or '133.11.229.168',
                port = url.port or 3306,
                user = url.username or 'fisher',
                password = url.password or 'atlantic',
                database = url.path[1:],
            )
            cur = conn.cursor()
            out=[]
            datab=["CTLv8_has_key"]
            i=0
        while len(out)==0:
            com='SELECT rad,mass,Teff,logg,Vmag,Hmag,ra,`dec`,Tmag FROM '+datab[i]+' where ID='+str(TID)
            print(com)
            cur.execute(com)
            out=cur.fetchall()
            i=i+1
        out=out[0]
        ax5.text(0.1,0.6,str(out[0])+"Rs,  "+str(out[1])+"Ms",fontsize=12)
        ax5.text(0.1,0.45,str(out[2])+"K,  logg="+str(out[3]),fontsize=12)
        ax5.text(0.1,0.3,"V="+str(out[4])+", H="+str(out[5])+", Tmag="+str(out[8]),fontsize=12)
        ax5.text(0.1,0.15,"RA="+str(out[6])+",   DEC="+str(out[7]),fontsize=12)

    except:
        ax5.text(0.1,0.5,"Information Unavailable.")

    ax5.xaxis.set_ticklabels([])
    ax5.yaxis.set_ticklabels([])
    ax5.axes.get_xaxis().set_ticks([])
    ax5.axes.get_yaxis().set_ticks([])


    builtins.ax6=fig.add_subplot(4,8,20)
#    ax6.imshow(ap.astype(np.float)-apbkg.astype(np.float),alpha=0.3, cmap="bwr")
    ax6.imshow(np.log10(medcnts))#,cmap="gray")
    try:
        ax6 = draw_contours(ax6, ap, color="orange", lw=1.3,ls="dashed")
        ax6 = draw_contours(ax6, apbkg, color="green", lw=1.3,ls="dashed")
    except:
        ax6.imshow(ap.astype(np.float)-apbkg.astype(np.float),alpha=0.3, cmap="bwr")
        
    ax6.title.set_text("default aperture")

    builtins.ax7=fig.add_subplot(4,8,19)
    ax7.imshow(np.log10(medcnts))
    ax7.title.set_text("median (log)")

    cid = fig.canvas.mpl_connect('key_press_event', onp.oncpaint)

    plt.show()
