#/usr/bin/python
import pylab
import numpy as np
import matplotlib.pyplot as plt
import h5py
import argparse
from scipy.signal import medfilt 
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy import stats
from aperture_contour import draw_contours
from aperture_contour import make_contours
import glob
import os
import sys

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='VI for TESS (h5)')
    parser.add_argument('-d', nargs=1, default=[""], help='directory', type=str)
    parser.add_argument('-n', help='do not connect mysql server', action='store_true')

    args = parser.parse_args()    
    flist=glob.glob(os.path.join(args.d[0],"*.h5"))
    for h5path in flist:
        
        with h5py.File(h5path, "r") as f:
            print(h5path)

            try:
                TID=f["header"]["TID"].value
            except:
                TID=-1 #guest TID

                #        for k in (f["header"]):
                #            print(k)
            sector=f["header"]["sector"].value
            camera=f["header"]["camera"].value
            try:
                chip=f["header"]["chip"].value
            except:
                #OLD format
                chip=f["header"]["CCD"].value
        
            tagn=str(TID)+"_"+str(sector)+"_"+str(camera)+"_"+str(chip)
            #x = f["header"]["cx"].value
            #y = f["header"]["cy"].value
            time = f["LC"]["TIME"].value
            flux_orig = f["LC"]["SAP_FLUX"].value
            quality = f["LC"]["QUALITY"].value
            tpft=f["TPF"]["TIME"].value
            cnts=f["TPF"]["ROW_CNTS"].value
            medcnts=(np.median(cnts,axis=0))
            ap=f["APERTURE_MASK"]["FLUX"].value
            apbkg=f["APERTURE_MASK"]["FLUX_BKG"].value
            
            #######################################
            nt=cnts.shape[0]
            nx=cnts.shape[1]
            ny=cnts.shape[2]
            cntsf=cnts.reshape(nt,nx*ny)
            medcntsf=medcnts.reshape(nx*ny)
            cbar = [1]
            
            

        fig=plt.figure(figsize=(12,5))
        ax=fig.add_subplot(2,1,1)
        mask=quality==0
        ax.plot(time[mask],flux_orig[mask],".")
        #explanation
        ax5=fig.add_subplot(2,2,3)
        ax5.text(0.1,0.8,tagn,color="C0",fontsize=12)
        #    if True:
        try:
            print("Try to connect mysql server...")
            if TID>0 and not args.n:
                import mysql.connector
                from urllib.parse import urlparse
                url = urlparse('mysql://fisher:atlantic@133.11.231.118:3306/TESS')
                conn = mysql.connector.connect(
                    host = url.hostname or '133.11.231.118',
                    port = url.port or 3306,
                    user = url.username or 'fisher',
                    password = url.password or 'atlantic',
                    database = url.path[1:],
                )
                cur = conn.cursor()
                out=[]
                datab=["TICv7s_has_key","TICv7n_has_key"]
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
        
        
        ax6=fig.add_subplot(2,2,4)
        #    ax6.imshow(ap.astype(np.float)-apbkg.astype(np.float),alpha=0.3, cmap="bwr")
        ax6.imshow(np.log10(medcnts))#,cmap="gray")
        try:
            ax6 = draw_contours(ax6, ap, color="orange", lw=1.3,ls="dashed")
            ax6 = draw_contours(ax6, apbkg, color="green", lw=1.3,ls="dashed")
        except:
            ax6.imshow(ap.astype(np.float)-apbkg.astype(np.float),alpha=0.3, cmap="bwr")
            
        plt.savefig("vi/vi"+str(TID)+".png")
