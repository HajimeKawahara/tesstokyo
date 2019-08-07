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

import sys

def numericalSort(value):
    numbers = re.compile(r'(\d+)')
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

def different_centroid():
    Xp=(np.array(range(ny))*(np.array([np.ones(nx)]).T)).reshape(nx*ny)
    Yp=((np.array(range(nx))*(np.array([np.ones(ny)]).T)).T).reshape(nx*ny)

    diff = cntsf[:,:]-medcntsf
    diffbkg = cntsf[:,maskbkgc]-medcntsf[maskbkgc]

    print(np.shape(Xp))
    cX=np.sum(Xp*diff,axis=1)/np.sum(diff,axis=1)
    cY=np.sum(Yp*diff,axis=1)/np.sum(diff,axis=1)

    cXp=np.sum(Xp[maskbkgc]*diffbkg,axis=1)/np.sum(diffbkg,axis=1)
    cYp=np.sum(Yp[maskbkgc]*diffbkg,axis=1)/np.sum(diffbkg,axis=1)

    medcX=np.median(cX)
    medcY=np.median(cY)
    ccXY = np.sqrt((cX - medcX)**2 + (cY - medcY)**2)
    ccXY=ccXY/np.median(ccXY)
    
    medcXp=np.median(cXp)
    medcYp=np.median(cYp)
    ccXYp = np.sqrt((cXp - medcXp)**2 + (cYp - medcYp)**2)
    ccXYp=ccXYp/np.median(ccXYp)
    return ccXYp

def explanation():
    exp=[]
    exp.append("   ")
    exp.append("************************************")
    exp.append("qlh5: Quick Analyzer of TESS FFI H5")
    exp.append("************************************")
    exp.append("   ")
    exp.append("PRESS KEY IN THE LEFT BOTTOM PANEL")
    exp.append("z: use as aperture (Z: use default)")
    exp.append("x: use as background (X: use default)")
    exp.append("v: do not use")
    exp.append("b: compare")
    exp.append("l: switch log/linear color")
    exp.append("1: smoothing")
    exp.append("   ")
    exp.append("PRESS KEY IN THE UPPER PANEL")
    exp.append("p: pixel image at location")
    exp.append("d: pixel difference image at location from the median image (D: move to right)")
    exp.append("y: pixel difference image at two locations. put twice. (Y: move to right)")
    exp.append("t: pixel difference image at two ranges. put 2 (inside) 2+2 (outside) times.")
    exp.append("m: nearest pixel difference image at location (M: move to right)")
    exp.append("   ")
    exp.append("-------------------------------------")
    return exp

def adaptive_bkgmask(diff,crit=3):
    maskWW=np.abs((diff - np.median(diff))/np.std(diff)) > crit
    maskW = np.copy(maskWW)
    maskW[1:,:]=maskWW[:-1,:]+maskW[1:,:]
    maskW[:-1,:]=maskWW[1:,:]+maskW[:-1,:]
    maskW[:,1:]=maskWW[:,:-1]+maskW[:,1:]
    maskW[:,:-1]=maskWW[:,1:]+maskW[:,:-1]
    
    maskW[1:,:-1]=maskWW[:-1,1:]+maskW[1:,:-1]
    maskW[:-1,1:]=maskWW[1:,:-1]+maskW[:-1,1:]
    maskW[1:,1:]=maskWW[:-1,:-1]+maskW[1:,1:]
    maskW[:-1,:-1]=maskWW[1:,1:]+maskW[:-1,:-1]        
    return maskW

def oncpaint(event):


    ix=int(event.xdata+0.5)
    iy=int(event.ydata+0.5)
    if event.key == "z":
        masklc[ix+iy*nx]=True
        maskbkg[ix+iy*nx]=False
    elif event.key == "x":
        maskbkg[ix+iy*nx]=True
        masklc[ix+iy*nx]=False
    elif event.key == "v":
        masklc[ix+iy*nx]=False
        maskbkg[ix+iy*nx]=False
    elif event.key == "Z":
        for i in range(len(masklc)):
            masklc[i]=masklcc[i]
            if masklc[i]:
                maskbkg[i]=False
    elif event.key == "X":
        for i in range(len(maskbkg)):
            maskbkg[i]=maskbkgc[i]
            if maskbkgc[i]:
                masklc[i]=False
#        masklc[maskbkgc]=False
        
    #setting pixel panel
    ax2.clear()
    if event.key == "p":
        ax.clear()
        ind = np.searchsorted(time[mask], event.xdata, side='left')
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.axvline(time[mask][ind],ls="dotted",color="green")
        ax.title.set_text("Slicing mode. t="+str(round(time[mask][ind]))+" File: "+tagn)

        ax2.imshow(np.log10(cnts[mask][ind]))
        ax2.title.set_text('Log color Image at t='+str(round(time[mask][ind])))
    elif event.key == "m":
        ax.clear()            
        ind = np.searchsorted(time[mask_qzero], event.xdata, side='left')
        ax.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),".",color="gray",alpha=0.3)
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")

        ax.axvline(time[mask_qzero][ind],ls="dotted",color="red")
        ax.title.set_text("Difference from Median. t="+str(round(time[mask_qzero][ind]))+" File: "+tagn)

        c=ax2.imshow(cnts[mask_qzero][ind]-cnts[mask_qzero][ind-1])
        ax2.title.set_text('Nearest Difference Image at t='+str(round(time[mask_qzero][ind],3)))

        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        try:
            cbar[0].remove()
            plt.draw()
        except:
            print("Set colorbar.")
        cbar[0] = plt.colorbar(c,ax=ax2,cax=cax)
        pind[0]=ind

    elif event.key == "M":
        ax.clear()            
        ind = pind[0]+1
        pind[0] = ind
        ax.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),".",color="gray",alpha=0.3)
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.axvline(time[mask_qzero][ind],ls="dotted",color="red")
        ax.title.set_text("Moving Nearest Difference. Put M to move. File: "+tagn)
        diff = cnts[mask_qzero][ind]-cnts[mask_qzero][ind-1]        
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        c=ax2.imshow(diff)
        try:
            cbar[0].remove()
            plt.draw()
        except:
            print("Set colorbar.")
        cbar[0] = plt.colorbar(c,ax=ax2,cax=cax)
        ax2.title.set_text('Difference Image at t='+str(round(time[mask_qzero][ind],3)))

        
    elif event.key == "d":
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        ax.clear()            
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

        ind = np.searchsorted(time[mask], event.xdata, side='left')
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.axvline(time[mask][ind],ls="dotted",color="red")
        ax.title.set_text("Difference from Median. t="+str(round(time[mask][ind],3))+" File: "+tagn)

        c=ax2.imshow(cnts[mask][ind]-medcnts)
        ax2.title.set_text('Difference Image at t='+str(round(time[mask][ind],3)))
        pind[0]=ind
    elif event.key == "D":
        ax.clear()            
        ind = pind[0]+1
        pind[0] = ind
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.axvline(time[mask][ind],ls="dotted",color="red")
        ax.title.set_text("Moving Difference from Median. Put D to move. File: "+tagn)

        diff = cnts[mask][ind]-medcnts

#        maskW=adaptive_bkgmask(diff)
#        apbkgx=np.copy(apbkg)
#        apbkgx[maskW]=False
#        ax6.clear()
#        ax6.imshow(ap,alpha=0.3, cmap="gray")
#        ax6.imshow(apbkgx,alpha=0.3, cmap="bwr")
#        ax6.imshow(maskW,alpha=0.1, cmap="bwr")
#        ax6.title.set_text("adaptive aperture")
        
        divider = make_axes_locatable(ax2)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        c=ax2.imshow(diff)
        try:
            cbar[0].remove()
            plt.draw()
        except:
            print("Set colorbar.")
        cbar[0] = plt.colorbar(c,ax=ax2,cax=cax)

        ax2.title.set_text('Difference Image at t='+str(round(time[mask][ind])))
    elif event.key == "y":
        if diffind[0] < 0:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            diffind[0]=np.searchsorted(time[mask], event.xdata, side='left')
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvline(time[mask][diffind[0]],ls="dashed",color="red")
            ax.title.set_text("Difference mode. Put y again for comparison. File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax2.imshow(np.log10(cnts[mask][diffind[0]]))
        else:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()

            diffind[1]=np.searchsorted(time[mask], event.xdata, side='left')
            c=ax2.imshow(cnts[mask][diffind[1]]-cnts[mask][diffind[0]])
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            try:
                cbar[0].remove()
                plt.draw()
            except:
                print("Set colorbar.")
            cbar[0] = plt.colorbar(c,ax=ax2,cax=cax)
            
            ax2.title.set_text('Difference Image t='+str(round(time[mask][diffind[0]],3))+" from t="+str(round(time[mask][diffind[1]],3)))
            ax.clear()            
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvline(time[mask][diffind[0]],color="red")            
            ax.axvline(time[mask][diffind[1]],ls="dashed",color="blue")
            ax.title.set_text("Difference of t="+str(round(time[mask][diffind[0]],3))+" from t="+str(round(time[mask][diffind[1]],3))+" File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            diffind[0]=-1
            diffind[1]=-1
    elif event.key == "Y":
        if diffind[0] == -1:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            diffind[0]=np.searchsorted(time[mask], event.xdata, side='left')
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvline(time[mask][diffind[0]],color="red")
            ax2.imshow(np.log10(cnts[mask][diffind[0]]))
        elif diffind[1] == -1:
            diffind[1]=np.searchsorted(time[mask], event.xdata, side='left')
            ax2.imshow(cnts[mask][diffind[1]]-cnts[mask][diffind[0]])
            ax2.title.set_text('Difference Image t='+str(round(time[mask][diffind[0]],3))+" from t="+str(round(time[mask][diffind[1]],3)))
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvline(time[mask][diffind[0]],color="red")            
            ax.axvline(time[mask][diffind[1]],ls="dashed",color="blue")
            pdiffind[0]=np.copy(diffind[0])
            diffind[0]=-2
        else:            
            diffind[1]=np.copy(diffind[1])+1
            ax2.imshow(cnts[mask][diffind[1]]-cnts[mask][pdiffind[0]])
            ax2.title.set_text('Difference Image t='+str(round(time[mask][pdiffind[0]],3))+" from t="+str(round(time[mask][diffind[1]],3)))
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvline(time[mask][pdiffind[0]],color="red")            
            ax.axvline(time[mask][diffind[1]],ls="dashed",color="blue")
    elif event.key == "t":
        if diffrangeind[0] < 0:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            diffrangeind[0]=np.searchsorted(time[mask], event.xdata, side='left')
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvline(time[mask][diffrangeind[0]],ls="dashed",color="red")
            ax.title.set_text("Difference Range mode. Put t  again  to determine the inside range. File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax2.imshow(np.log10(cnts[mask][diffrangeind[0]]))
        elif diffrangeind[1] < 0:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            diffrangeind[1]=np.searchsorted(time[mask], event.xdata, side='left')
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvspan(time[mask][diffrangeind[0]],time[mask][diffrangeind[1]],alpha=0.2,color="red")
            ax.title.set_text("Difference Range mode. Put t again to determine the outside range 1 starting point. File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax2.imshow(np.log10(cnts[mask][diffrangeind[1]]))
        elif diffrangeind[2] < 0:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            diffrangeind[2]=np.searchsorted(time[mask], event.xdata, side='left')
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvspan(time[mask][diffrangeind[0]],time[mask][diffrangeind[1]],alpha=0.2,color="red")
            ax.axvline(time[mask][diffrangeind[2]],ls="dashed",color="blue")
            ax.title.set_text("Difference Range mode. Put t again to determine the outside range 1 ending point. File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax2.imshow(np.log10(cnts[mask][diffrangeind[2]]))
        elif diffrangeind[3] < 0:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            diffrangeind[3]=np.searchsorted(time[mask], event.xdata, side='left')

            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvspan(time[mask][diffrangeind[0]],time[mask][diffrangeind[1]],alpha=0.2,color="red")
            ax.axvspan(time[mask][diffrangeind[2]],time[mask][diffrangeind[3]],alpha=0.2,color="blue")
            ax.title.set_text("Difference Range mode. Put t again to determine the outside range 2 starting point. File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax2.imshow(np.log10(cnts[mask][diffrangeind[3]]))

        elif diffrangeind[4] < 0:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()
            ax.clear()            
            diffrangeind[4]=np.searchsorted(time[mask], event.xdata, side='left')
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvspan(time[mask][diffrangeind[0]],time[mask][diffrangeind[1]],alpha=0.2,color="red")
            ax.axvspan(time[mask][diffrangeind[2]],time[mask][diffrangeind[3]],alpha=0.2,color="blue")
            ax.axvline(time[mask][diffrangeind[4]],ls="dashed",color="blue")
            ax.title.set_text("Difference Range mode. Put t again to determine the outside range 2 starting point. File: "+tagn)
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)

            ax2.imshow(np.log10(cnts[mask][diffrangeind[4]]))
        else:
            xlim=ax.get_xlim()
            ylim=ax.get_ylim()

            diffrangeind[5]=np.searchsorted(time[mask], event.xdata, side='left')
            iimax=np.max([diffrangeind[1],diffrangeind[0]])
            iimin=np.min([diffrangeind[1],diffrangeind[0]])
            inside = np.median(cnts[mask][iimin:iimax],axis=0)

            iimax=np.max([diffrangeind[3],diffrangeind[2]])
            iimin=np.min([diffrangeind[3],diffrangeind[2]])

            jimax=np.max([diffrangeind[5],diffrangeind[4]])
            jimin=np.min([diffrangeind[5],diffrangeind[4]])
            outside = np.median(np.concatenate([cnts[mask][iimin:iimax],cnts[mask][jimin:jimax]]),axis=0)
            c=ax2.imshow(outside-inside)
            divider = make_axes_locatable(ax2)
            cax = divider.append_axes("right", size="5%", pad=0.05)
            try:
                cbar[0].remove()
                plt.draw()
            except:
                print("Set colorbar.")
            cbar[0] = plt.colorbar(c,ax=ax2,cax=cax)
            
            ax2.title.set_text('Difference Range Image')
            ax.clear()            
            ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
            ax.axvspan(time[mask][diffrangeind[0]],time[mask][diffrangeind[1]],alpha=0.2,color="red")
            ax.axvspan(time[mask][diffrangeind[2]],time[mask][diffrangeind[3]],alpha=0.2,color="blue")
            ax.axvspan(time[mask][diffrangeind[4]],time[mask][diffrangeind[5]],alpha=0.2,color="blue")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            for i in range(0,6):
                diffrangeind[i]=-1

            
            
    else:
        if event.key=="l" and not logsw[0]:
            logsw[0]=True
        elif event.key=="l" and logsw[0]:
            logsw[0]=False

        if logsw[0]:
            ax2.imshow(np.log10(medcnts))
        else:
            ax2.imshow(medcnts)

        try:
            points_list = make_contours(masklc.reshape(ny,nx))
            for points in points_list:
                for i in range(len(points) - 1):
                    ax2.plot((points[i][0], points[i+1][0]), (points[i][1], points[i+1][1]),  "-",color="orange",lw=1.3,ls="dashed")
        except:
            lcind=np.where(masklc.reshape(ny,nx) == True)
            ax2.plot(lcind[1],lcind[0],"*",color="cyan")

        try:
            points_list = make_contours(maskbkg.reshape(ny,nx))
            for points in points_list:
                for i in range(len(points) - 1):
                    ax2.plot((points[i][0], points[i+1][0]), (points[i][1], points[i+1][1]),  "-",color="green",lw=2,ls="dotted")
        except:
            bkgind=np.where(maskbkg.reshape(ny,nx) == True)
            ax2.plot(bkgind[1],bkgind[0],"X",color="red")


    nlc=np.sum(masklc)
    nbkg=np.sum(maskbkg)

    lc=np.sum(cntsf[:,masklc],axis=1)
    if nbkg>0:
        ax4.clear()
        bkg=np.sum(cntsf[:,maskbkg],axis=1)
        lc = lc - bkg/nbkg*nlc
        ax4.plot(tpft[mask],bkg[mask],".",label="background")
        ax4.legend()
    if  nlc>0:
        ax3.clear()
        ax3.plot(tpft[mask],lc[mask],".",label="cnts")
        ax3.legend()

    if event.key == "b":
        ax.clear()
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.plot(time[mask],lc[mask]/np.median(lc[mask]),".")
        ax.title.set_text("STD = "+str(np.std(lc[mask]/np.median(lc[mask]))))
    if event.key == "1":
        ax.clear()
        fsmt=medfilt(flux[mask]/np.median(flux[mask]),kernel_size=31)
        ax.plot(time[mask],flux[mask]/np.median(flux[mask])/fsmt,".")
        if nbkg > 0:
            fsmtx=medfilt(lc[mask]/np.median(lc[mask]),kernel_size=51)
            ax.plot(time[mask],lc[mask]/np.median(lc[mask])/fsmtx,".")


    if event.key == "h":

        tip=0.06
        tip2=0.05
        for i,line in enumerate(expl[6:13]):
            ax2.annotate(line, xy=(0.1,0.9-tip2*i), xycoords='axes fraction', fontsize=12,color="white",bbox=dict(boxstyle='round,pad=0.2', color="None",fc='black', alpha=0.3))
        for i,line in enumerate(expl[14:-2]):
            ax.annotate(line, xy=(0.4,0.7-tip*i), xycoords='axes fraction', fontsize=12,color="green",bbox=dict(boxstyle='round,pad=0.2', color="None",fc='white', alpha=0.7))


    fig.canvas.draw()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='QL for TESS (h5)')
    parser.add_argument('-f', nargs=1, default=[""], help='tic', type=str)
    parser.add_argument('-n', help='do not connect mysql server', action='store_true')

    expl=explanation()
    for line in expl:
        print(line)
    args = parser.parse_args()    
    h5path = args.f[0]
    with h5py.File(h5path, "r") as f:
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


    ###########################################
    #
    # make lc from the original aperture
    #
    ###########################################
    masklcc=ap.reshape(nx*ny).astype(np.bool)
    maskbkgc=apbkg.reshape(nx*ny).astype(np.bool)
    nlcc=np.sum(masklcc)
    nbkgc=np.sum(maskbkgc)
    mask_qzero=(quality==0) + (quality==2)

    lcc=np.mean(cntsf[:,masklcc],axis=1)
    bkgc=np.mean(cntsf[:,maskbkgc],axis=1)
    flux = lcc - bkgc
    
    ######## median filtered BKG (not good)#######
    #    bkgcmed=medfilt(bkgc,kernel_size=21)
    #    flux = lcc - bkgcmed
    ######### ADAPTIVE BKG (not good) ######################
    #    flux=[]
    #    for i in range(0,nt):
    #        diff = cnts[i,:,:]-medcnts
    #        maskW=adaptive_bkgmask(diff)
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
    ndiff = cnts[1:,:,:]-cnts[:-1,:,:]
    ndiff = np.concatenate([np.zeros((1,nx,ny)),ndiff])
    ndiff = ndiff.reshape(nt,nx*ny)
    ndmax = np.max(np.abs(ndiff[:,maskbkgc]),axis=1)
    ndmax = ndmax/medfilt(ndmax,kernel_size=51)

    Q2 = stats.scoreatpercentile(ndmax, 50+34.1)
    Q1 = stats.scoreatpercentile(ndmax, 50-34.1)
    Qsigma = (Q2-Q1)/2.0
    crit = 4.0
    mask_asteroid = (ndmax - np.nanmedian(ndmax)) < crit*Qsigma
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
    mask = mask_qzero*mask_asteroid


    print("FLUX ORIGINAL STD",np.std(flux_orig[mask]/np.median(flux_orig[mask])))
    print("FLUX STD",np.std(flux[mask]/np.median(flux[mask])))
    ###############################################
    masklc=np.zeros(nx*ny,dtype=bool)
    maskbkg=np.zeros(nx*ny,dtype=bool)
    logsw=np.ones(1,dtype=bool)
    pind=np.zeros(1,dtype=int)
    diffind=[-1,-1]
    pdiffind=[-1]
    diffrangeind=[-1,-1,-1,-1,-1,-1]

    ######################################

    fig=plt.figure(figsize=(25,10))
    ax=fig.add_subplot(211)
    ax.title.set_text("qlh5 (*_*)/ "+" File: "+tagn)
#    ax.plot(time,quality)
    ax.plot(time[mask],flux_orig[mask]/np.median(flux_orig[mask]),".",label="LC (original)",color="C1",alpha=0.2)
    ax.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),".",color="gray",alpha=0.5)
    ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".",label=tagn,color="C0")
    ax.title.set_text(tagn)
    ax.legend()

    ax2=fig.add_subplot(245)
    ax2.title.set_text("PRESS h for key binds.")
    ax2.imshow(np.log10(medcnts))

    ax3 = fig.add_subplot(426)
    ax3.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),".",color="black")
#    ax3.plot(time[mask_qzero],flux[mask_qzero]/np.median(flux[mask]),color="black")

#    ax3.plot(time[mask],flux[mask]/np.median(flux[mask]),".",label="flux")
    ax3.legend()

    ax4 = fig.add_subplot(428)
#    ax4.plot(time,ndsum,".",label="asteroid indicator")
    ax4.plot(time[mask_qzero],ndmax[mask_qzero],".",label="asteroid indicator (masked)",color="red")
    ax4.plot(time[mask],ndmax[mask],".",label="asteroid indicator",color="C0")

    ax4.legend()

    #explanation
    ax5=fig.add_subplot(4,4,14)
    ax5.text(0.1,0.8,tagn,color="C0",fontsize=12)
#    if True:
    try:
        print("Try to connect mysql server...")
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


    ax6=fig.add_subplot(4,8,20)
#    ax6.imshow(ap.astype(np.float)-apbkg.astype(np.float),alpha=0.3, cmap="bwr")
    ax6.imshow(np.log10(medcnts))#,cmap="gray")
    try:
        ax6 = draw_contours(ax6, ap, color="orange", lw=1.3,ls="dashed")
        ax6 = draw_contours(ax6, apbkg, color="green", lw=1.3,ls="dashed")
    except:
        ax6.imshow(ap.astype(np.float)-apbkg.astype(np.float),alpha=0.3, cmap="bwr")
        
    ax6.title.set_text("default aperture")

    ax7=fig.add_subplot(4,8,19)
    ax7.imshow(np.log10(medcnts))
    ax7.title.set_text("median (log)")

    cid = fig.canvas.mpl_connect('key_press_event', oncpaint)

    plt.show()
