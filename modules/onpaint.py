import numpy as np
from scipy.signal import medfilt 
from aperture_contour import draw_contours
from aperture_contour import make_contours

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

    elif event.key=="C":
        crossbkgmask=mask.crossbkgmask()
        if swicb[0]:
            crossbkgmask=(crossbkgmask).T
            swicb[0]=False
        else:
            swicb[0]=-1
            swicb[0]=True
        cbm=crossbkgmask.reshape(nx*ny).astype(np.bool)
        icbm=np.invert(cbm)
                    
        for i in range(len(masklc)):
            masklc[i]=cbm[i] and maskbkgc[i]
        for i in range(len(maskbkg)):
            maskbkg[i]=icbm[i] and maskbkgc[i]
        

    ######### LC PANEL #####
    #setting pixel panel
    #pk m/M d/D y/Y t 1 b l h B
    if event.key == "p":
        ax2.clear()
        ax.clear()
        ind = np.searchsorted(time[mask], event.xdata, side='left')
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.axvline(time[mask][ind],ls="dotted",color="green")
        ax.title.set_text("Slicing mode. t="+str(round(time[mask][ind]))+" File: "+tagn)

        ax2.imshow(np.log10(cnts[mask][ind]))
        ax2.title.set_text('Log color Image at t='+str(round(time[mask][ind])))
    elif event.key == "k":
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        ax.clear()            
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        
        for j,iq in enumerate(np.unique(quality)): 
            qmask=quality==iq
            ax.plot(time[qmask],flux[qmask]/np.median(flux[mask]),".",color="C"+str(j),label="Q="+str(iq))
        ax.legend()


    elif event.key == "m":
        ax2.clear()
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
        ax2.clear()
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
        ax2.clear()
           
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
        ax2.clear()

        ind = pind[0]+1
        pind[0] = ind
        ax.plot(time[mask],flux[mask]/np.median(flux[mask]),".")
        ax.axvline(time[mask][ind],ls="dotted",color="red")
        ax.title.set_text("Moving Difference from Median. Put D to move. File: "+tagn)

        diff = cnts[mask][ind]-medcnts

#        maskW=mask.adaptive_bkgmask(diff)
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
        ax2.clear()

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
        ax2.clear()

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
        ax2.clear()

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
        ax2.clear()

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
        #### smooth bkg XXX ###
        smtbkg=medfilt(bkg,kernel_size=151)
        lc = lc - bkg/nbkg*nlc
        #lc = lc - smtbkg/nbkg*nlc
        
        ax4.plot(tpft[mask],bkg[mask],".",label="background")
        ax4.plot(tpft[mask],smtbkg[mask],".",label="background")

        ax4.legend()
    if  nlc>0:
        ax3.clear()
        ax3.plot(tpft[mask],lc[mask],".",label="cnts")
        ax3.legend()

    if event.key == "b":
        xlim=ax.get_xlim()
        ylim=ax.get_ylim()
        ax.clear()

        ax.plot(time[mask],flux[mask]/np.abs(np.median(flux[mask])),".",label="original")
        ax.plot(time[mask],lc[mask]/np.abs(np.median(lc[mask])),".",label="new aperture")
#        ax.set_xlim(xlim)
#        ax.set_ylim(ylim)

        ax.legend()
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
