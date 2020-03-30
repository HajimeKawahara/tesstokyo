def different_centroid(nx,ny,cntsf,medcntsf,maskbkgc):
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
