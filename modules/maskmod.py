def crossbkgmask():
    crossbkgmask=np.array([[True,True,True,True,True,True,True,True,True,True,True,True,True],[True,True,True,True,True,True,True,True,True,True,True,True,True],[True,True,True,True,True,True,True,True,True,True,True,True,True],[True,True,True,True,True,True,True,True,True,True,True,True,True],[True,True,True,True,True,True,True,True,True,True,True,True,True],[True,True,True,True,True,True,True,True,True,True,True,True,True],[True,True,True,True,True,True,False,False,False,False,False,False,False],[False,False,False,False,False,False,False,False,False,False,False,False,False],[False,False,False,False,False,False,False,False,False,False,False,False,False],[False,False,False,False,False,False,False,False,False,False,False,False,False],[False,False,False,False,False,False,False,False,False,False,False,False,False],[False,False,False,False,False,False,False,False,False,False,False,False,False],[False,False,False,False,False,False,False,False,False,False,False,False,False]])
    return crossbkgmask

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
