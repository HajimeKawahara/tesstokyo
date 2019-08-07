import argparse
import pandas as pd
import tarfile
import os
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='pack tic files')
    parser.add_argument('-f', nargs=1, default=["pack.txt"], help='filename of the list of TICINT___', type=str)
    parser.add_argument('-d', nargs=1, default=["/pike/pipeline/"], help='dir', type=str)

    args = parser.parse_args()    

    filename=args.f[0]
    dat=pd.read_csv(filename,names=('TIC',))
    print("mkdir pack")
    for f in dat["TIC"]:        
        path=os.path.join(args.d[0],"step3/tess_"+f+".h5")
        if not os.path.exists(path):
            path=os.path.join(args.d[0],"TIC3/tess_"+f+".h5")

        print("cp ",path,"pack")
    print("tar cvfz pack.tar.gz pack")
    
