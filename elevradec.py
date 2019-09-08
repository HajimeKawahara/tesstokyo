import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import get_sun
import matplotlib.pyplot as plt
from astroquery.simbad import Simbad
import argparse

import sys
#Seq|Star|RAJ2000|DEJ2000|Dist|n_Dist|Age|Hmag|SpType|MA|l_MB|MB|l_rho|rho|l_Per|Per|x_Per|l_Ecc|Ecc|l_acrit|acrit|Ref|SimbadName

telpos={}
telpos["Subaru"]=[19.828611,360.0-155.48055,4139.,-10.0]
telpos["Ishigaki"]=[24.3666,124.1333,197.,9.0]
telpos["Kiso"]=[35.7972222,137.6254166,1130.0,9.0]


def is_float_expression(s):
    try:
        val = float(s)
        return True
    except ValueError:
        return False
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='visibility')
    parser.add_argument('-n', nargs="+", help='name',type=str)
    parser.add_argument('-a', nargs=1, default=[40.0],help='maximum elevation [deg]',type=float)
    parser.add_argument('-p', nargs=3, help='Manual latitude/longitude/height(m) of the observatory.',type=float)
    parser.add_argument('-t', nargs=1, help='Telescope name Subaru,Ishigaki',type=str)

    parser.add_argument('-d', nargs=1, default=["2019-9-8"],help='observation date',type=str)
    parser.add_argument('-r', nargs=1, help='ra dec',type=str)

    args = parser.parse_args()
    radec = args.r[0]
    print(radec)
    
    maxalt = args.a[0]#30.0 #maximum altitude

    if args.p:
        lat=args.p[0]#19.0 + 49/60.+ 43/3600.0
        lon=args.p[1]#155.0 + 28/60.0 + 50/3600.0 
        height = args.p[2]#4139.0
    elif args.t:
        lat,lon,height,offt=telpos[args.t[0]]
    else:
        sys.exit("Provide the telescope name or position.")
        
    utcoffset = -offt*u.hour
    midlocal=Time(args.d[0])+(1*u.d)
    midnight = midlocal + utcoffset
    print("Midnight (UT)=",midnight)

    print(midlocal.iso)
    location=EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=height*u.m)
    delta_midnight = np.linspace(-12, 12, 1000)*u.hour
    times_obs = midnight + delta_midnight    
    frame=AltAz(obstime=midnight+delta_midnight,location=location)
    sunaltazs = get_sun(midnight+delta_midnight).transform_to(frame)

    #### OBS mode ####
    maskobs=np.ones(len(delta_midnight),dtype=bool)
    ichangepoint=np.argmin(sunaltazs.alt.degree)
#    nightmask=(sunaltazs.alt < -0*u.deg)
#    print(len(delta_midnight[nightmask]))

    nightmask=(sunaltazs.alt < -0*u.deg)*maskobs

    fig=plt.figure(figsize=(7,7))
    ax=fig.add_subplot(111)
    ic=0
    lsarr=["solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted","solid","dashed","dotted"]

    for i,name in enumerate(args.r):

        if len(radec.split(":"))>1:
            c = SkyCoord(radec, unit=(u.hourangle, u.deg))
        else:
            c = SkyCoord(radec, unit=(u.deg, u.deg))
        
        altitude = c.transform_to(frame)
        iic=np.mod(ic,7)
        lab=""
        plt.plot(delta_midnight,altitude.alt,label=lab,color="C"+str(iic),ls=lsarr[int(ic/7)])
        ic=ic+1
        
        plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                         sunaltazs.alt < -0*u.deg, color='0.5', zorder=0)
        plt.fill_between(delta_midnight.to('hr').value, 0, 90,
                         sunaltazs.alt < -18*u.deg, color='k', zorder=0)
        plt.ylim(0.,90)
        plt.xlim(-12, 12)
        plt.xticks(np.arange(13)*2 -12)
        plt.ylim(10, 90)
        plt.axhline(30.0,color="gray",lw=1)
        plt.xlabel('Hours from Midnight = '+midlocal.iso[0])
        plt.ylabel('Altitude [deg]')
        plt.title("")
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize=10)
    plt.show()
