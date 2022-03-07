#   Maximum potential intensity, changed to python according to the code
#   from Emanuel's website
# 
#  INPUT:   SST: Sea surface temperature in C
#
#           PSL: Sea level pressure (mb)
#
#           P,T,R: One-dimensional arrays
#             containing pressure (mb), temperature (C),
#             and mixing ratio (g/kg). The arrays MUST be
#             arranged so that the lowest index corresponds
#             to the lowest model level, with increasing index
#             corresponding to decreasing pressure. The temperature
#             sounding should extend to at least the tropopause and
#             preferably to the lower stratosphere, however the
#             mixing ratios are not important above the boundary
#             layer. Missing mixing ratios can be replaced by zeros.
#
#
#  OUTPUT:  PMIN is the minimum central pressure, in mb
#
#           VMAX is the maximum surface wind speed, in m/s
#                  (reduced to reflect surface drag)
#
#           TO is the outflow temperature (K)
#
#           IFL is a flag: A value of 1 means OK; a value of 0
#              indicates no convergence (hypercane); a value of 2
#              means that the CAPE routine failed.
#
#-----------------------------------------------------------------------------
#
import numpy as np
from cape import cape
#import cape #swy 1/3/2022
#
def pcmin(SST,PSL,P,T,R):
   #
   #   ***   Adjustable constant: Ratio of C_k to C_D    ***
   #
   CKCD=0.9
   #
   #   ***   Adjustable constant for buoyancy of displaced parcels:  ***
   #   ***    0=Reversible ascent;  1=Pseudo-adiabatic ascent        ***
   #
   SIG=0.0
   #
   #   ***  Adjustable switch: if IDISS = 0, no dissipative heating is   ***
   #   ***     allowed; otherwise, it is                                 ***
   #
   IDISS=1
   #
   #   ***  Exponent, b, in assumed profile of azimuthal velocity in eye,   ***
   #   ***   V=V_m(r/r_m)^b. Used only in calculation of central pressure   ***
   #
   b=2.0
   #
   #   *** Set level from which parcels lifted   ***
   #
   NK=0
   #
   #   *** Factor to reduce gradient wind to 10 m wind
   #
   VREDUC=0.8
   #
   #
   #--------------------------------------------------------------------------
   #
   RAT = np.nan
   CAPEMS = np.nan
   CAPEM = np.nan
   FAC = np.nan
   #
   SSTK=SST+273.15
   TOMS=230.0
   ES0=6.112*np.exp(17.67*SST/(243.5+SST))
   R=R*0.001
   T=T+273.15
   #
   
   if (SST <= 5.0) | (np.min(T) <= 100.0) :
      VMAX = np.nan
      PMIN = np.nan
      TO   = np.nan
      IFL  = 0
      #print("Top nan")
      return PMIN,VMAX,TO,IFL
      #
   IFL=1
   NP=0
   PM=970.0
   PMOLD=PM
   PNEW=0.0
   #
   #   ***   Find environmental CAPE ***
   #
   TP=T[NK]
   RP=R[NK]
   PP=P[NK]
   CAPEA, tmp, IFLAG= cape(TP,RP,PP,T,R,P,SIG)
   if (IFLAG != 1) :
      IFL=2
      #
   #   ***   Begin iteration to find mimimum pressure   ***
   #
   rec = 0
   while (abs(PNEW-PMOLD)) > 0.2 :
      #
      #   ***  Find CAPE at radius of maximum winds   ***
      #
      TP=T[NK]
      PP=min(PM,1000.0)
      RP=0.622*R[NK]*PSL/(PP*(0.622+R[NK])-R[NK]*PSL)
      #print("pcmin_fullterm: P = ", P) # swy 2/3/2022
      CAPEM, TOM, IFLAG=cape(TP,RP,PP,T,R,P,SIG)
      if (IFLAG != 1) :
         IFL=2
         #
      #  ***  Find saturation CAPE at radius of maximum winds   ***
      #
      TP=SSTK
      PP=min(PM,1000.0)
      RP=0.622*ES0/(PP-ES0)
      CAPEMS, TOMS, IFLAG=cape(TP,RP,PP,T,R,P,SIG)
      TO=TOMS
      if (IFLAG != 1) :
         IFL=2
         #
      RAT=SSTK/TOMS
      if (IDISS == 0) :
         RAT=1.0
         #
      #  ***  Initial estimate of minimum pressure   ***
      #
      RS0=RP
      TV1=T[0]*(1.+R[0]/0.622)/(1.+R[0])
      TVAV=0.5*(TV1+SSTK*(1.+RS0/0.622)/(1.+RS0))
      CAT=CAPEM-CAPEA+0.5*CKCD*RAT*(CAPEMS-CAPEM)
      CAT=max(CAT,0.0)
      PNEW=PSL*np.exp(-CAT/(287.04*TVAV))
      #
      #   ***  Test for convergence   ***
      #
      PMOLD=PM
      PM=PNEW
      NP=NP+1
      
      if (NP > 200) | (PM < 400) :
         PMIN=np.nan
         VMAX=np.nan
         TO=np.nan
         IFL=0
         #print("Bottom nan")
         return PMIN,VMAX,TO,IFL
         # 
   CATFAC=0.5*(1.+1./b)
   CAT=CAPEM-CAPEA+CKCD*RAT*CATFAC*(CAPEMS-CAPEM)
   CAT=max(CAT,0.0)
   PMIN=PSL*np.exp(-CAT/(287.04*TVAV))
   
   # if np.isnan(PMIN) == True:       # swy 4/3/2022
   #     print("unknown nan???")
   
    
   FAC=max(0.0,(CAPEMS-CAPEM))
   VMAX=VREDUC*np.sqrt(CKCD*RAT*FAC)
         #
   return(PMIN,VMAX,TO,IFL,RAT,CAPEMS,CAPEM,FAC,CAPEA)
   #
