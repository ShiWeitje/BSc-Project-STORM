#======-------------------------------------------------------------------
#
#     This function calculates the CAPE of a parcel with pressure PP (mb),
#       temperature TP (K) and mixing ratio RP (gm/gm) given a sounding
#       of temperature (T in K) and mixing ratio (R in gm/gm) as a function
#       of pressure (P in mb). CAPED is
#       the calculated value of CAPE and TOB is the temperature at the
#       level of neutral buoyancy.  IFLAG is a flag
#       integer. If IFLAG = 1, routine is successful; if it is 0, routine did
#       not run owing to improper sounding (e.g.no water vapor at parcel level).
#       IFLAG=2 indicates that routine did not converge.
#
#======-------------------------------------------------------------------
#
import numpy as np
#
def cape(TP,RP,PP,T,R,P,SIG) :
   #
   #====== Change to a float type array
   #
#   TP = TP.astype('float')
#   PP = PP.astype('float')
#   RP = RP.astype('float')
#   T = T.astype('float')
#   P = P.astype('float')
#   R = R.astype('float')
   #
   #====== Pressure below which sounding is ignored
   ptop=59
   #
   Nold=np.max(np.size(P))
   N=0
   for i in np.arange(Nold-1,-1,-1) :
      if (P[i] > ptop) :
         N=max(N,i)
         break
         #
   if (N+1 < Nold) :
      P=P[:N+1]
      T=T[:N+1]
      R=R[:N+1]
      #
   TVRDIF=np.zeros(N+1)
   #
   #====== Get minimum sounding level at or above PP 
   #
   #print("N: ", N) # swy 2/3/2022
   for i in range(0,N+1) :
      #print("JMIN %i ?" % i)   #swy
      #print("P[i] = ", P[i])   #swy
      if (P[i] <= PP) :
         JMIN=i
         #print("JMIN = %i succesful, now break" % i)
         break
         #
   #====== Default values  
   #
   CAPED=0.0
   TOB=T[0]
   IFLAG=1
   #
   #====== Check that sounding is suitable
   #
   if (RP < 1e-6) | (TP < 200) :
      IFLAG=0
      CAPED = np.nan
      TOB = np.nan
      return CAPED,TOB,IFLAG
      #
   #====== Assign values of thermodynamic constants
   #
   CPD=1005.7
   CPV=1870.0
   CL=4190.0
   CL=2500.0
   CPVMCL=CPV-CL
   RV=461.5
   RD=287.04
   EPS=RD/RV
   ALV0=2.501e6
   #
   #====== Define various parcel quantities, including reversible 
   #======                      entropy, S.                       
   #
   TPC=TP-273.15
   ESP=6.112*np.exp(17.67*TPC/(243.5+TPC))
   EVP=RP*PP/(EPS+RP)
   RH=EVP/ESP
   RH=min(RH,1.0)
   ALV=ALV0+CPVMCL*TPC
   S=(CPD+RP*CL)*np.log(TP)-RD*np.log(PP-EVP)+ALV*RP/TP-RP*RV*np.log(RH)
   #
   #====== Find lifted condensation pressure, PLCL 
   #
   CHI=TP/(1669.0-122.0*RH-TP)
   PLCL=PP*(RH**CHI)
   #
   #====== Begin updraft loop 
   #
   NCMAX=0
   #
   for J in range(JMIN,N+1) :
      #
      #====== Parcel quantities below lifted condensation level 
      #
      if (P[J] >= PLCL) :
         TG=TP*(P[J]/PP)**(RD/CPD)
         RG=RP
         # 
         #====== Calculate buoyancy 
         # 
         TLVR=TG*(1.+RG/EPS)/(1.+RG)
         TVRDIF[J]=TLVR-T[J]*(1.+R[J]/EPS)/(1+R[J])
      else :
         # 
         #====== Parcel quantities above lifted condensation level  ***
         # 
         TGNEW=T[J]
         TJC=T[J]-273.15
         ES=6.112*np.exp(17.67*TJC/(243.5+TJC))
         RG=EPS*ES/(P[J]-ES)
         #
         #====== Iteratively calculate lifted parcel temperature and mixing 
         #======                ratio for reversible ascent                  
         #
         NC=0
         TG=0.0
         #
         while ((np.abs(TGNEW-TG)) > 0.001) :
            #
            TG=TGNEW
            TC=TG-273.15
            ENEW=6.112*np.exp(17.67*TC/(243.5+TC))
            RG=EPS*ENEW/(P[J]-ENEW)
            #
            NC=NC+1
            #
            #====== Calculate estimates of the rates of change of the entropy  
            #====== with temperature at constant pressure    
            #
            ALV=ALV0+CPVMCL*(TG-273.15)
            SL=(CPD+RP*CL+ALV*ALV*RG/(RV*TG*TG))/TG
            EM=RG*P[J]/(EPS+RG)
            SG=(CPD+RP*CL)*np.log(TG)-RD*np.log(P[J]-EM)+ALV*RG/TG
            #
            if (NC < 3) :
               AP=0.3
            else :
               AP=1.0
               #
            TGNEW=TG+AP*(S-SG)/SL;
            #
            #------ Bail out if things get out of hand 
            #
            if (NC > 500) | (ENEW > (P[J]-1)) :
               IFLAG=2
               #
         NCMAX = max(NC,NCMAX)
         #
         #------ Calculate buoyancy 
         #
         RMEAN=SIG*RG+(1-SIG)*RP
         TLVR=TG*(1.+RG/EPS)/(1.+RMEAN)
         TVRDIF[J]=TLVR-T[J]*(1.+R[J]/EPS)/(1.+R[J])
         #
   #====== Begin loop to find NA, PA, and CAPE from reversible ascent
   #
   NA=0.0
   PA=0.0
   #
   #====== Find maximum level of positive buoyancy, INB 
   #
   INB=1
   for J in np.arange(N,JMIN-1,-1) :
      if TVRDIF[J] > 0 :
         INB=max(INB,J)
         #
   if (INB == JMIN) :
      return CAPED,TOB,IFLAG
      #
   #====== Find positive and negative areas and CAPE 
   #
   if (INB > 1) :
      for J in range((JMIN+1),INB+1) :
         PFAC=RD*(TVRDIF[J]+TVRDIF[J-1])*(P[J-1]-P[J])/(P[J]+P[J-1])
         PA=PA+max(PFAC,0.0)
         NA=NA-min(PFAC,0.0)
         #
      #====== Find area between parcel pressure and first level above it
      #
      PMA=(PP+P[JMIN]) 
      PFAC=RD*(PP-P[JMIN])/PMA
      PA=PA+PFAC*max(TVRDIF[JMIN],0.0)
      NA=NA-PFAC*min(TVRDIF[JMIN],0.0)
      #
      #====== Find residual positive area above INB and TO 
      #
      PAT=0.0
      TOB=T[INB]
      if (INB < N) :
         PINB=(P[INB+1]*TVRDIF[INB]-P[INB]*TVRDIF[INB+1])/(TVRDIF[INB]-TVRDIF[INB+1])
         PAT=RD*TVRDIF[INB]*(P[INB]-PINB)/(P[INB]+PINB)
         TOB=(T[INB]*(PINB-P[INB+1])+T[INB+1]*(P[INB]-PINB))/(P[INB]-P[INB+1])
         #
      #====== Find CAPE
      #
      CAPED=PA+PAT-NA
      CAPED=max(CAPED,0.0)
      #
   return CAPED,TOB,IFLAG
   #
#TP  = 26.0+273.15
#RP  = 17.6*1e-3
#PP  = 1000
#T   = np.asarray([26.0,23.0,19.8,17.3,14.6,11.8, 8.6, 5.1, 1.4,-2.5,-6.9,-11.9,-17.7,-24.8,-33.2,-43.3,-55.2,-61.5,-67.6])+273.15
#R   = np.asarray([17.6,15.3,13.0,11.0, 8.4, 7.1, 5.8, 4.6, 3.6, 3.2, 2.1,  1.4,    0,    0,    0,    0,    0,    0,    0] )*1e-3
#P   = np.asarray([1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500,  450,  400,  350,  300,  250,  200,  150,  100] )
#SIG = 0.0
#print( cape(TP,RP,PP,T,R,P,SIG) )











