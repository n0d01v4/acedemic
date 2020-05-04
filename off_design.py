# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:00:09 2019

@author: NIHARIKAPRIYADARSHI
"""
from array import *
import numpy as np
import math
import xlwt 
from xlwt import Workbook
pi=3.14159
Rad=pi/180
alpha1=72.47       #stator exit blade angle
A0=.0120176        #area at turbine inlet
A1=.0038008        #area upstream of stator exit
A3=.0039789        #area upstream of rotor inlet
A5=.0044543        #area downsream of rotor exit
BL1=0.95168        #blockage factor at stator exit
BL4=0.93718        #blockage factor at rotor exit
B4=-56.86          #rotor exit blade angle 
DELV=.05           #incremental value of stator exit critical velocity ratio
DELY=.05           #incremental valu of rotor exit critical velocity ratio
D2=12.9997         #stator exit diameter
D3=12.6238         #rotor inlet diameter
eff=0.913          #design value efficiency
G=1.6667           #gamma
G1=0.5*(G-1)
G2=G/(G-1)
G3=1/G2
G4=-(G+1)/(2*(G-1))
G5=G/(G+1)
G6=(G-1)/(G+1)
PD=0.98783         #stator total ressure ratio
P00=172368.9        #inlet total pressure
R=99.1976          #gas constant
T00=1144.44         #inlet total temerature
U3=237.9524        #rotor inlet tip speed
U4U3=.567765       #ratio of rotor exit mean speed to rotor inlet tip speed
VMAX=0.8           #final value of stator exit critical velocity ratio
V1=.30             #initial value of stator exit critical velocity ratio
w=0.338790         #design value of mass flow rate
XK=.28699          #rotor loss coefficient
YMAX=1.40          #final value of rotor exit critical velocity ratio
ZZ=22              #number of blades at rotor inlet
P0P5=1.74          #totol to total pressure ratio

#stator analysis
Vcr0=(2*G5*R*T00)**0.5
RHOVcr0=P00*(2*G5/(R*T00))**0.5
Vcr3=Vcr2=Vcr1=Vcr0
##
VVCR=0.5
while VVCR>0
    RVRVCR1=VVCR*(1-G6*VVCR**2)**(1/(G-1))
    RVRVCR0=PD*RVRVCR1*(A1/A0)*math.cos(alpha1.Rad)
    RVRVCR1i=w/(PD*RHOVcr0*A1*math.cos(alpha1.Rad))
    if RVRVCR1i==RVRVCR1:
        break
        es=e3ds
## equaion 49 is substituted into equation 50 which can then be solved iteratively for (v/vcr)0
w_cal=PD*RVRVCR1*RHOVcr0*A1*math.cos(alpha1*Rad)
#for value of v/vcr greator than 1.0 the choking value of mass flow rate is used to calculate a new stator exit flow angle
alpha1_=(math.acos(wcr/(PD*RVRVCR1*RHOVcr0*A1)))/Rad    #wcr not known here
#station 1
VuVcr1=VVcr1*math.sin(alpha1*Rad)
VrVcr1=VVcr1*math.cos(alpha1*Rad)
#station 2
RHOVcr2=RHOVcr1    #RHOVcr1 not defined
D2=D1
VuVcr2=VuVcr1
X=RVRVCR2*math.cos(alpha2*Rad)
Y=RVRVCR1*math.cos(alpha1*Rad)*(A1/A2)
Z=(1-G6*(VrVcr2**2+VuVcr2**2))**(1/(G-1))*VrVcr2
VVcr2=(VrVcr2**2+VuVcr2**2)**0.5
alpha2=(math.acos(VrVcr2/VVcr2))/Rad
#the above 5 equation is solved iteratively to determine vrvcr2,vvcr2,alpha2

#station 3
RHOVcr3=RHOVcr2
VuVcr3=VuVcr2*(D2/D3)
X=RVRVCR3*math.cos(alpha3*Rad)
Y=RVRVCR2*math.cos(alpha2*Rad)*(D2/D3)
Z=(1-G6*(VrVcr3**2+VuVcr3**2))**(1/(G-1))*VrVcr3
VVcr3=(VrVcr3**2+VuVcr3**2)**0.5
alpha3=(math.acos(VrVcr3/VVcr3))/Rad
#the above 5 equation is solved iteratively to determine vrvcr3,vvcr3,alpha3

#Rotor analysis
TR0T03=(1-G6*((2*U3*Vu3)/(Vcr3**2))-(U3/Vcr3)**2)
PR0P03=TR0T03**G2
WcrVcr3=TR0T03**0.5
RR0WcR0Vc3=TR0T03**(-G4)
# velocity diagram geometry
WuWcr3=(VuVcr3-U3/Vcr3)*(Vcr3/Wcr3)
WWcr3=(WuWcr3**2+VrVcr3**2)**0.5*(Vcr3/Wcr3)
WrWcr3=(WWcr3**2-WuWcr3**2)**0.5
Beta3=(math.asin(WuWcr3/WWcr3))/Rad
#optimum rotor inlet flow angle fi is calculated as follows
VuoptU3=1-(1.98/ZZ)
Vu3opt=U3*VuoptU3
Wu3opt=Vu3opt-U3
FI=(math.atan((Wu3opt/Vcr3)/VrVcr3))/Rad
I3=Beta3-FI                                   #rotor incidence angle
T4RT3R=1-G6*(U/Wcr3)**2*(1-(U4/U3)**2)
P4RidP3R=T4RT3R**(G2)
Wcr4Wcr3=T4RT3R**0.5
RW4RW3crid=T4RT3R**(-G4)
#rotor exit conditions
RAWx4=RAWr3
RWRRWcr4=(RWRRWcr3*(A3/A4)*(math.cos(Beta3*Rad)/math.cos(Beta4*Rad))/((PRPRid4)*RW4RW3crid)
A4=B4*A5
PRPRid4=(1-(G6*(K*WWcr4**2+(W3/Wcr4)**2*(K*(math.cos(I3*Rad))**2+(math.sin(I3*Rad))**n)))/(1-G6*WWcr4**2))**G2
#using above 3 equations w/wcr4 and pr/prid4 are determined by ieration process
#as the velocity raio w/wcr4 is increased beyond the choking value, the exit flow angle b4 is adjusted using below equation
Beta4=(math.acos((RWRRWcr3*(A3/A4)*math.cos(Beta3*Rad))/(PRPRid4*RW4RW3crid*RWRRWcr4)))
WuWcr4=WWcr4*math.sin(Beta4*Rad)
WxWcr4=WWcr4*math.cos(Beta4*Rad)
RRWcr5=RRWcr4
WuWcr5=WuWcr4
X=RWRRWcr5*math.cos(Beta5*Rad)
Y=RWRRWcr4*math.cos(Beta4*Rad)*B4
B4=A4/A5
Z=((1-(G6*(WxWcr5**2+WuWcr4**2)))**(1/(G-1)))*WxWcr5
WWcr5=(WxWcr5**2+WuWcr4**2)**0.5
Beta5=(math.sin(WuWcr5/WWcr5))/Rad
#the above 5 equations solved iteratively for wx/wcr5,w/wcr5 and b5
T0TR4=T0TR5=1-G6*WuWcr4**2+G6*VuWcr4**2
VuWcr4=WuWcr4+(U4/Vcr3)*(Vcr3/Wcr3)/Wcr4Wcr3
P0PR5=T0TR5**G2
VcrWcr5=T0TR5**0.5
VuWcr5=VuWcr4
WxWcr5=VxWcr4
VxVcr5=WxWcr5*WcrVcr5
VuVcr5=VuWcr4*WcrVcr5
VVcr5=(VxVcr5**2+VuVcr5**2)**0.5
alpha5=(math.asin(VuVcr5/VVcr5))
#overall turbine performance
T05T00=1-2*G6*((U3*Vu3/Vcr3**2)-(U4*Vu4/Vcr3**2))
Vu4Vcr3=WuWcr4*Wcr4Wcr3*WcrVcr3+(U4/Vcr3)
Vcr5=(T05T00**0.5)*Vcr0
P05P00=P01P00*PRP03*P4RidP3R*PRPRid4*P0PR5
P5P05=(1-G6*VVcr5**2)*G2
P5P00=P05P00*P5P05
nt=(1-T05T00)/(1-P05P00**G3)
ns=(1-T05T00)/(1-P5P00**G3)

#additional performance parameter
Thetacr=(G5*R*T00)/(G5*R*T0)*      #Need to see
delta=P00/P0*
epsilon=(0.7395945/G)*((G+1)/2)**G2
Neq=(C1*U3)/(D3*Thetacr**0.5)
DHeq=G2*R*T00*(1-T05T00)/Thetacr
Weq=w*epsilon*(Thetacr**0.5)*delta
Leq=(C2*w*DH0*epsilon)/(N*delta)
wNeq=w*N*epsilon/delta
DHids=(1-P5P00**G3)*G2*R*T00
DHidt=(1-P05P00**G3)*G2*R*T00
P00P5eq=(1-(G3/(R*T0*))*DHids/Thetacr)**(-G2)
P00P05eq=(1-(G3/(R*T0*))*DHidt/Thetacr)**(-G2)
P5P00eq=1/P00P5eq
DH0=DHeq*Thetacr
N=Neq*Thetacr**0.5
L=Leq*delta/epsilon
MU=U3/(2*DHids)**0.5
WF=DH0/U3**2
P=DH0*w/C3
WTP=w*(T00**0.5)/P00
DH0T00=DH0/T00
NT=N/T00**0.5
TORP=L/P00
T05=T00*(T05T00)
P05=P00*(P05P00)
RHO05=P05/(R*T05)
RHO5=RHO05*P5P05**(1/G)
Ns=(N*(w/RHO5)**0.5)/DHidt**0.75
