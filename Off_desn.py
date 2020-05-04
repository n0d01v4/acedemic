# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 23:17:51 2019

@author: NIHARIKAPRIYADARSHI
"""
from array import *
import numpy as np
import math
import xlwt 
from xlwt import Workbook

def func(x,G6,const,k): 
    return x*(1-G6*(x**2+k**2))**(1/(G-1))-const
def derivFunc( x,G6,k): 
    return (1-G6*(x**2+k**2))**(1/(G-1))-((2*G6*x**2)/(G-1))*(1-G6*(x**2+k**2))**((1/(G-1))-1) 
def newtonRaphson( x,G6,const,k ): 
    h = func(x,G6,const,k) / derivFunc(x,G6,k) 
    while abs(h) >= 0.0001: 
        h = func(x,G6,const,k)/derivFunc(x,G6,k) 
        x = x - h 
    print("The value of the root is : ", x)
    return x





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
RVRVCR1=w/(PD*RHOVcr0*A1*math.cos(alpha1*Rad))
RVRVCR0=PD*RVRVCR1*(A1/A0)*math.cos(alpha1*Rad)
x=2
x=abs(newtonRaphson( x,G6,RVRVCR1,0 ))
VVcr1=x

if VVCR>=1:
    wc=0.227
    alpha1=wc/(PD*RVRVCR1*RHOVcr0*A1)
VuVcr1=VVcr1*math.sin(alpha1*Rad)
VrVcr1=VVcr1*math.cos(alpha1*Rad)
P01=PD*P00
T01=T00
RHOVcr1=P01*(2*G5/(R*T01))**0.5
RHOVcr2=RHOVcr1
D1=D2
VuVcr2=VuVcr1
#X=RVRVCR2*math.cos(alpha2*Rad)
A2=A1
Y=RVRVCR1*math.cos(alpha1*Rad)*(A1/A2)
X=Y
x=2
VrVcr2=abs(newtonRaphson( x,G6,X,VuVcr2 ))
#Z=(1-G6*(VrVcr2**2+VuVcr2**2))**(1/(G-1))*VrVcr2
VVcr2=(VrVcr2**2+VuVcr2**2)**0.5
alpha2=(math.acos(VrVcr2/VVcr2))/Rad
print(Y)
X=RVRVCR2*math.cos(alpha2*Rad)
#the above 5 equation is solved iteratively to determine vrvcr2,vvcr2,alpha2
