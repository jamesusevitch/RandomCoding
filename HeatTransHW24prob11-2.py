# -*- coding: utf-8 -*-
"""
Created on Sun Mar  6 02:03:51 2016

@author: James Usevitch

Heat Transfer HW 24, Problem 11.2
"""

import numpy as np
import math
import matplotlib.pyplot as plt

#------------
# Part (a)
#------------


global Ai, Ao, Di, Do, Rfi, Rfo, \
       vw, va, Proa, Prob, Pri, ka, kw, hi

To = 15 + 273.15
Tmi = 75 + 273.15

umi = .5
Ai = math.pi*.022
Ao = math.pi*.027

Di = .022
Do = .027

Rfi = .0004
Rfo = .0002

hi = 72.1060136

vw = .389e-6
va = 1.483535e-5
Proa = .710081
Prob = 8.0225
Pri = 2.3492
ks = 15.1
kw = .66652
ka = .025352

# Re for part (c)
def Re(Vo, fluid='air'):
    global va, vw
    if fluid == 'air':
        return .027*Vo/va
    elif fluid == 'water':
        return .027*Vo/vw
    
def ho(Vo, fluid='air'):
    global Proa, Prob, ka, kw
    
    if fluid == 'air':
        return ka*.193*(Re(Vo, fluid))**(.618)*Proa**(1/3)
    elif fluid == 'water':
        # Change in a second        
        return kw*.027*(Re(Vo, fluid))**.805*Prob**(1/3)
    
def heatCoef(Vo, fluid='air'):
    global ks, Ao, Ai, hi, Rfi, Rfo, Do, Di
    
#    if fluid == 'air':
    return (1.0/Ao)*((1.0/(hi*Ai)) + Rfi/Ai + \
    math.log(Do/Di)/(2.0*math.pi*ks) + Rfo/Ao + \
    1.0/(ho(Vo, fluid)*Ao))**(-1)
#    elif fluid == 'water':

# Re for part (d)
def Re2(umi):
    global vw, Di
    return Di*umi/vw

def hi2(umi):
    global kw, Pri
    return kw*.023*(Re2(umi))**(4.0/5.0)*Pri**(.3)

def heatCoef2(umi, Vo=1):
    global ks, Ao, Ai, Rfi, Rfo, Do, Di
    
#    if fluid == 'air':
    return (1.0/Ao)*((1.0/(hi2(umi)*Ai)) + Rfi/Ai + \
    math.log(Do/Di)/(2.0*math.pi*ks) + Rfo/Ao + \
    1.0/(ho(Vo, 'water')*Ao))**(-1)
#    elif fluid == 'water':

# Calculate part (c) of the Homework
Vo = np.linspace(5, 30)
Coef = np.zeros(Vo.shape[0])

for i in range(Vo.shape[0]):
    Coef[i-1] = heatCoef(float(Vo[i-1]), 'air')

# Calculate part (d) of the Homework

umi = np.linspace(.5, 2.5)
Coef2 = np.zeros(umi.shape[0])

for i in range(umi.shape[0]):
    Coef2[i-1] = heatCoef2(float(umi[i-1]))

#print 'Coef2'
#print Coef2
#print 'umi'
#print umi
#print Coef
#print ' '
#print Vo

plt.figure(1, (12,7))

plt.title('Problem 11.2 Part (c)')
plt.plot(Vo, Coef, 'b')
plt.xlabel('Velocity (m/s)')
plt.ylabel('Overall Heat Transfer Coefficient (W/(m^2*K))')

plt.figure(2, (12,7))

plt.title('Problem 11.2 Part (d)')
plt.plot(umi, Coef2, 'r')
plt.xlabel('Mean Internal Velocity (m/s)')
plt.ylabel('Overall Heat Transfer Coefficient (W/(m^2*K))')