#!/usr/bin/python

import numpy as np
from scipy.integrate import odeint

# number of timesteps
nTimeSteps = 15

# starting redshift

z_start = 50
final_redshift = 0

print 'start redshift: ', z_start
print 'end redshift: ', final_redshift

# cosmological parameters
H0 = 70.5                       # Hubble constant in km/s/Mpc
Mpc_in_m = 3.08568025E22
H0_in_SI = H0/(Mpc_in_m/1000)   # Hubble constant in SI units
omega_m=0.2736                  # Omega Matter
omega_L=0.726                   # Omega Lambda
omega_r=8.5714E-5               # Omega radiation, actually irrelevant?

# helper function to get the index of an array where it has the value "value"
def find_nearest(array,value):
    idx=(np.abs(array-value)).argmin()
    return idx, array[idx]

# in order to get the scale factor evolution with respect to time
# we have to solve the differential equation on page 153 of Schneider,
# Einfuehrung in die Extragalaktische Astronomie und Kosmologie
def integrand(a,t):
    return np.sqrt( H0_in_SI**2*(omega_r/a**2 + omega_m/a+
                          (1-omega_m-omega_L)+omega_L*a**2))*31536000*1E9
a0 = 1E-11                      # initial condition at t=0
t = np.linspace(0,14,2000000)
# now solve the ode using odeint from scipy
a =  odeint(integrand,a0,t,rtol=1E-10)

# get the time at the end redshift
idx_end   = find_nearest(a,1.0/(final_redshift+1.0))[0]
timeEnd   = t[idx_end]
# get the time at the starting redshift
a_start   = 1./(z_start+1.0)
idx_start = find_nearest(a,a_start)[0]
timeStart = t[idx_start]

# specify the time vector for which we want to have the scale factor
t = np.linspace(timeStart,timeEnd,nTimeSteps)
a =  odeint(integrand,a_start,t,rtol=1E-10)

redshift = 1./a-1.0

# output a gadgettimerfile
np.savetxt('./gadget_timer_'+str(nTimeSteps)+'.txt',a,fmt='%10.7f')


# in the commented sections below we find an earlier version
# of this file where the differential equation was solved
# "by hand"

# def time2scaleFactor(lookback_time):

#     n = 300000                           # number of integration steps
#     lookback_time *=31536000*1E9         # *Gyr in seconds
#     h = float(lookback_time)/n           # t_incr
#     a=1.0
#     for i in range(0,n):
#         decr = np.sqrt( H0_in_SI**2*(omega_r/a**2 + omega_m/a+
#                                      (1-omega_m-omega_L)+omega_L*a**2))
#         a = a - decr*h
#     return a

# def scaleFactor2time(scaleFactor):
#     n = 300000
#     t = 0.
#     h = float(scaleFactor)/n
#     a = 0.00000000000000000001
#     for i in range(0,n):
#         t = t+h/np.sqrt(1./a**2*omega_r+1./a*omega_m+
#                         (1.-omega_m-omega_L)+a**2*omega_L)
#         a = a+h
#     t = t/H0_in_SI/31536000/1E9
#     return t

# make list of output times

# startTime = scaleFactor2time(1./(z_start+1.))
# endTime   = scaleFactor2time(1.0)
# timeSpan = endTime - startTime
# print startTime, endTime
# print 'The simulation covers a timespan of ', timeSpan, ' Gyrs'
# tstep = timeSpan/nTimeSteps

# # make array
# timeArray = np.zeros((nTimeSteps+1,4))
# # setting first value
# timeArray[0,0] = 0
# timeArray[0,1] = startTime
# timeArray[0,2] = 1./(z_start+1)
# timeArray[0,3] = z_start
# for i in range(1,nTimeSteps+1):
#     print 'Calculating scale factor for timestep ', i
#     timeArray[i,0] = i
#     timeArray[i,1] = startTime+i*tstep
#     timeArray[i,2] = time2scaleFactor(endTime-timeArray[i,1])
#     timeArray[i,3] = 1./timeArray[i,2]-1.

# print timeArray
