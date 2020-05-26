# File : Markov_SIR_validation.py
# Author : Jake Rugh
# Date : May 18, 2020
# Description : Use to confirm appropriate interval size.
# Note : step interval must be short enough so that only one transition can
#   happen, i.e. the sum of chances of a recovery or infection must never
#   exceed 100% for any step

import matplotlib.pyplot as plt
import random
import numpy as np

# Only edit these parameters:
#----------------------------------------------------
t0 = 24 * 60 * 60 * 60 # step interval (currently 60 steps per second)
recovery = 4 * t0 # recovery time (taken from various sources)
R0 = 3.226 # ratio of infections to recoveries
N = 4886000 # total population
I0 = 10 # amount of initially infected people
steps = 150000 * 60 * 30
#----------------------------------------------------

k = R0 / recovery # rate of infection
gamma = 1 / recovery # rate of recovery
S = N - I0 # susceptible
I = I0 # infected
R = 0 # recovered

S_arr = [S] # arrays of S,I,R,t values for plotting
I_arr = [I]
R_arr = [R]

t_arr = [0]
delta_I_arr = [] # testing purposes
delta_R_arr = []

for n in range(0, steps):
    delta_I = k * S * I / N # chance of infection
    delta_R = gamma * I # chance of recovery
    
    delta_I_arr.append(delta_I) # testing purposes
    delta_R_arr.append(delta_R)
    
    step = random.random() # picks random float between 0 and 1
    
    if step > 0 and step <= delta_I: # someone is infected
        I += 1
        S -= 1
    elif step > delta_I and step <= delta_R + delta_I: # someone recovers
        R += 1
        I -= 1
    else: # nothing happens
        pass
    
    S_arr.append(S) # arrays for plotting
    I_arr.append(I)
    R_arr.append(R)
    t_arr.append(t_arr[-1] + 1 / t0)
    
del t_arr[-1]
    
plt.plot(t_arr, delta_I_arr, "g", label="chance of infection") # testing purposes
plt.plot(t_arr, delta_R_arr, "r", label="chance of recovery")
plt.plot(t_arr, np.add(delta_I_arr, delta_R_arr), "y", label="sum of both")
plt.xlabel("time (days)")
plt.ylabel("probability")

plt.legend()
plt.show()