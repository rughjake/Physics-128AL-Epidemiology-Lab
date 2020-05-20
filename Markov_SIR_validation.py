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
t0 = 24 * 60 # step interval (currently 1 minute)
recovery = 14 * t0 # recovery time (currently 14 days)
R0 = 3 # ratio of infections to recoveries
N = 10000 # total population
I0 = 10 # amount of initially infected people
steps = 200000
#----------------------------------------------------

k = R0 / recovery # rate of infection
gamma = 1 / recovery # rate of recovery
S = N - I0 # susceptible
I = I0 # infected
R = 0 # recovered

S_arr = [S] # arrays of S,I,R,t values for plotting
I_arr = [I]
R_arr = [R]

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
    
plt.plot(delta_I_arr, "g", label="chance of infection") # testing purposes
plt.plot(delta_R_arr, "r", label="chance of recovery")
plt.plot(np.add(delta_I_arr, delta_R_arr), "y", label="sum of both")
plt.xlabel("time (minutes)")
plt.ylabel("probability")

plt.legend()
plt.show()