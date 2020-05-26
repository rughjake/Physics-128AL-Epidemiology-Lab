# File : Markov_SIR.py
# Author : Jake Rugh
# Date : May 20, 2020
# Description : Models S, I, R counts over time for certain starting conditions
#   Set parameters such as interval time, recovery time, R0, total population, 
#   initially infected population, number of intervals in each simulation, and
#   number of simulations to run. Takes average and standard deviation over all
#   simulations.
# Note : step interval must be short enough so that only one transition can
#   happen, i.e. the sum of chances of a recovery or infection must never
#   exceed 100% for any step

import matplotlib.pyplot as plt
import random
import numpy as np

# SIR/simulation parameters
#----------------------------------------------------
t0 = 24 * 60 # step interval (currently 60 steps per second)
recovery = 4 * t0 # recovery time (~4 given by various sources)
R0 = .308 # extracted from New Zealand's case data
N = 4886000 # total population (New Zealand's population)
I0 = 20 # amount of initially infected people
steps = 100000 # amount of steps to take in each simulation
sims = 20 # amount of simulations to average over
#----------------------------------------------------

# Self-isolation parameters:
#----------------------------------------------------
iso_start = 21 # days after start when self-isolation starts
self_isolation = .10 # infection reduction after self-isolation
partial_start = 23 # start of partial lockdown
partial_lockdown = .15 # infection reduction after partial lockdown
full_start = 25 # start of full lockdown
full_lockdown = .40 # infection reduction after total lockdown
#----------------------------------------------------

k = R0 / recovery # rate of infection
gamma = 1 / recovery # rate of recovery

S_lists = []
I_lists = []
R_lists = []

for a in range(0, sims):
    S = N - I0 # susceptible
    I = I0 # infected
    R = 0 # recovered

    S_arr = [S] # arrays of S,I,R,t values for plotting
    I_arr = [I]
    R_arr = [R]
    t_arr = [0]
    
    for t in range(0, steps):
        delta_I = k * S * I / N # chance of infection
        delta_R = gamma * I # chance of recovery
        if t >= iso_start * t0 and t < partial_start * t0:
            delta_I *= (1 - self_isolation)
        elif t >= partial_start * t0 and t < full_start * t0:
            delta_I *= (1 - partial_lockdown)
        elif t >= full_start * t0:
            delta_I *= (1 - full_lockdown)
    
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

    S_lists.append(S_arr)
    I_lists.append(I_arr)
    R_lists.append(R_arr)


S_avg = [] # averages
I_avg = []
R_avg = []

S_dev = [] # standard deviations
I_dev = []
R_dev = []

for a in range(0, steps + 1):
    S_vals = []
    I_vals = []
    R_vals = []
    
    for b in range(0, sims):
        S_vals.append(S_lists[b][a])
        I_vals.append(I_lists[b][a])
        R_vals.append(R_lists[b][a])
        
    S_avg.append(np.average(S_vals))
    I_avg.append(np.average(I_vals))
    R_avg.append(np.average(R_vals))
    
    S_dev.append(np.std(S_vals))
    I_dev.append(np.std(I_vals))
    R_dev.append(np.std(R_vals))
    
S_max = np.add(S_avg, S_dev)
I_max = np.add(I_avg, I_dev)        
R_max = np.add(R_avg, R_dev)

S_min = np.subtract(S_avg, S_dev)
I_min = np.subtract(I_avg, I_dev)        
R_min = np.subtract(R_avg, R_dev)
        
#plt.plot(t_arr, S_avg, "b", label="susceptible")
#plt.plot(t_arr, I_avg, "g", label="infected")
#plt.plot(t_arr, R_avg, "r", label="recovered")
plt.plot(t_arr, np.add(I_avg, R_avg), "y", label="confirmed cases")

#plt.fill_between(t_arr, S_max, S_min, facecolor="blue", alpha=0.5)
#plt.fill_between(t_arr, I_max, I_min, facecolor="green", alpha=0.5)
#plt.fill_between(t_arr, R_max, R_min, facecolor="red", alpha=0.5)
plt.fill_between(t_arr, np.add(I_max, R_max), np.add(I_min, R_min), facecolor="yellow", alpha=0.5)

plt.xlabel("time (days)")
plt.ylabel("cases")
plt.legend()
plt.show()