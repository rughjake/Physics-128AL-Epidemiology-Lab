# File : Markov_Q_comp.py
# Author : Jake Rugh
# Date : May 20, 2020
# Description : Compares different levels of quarantine measures and how they 
#   affect the rate of infection

import matplotlib.pyplot as plt
import random
import numpy as np

# SIR/simulation parameters
#----------------------------------------------------
t0 = 24 * 60
serial_interval = 4
R0 = .5 * serial_interval
N = 10000
I0 = 10 # amount of initially infected people
steps = 100000 # amount of steps to take in each simulation
sims = 20 # amount of simulations to average over
#----------------------------------------------------

# Quarantine parameters:
#----------------------------------------------------
Q_start = 5 # days after start when self-isolation starts
#----------------------------------------------------

gamma = 1 / serial_interval / t0 # rate of recovery
k = R0 * gamma # rate of infection

case_stacks = []

for b in range(0, 5):
    case_lists = []
    Q = b * .15
    
    for a in range(0, sims):
        S = N - I0 # susceptible
        I = I0 # infected
        R = 0 # recovered

        case_arr = [I]
        t_arr = [0]
    
        for t in range(t0, steps + t0):
            delta_I = k * S * I / N # chance of infection
            delta_R = gamma * I # chance of recovery
            if t >= Q_start * t0:
                delta_I *= (1 - Q)
    
            step = random.random() # picks random float between 0 and 1
    
            if step > 0 and step <= delta_I: # someone is infected
                I += 1
                S -= 1
            elif step > delta_I and step <= delta_R + delta_I: # someone recovers
                R += 1
                I -= 1
            else: # nothing happens
                pass
    
            case_arr.append(I + R) # arrays for plotting
            t_arr.append(t_arr[-1] + 1 / t0)

        case_lists.append(case_arr)
        
    case_stacks.append(case_lists)

plotting = []
for l in case_stacks:
    case_avg = []

    for a in range(0, steps + 1):
        case_vals = []
        
        for b in range(0, sims):
            case_vals.append(l[b][a])
        
        case_avg.append(np.average(case_vals))
    plotting.append(case_avg)
        
plt.plot(t_arr, plotting[0], "b", label="None")
plt.plot(t_arr, plotting[1], "g", label="Minimal")
plt.plot(t_arr, plotting[2], "r", label="Moderate")
plt.plot(t_arr, plotting[3], "k", label="Strict")

plt.xlabel("time (days)")
plt.ylabel("cases")
plt.title("Effects of Differing Levels of Quarantine")
plt.legend()
plt.show()