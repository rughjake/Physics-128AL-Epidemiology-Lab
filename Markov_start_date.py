# File : Markov_start_date.py
# Author : Jake Rugh
# Date : May 27, 2020
# Description : Explores the true start date of COVID-19 in New Zealand using
#   data from Markov_SIR.py

import matplotlib.pyplot as plt
import random
import numpy as np

final_I = 46 # people infected on March 1 (data from Markov_SIR.py)

# SIR/simulation parameters
#----------------------------------------------------
t0 = 24 * 60 # step interval (currently 60 steps per second)
serial_interval = 4 # recovery time (~4 given by various sources)
R0 = .311 * serial_interval # extracted from New Zealand's case data
N = 4886000 # total population (New Zealand's population)
I0 = 1 # amount of initially infected people
steps = 100000 # amount of steps to take in each simulation
sims = 1000 # amount of simulations to average over
#----------------------------------------------------

gamma = 1 / serial_interval / t0 # rate of recovery
k = R0 * gamma # rate of infection

end_times = []
end_times_hist = []

for a in range(0, sims):
    S = N - I0 # susceptible
    I = I0 # infected
    R = 0 # recovered
    
    for t in range(0, steps):
        delta_I = k * S * I / N # chance of infection
        delta_R = gamma * I # chance of recovery
    
        step = random.random() # picks random float between 0 and 1
    
        if step > 0 and step <= delta_I: # someone is infected
            I += 1
            S -= 1
        elif step > delta_I and step <= delta_R + delta_I: # someone recovers
            R += 1
            I -= 1
        else: # nothing happens
            pass
    
        if I + R == final_I:
            end_times.append(t / t0)
            end_times_hist.append(int(t / t0))
            break

plt.hist(end_times_hist)
plt.xlabel("time (days)")
plt.ylabel("occurrences")
plt.title("Simulated Date of First COVID-19 Case \n (given in days before Feb. 28, first confirmed case)")
plt.show()

print("mean : ", np.mean(end_times))
print("standard deviation : ", np.std(end_times))