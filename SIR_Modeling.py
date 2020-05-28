# -*- coding: utf-8 -*-
"""
Phys 128L: Epidemiology Lab
SIR Model
Written by Hunter Lukk

WARNING: Code is quite sloppy, and code changes are required to access different graphs/simulations. It ain't pretty but it works

The basic system obeys the following set of equations:
    S' = -k * SI, where k is beta/N
    I' = k * SI - gamma * I
    R' = gamma * I
    
With the additional restriction that S + I + R = N
    
"""

import matplotlib.pyplot as plt
import csv
import numpy as np
import scipy.stats as sp


class SIR:
    def __init__(self, ICs, params):
        self.ICs = ICs
        self.beta = float(params[0])
        self.gamma = float(params[1])
        self.drate = float(params[2])
        self.q1 = float(params[3])
        self.q2 = float(params[4])
        self.q3 = float(params[5])
        self.q4 = float(params[6])
        self.t0 = ICs[4]
        N = 0.0
        
        for i in range(4):
            N += float(ICs[i])
            
        self.N = N
        self.k = self.beta/N
    
    def advance(self, vals, day, dt=1):
        #takes in the current SIR values and linearizes to produce a change in the values for a time step dt
        #returns the incremented SIR values
        
        #This if statement checks to see if the quarantine condition has been reached
        #It makes sure that the quarantine flag qd has not been set to True yet
        if day >= self.q1 and day < self.q2:
            newS = vals[0] + dt *  (-0.8 * self.k * vals[0] * vals[1])
            newI = vals[1] + dt * (0.8 * self.k * vals[0] - self.gamma - self.drate) * vals[1]
            newR = vals[2] + dt * (self.gamma * vals[1])
            newD = vals[3] + dt * (self.drate * vals[1])
            
        elif day >= self.q2 and day < self.q3:
            newS = vals[0] + dt *  (-0.7 * self.k * vals[0] * vals[1])
            newI = vals[1] + dt * (0.7 * self.k * vals[0] - self.gamma - self.drate) * vals[1]
            newR = vals[2] + dt * (self.gamma * vals[1])
            newD = vals[3] + dt * (self.drate * vals[1])
            
        elif day >= self.q3 and day < self.q4:
            newS = vals[0] + dt *  (-0.6 * self.k * vals[0] * vals[1])
            newI = vals[1] + dt * (0.6 * self.k * vals[0] - self.gamma - self.drate) * vals[1]
            newR = vals[2] + dt * (self.gamma * vals[1])
            newD = vals[3] + dt * (self.drate * vals[1])
            
        elif day >= self.q4:
            newS = vals[0] + dt *  (-0.5 * self.k * vals[0] * vals[1])
            newI = vals[1] + dt * (0.5 * self.k * vals[0] - self.gamma - self.drate) * vals[1]
            newR = vals[2] + dt * (self.gamma * vals[1])
            newD = vals[3] + dt * (self.drate * vals[1])
            
        else:
            newS = vals[0] + dt *  (-self.k * vals[0] * vals[1])
            newI = vals[1] + dt * (self.k * vals[0] - self.gamma - self.drate) * vals[1]
            newR = vals[2] + dt * (self.gamma * vals[1])
            newD = vals[3] + dt * (self.drate * vals[1])
        
        
        newt = vals[4] + dt
        return [newS, newI, newR, newD, newt]
        
    def getSets(self, nsteps, dt=1):
        #this function gets the sets of values for SIR over a specified time period (nsteps * dt)
        #it does this by calling advance several times and adding the returns to a list
        #it then returns the 3x(nsteps+1) array of values for the functions
        
        Sset = [0 for j in range(nsteps + 1)]
        Iset = [0 for j in range(nsteps + 1)]
        Rset = [0 for j in range(nsteps + 1)]
        Dset = [0 for j in range(nsteps + 1)]
        tset = [0 for j in range(nsteps + 1)]
        Sset[0], Iset[0], Rset[0], Dset[0], tset[0] = self.ICs
        
        sets = [Sset, Iset, Rset, Dset, tset]
        for i in range(nsteps):
            sets[0][i + 1], sets[1][i + 1], sets[2][i + 1], sets[3][i + 1], sets[4][i + 1] = self.advance([sets[0][i], sets[1][i], sets[2][i], sets[3][i], sets[4][i]], i*dt, dt)
        return sets
    
    def graph(self, nsteps, dt=1, total = False, raw = False, log = False):
        sets = self.getSets(nsteps, dt)
        
        if total:
            for i in range(len(sets[0])):
                sets[1][i] = self.N - sets[0][i]
                
            fig = plt.figure()
            ax1 = fig.add_subplot(111)
            if not raw:
                ax1.plot(sets[4], sets[1])
                print('No Q Total:',sets[1][len(sets[1])-1])
                plt.title('No Quarantine: Total Cases')
                plt.axhline(1200, color = 'r', linestyle = '--', alpha = 0.7, label = 'Leveling of Cases with Quarantine')
                plt.xlabel('Days')
                plt.ylabel('Count')
                plt.legend()
                plt.savefig('NoQuarantine.pdf')
                plt.show()
                return
            plt.title('Cases with Quarantine Delayed by 1 Week')
        
        if raw:
                data = []
                day = []
                x = 49
                
                with open('NewZealand.csv','r') as file:
                    reader = csv.reader(file)
                    for line in reader:
                        data.append(float(line[2]))
                        day.append(float(line[3]))
                        
                for i in range(len(day)):
                    day[i] = day[i] - x
                        
                print('SIR Model Final count:',sets[1][len(sets[1])-1])
                        
                ax1.plot(sets[4], sets[1], color = 'g', label = 'SIR Model w/ delayed Q.')
                ax1.scatter(day[x:], data[x:], marker = '.', label = 'Raw Data')
                plt.legend()
            
        else:
            plt.plot(sets[4], sets[0], label = 'S')
            plt.plot(sets[4], sets[1], label = 'I')
            plt.plot(sets[4], sets[2], label = 'R')
            if self.drate != 0:
                plt.plot(sets[4], sets[3], label = 'D')
        
        plt.ylabel('Count')
        plt.xlabel('Days')
        plt.legend()
        if log:
            plt.ylabel('Count (logarithmic)')
            plt.yscale('log')
        plt.savefig("NewZealandDelay.pdf")
        plt.show()
    
        
if __name__ == "__main__":
    RealIC = [4886000, 15, 0, 0, 0]
    RealParam = [0.308, 1.0/4, 0, 41, 43, 45, 47]
    
    noQParam = [0.308, 1.0/4, 0, 999, 999, 999, 999]
    '''
    data = []
    day = []
    
    with open('NewZealand.csv','r') as file:
        reader = csv.reader(file)
        for line in reader:
            data.append(float(line[2]))
            day.append(float(line[3]))
    
    
    for i in range(len(data)):
        if data[i]>0:
            data[i] = np.log(data[i])
            
    
    
    m,b = np.polyfit(day[70:80],data[70:80],1)
    linx = np.linspace(0,100,100)
    liny = m * linx + b
    #print(m)
    #print(b)
    
    for i in range(100):
        if (np.abs(1 - (m * linx[i] + b))) <= 0.1:
            print(i)
    '''
    '''
    m, b, r, p, err = sp.linregress(day[70:80], data[70:80])
    linx = np.linspace(65, 90, 40)
    liny = m * linx + b
    print('Beta = ',m)
    print('R-squared = ',r*r)
    '''
    '''
    plt.scatter(day[49:], data[49:], marker = '.')
    #plt.plot(linx, liny, color = 'r', linestyle = '--', alpha = 0.5, label = 'Linear Regression: beta = 0.308, R-squared = 0.996')
    plt.title('New Zealand Case Count')
    plt.ylabel('Count')
    plt.xlabel('Days')
    #plt.legend()
    plt.savefig('NewZealand.pdf')
    plt.show()
    '''
    
    sir = SIR(RealIC, [0.308, 0.25, 0, 48, 50, 52, 54])
    sir.graph(9500, 0.01, True, True)
    sets = sir.getSets(9500, 0.01)
    print('One Week Delay of Quarantine gives', 100 * ((1829.195 / 1154) - 1),'% increase in total cases')
    #noQ = SIR(RealIC, noQParam)
    #noQ.graph(30000, 0.01, True)
    
    #print('Quarantine case reduction: approx ', 100 * (1 - (1154.0 / 1750000)),'%')