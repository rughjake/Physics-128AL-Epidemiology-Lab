# File : R0_extraction.py
# Author : Jake Rugh
# Date : May 20, 2020
# Description : Extract R0 from country data by finding slope of logarithmic
#   plot of confirmed cases. We have chosen New Zealand as our country.
# Note : I am studying only the linear part portion of the logarithmic data from
#   New Zealand.

import numpy as np
import matplotlib.pyplot as plt


file = open("total-cases-covid-19.csv", "r") # file containing NZ case data
entries = file.readlines()
file.close()
data = []
for i in entries:
    if i[0:11] == "New Zealand":
        data.append(i)

n = 0        
while data[n][17:20] != "Mar":
    n += 1
n += 6

dates = []
log_cases = []
for m in range(n, n + 13):
    line = data[m]
    line = line.split(",")
    dates.append(int(line[2][5:]))
    log_cases.append(np.log(int(line[-1])))
    
sum_x = 0
sum_x_squared = 0
sum_y = 0
sum_xy = 0
length = len(dates)
for n in range(0, length):
    sum_x += dates[n]
    sum_x_squared += dates[n]**2
    sum_y += log_cases[n]
    sum_xy += log_cases[n] * dates[n]

delta = length * sum_x_squared - sum_x**2

A = (sum_x_squared * sum_y - sum_x * sum_xy) / delta
B = (length * sum_xy - sum_x * sum_y) / delta
   
x = np.linspace(1, 29, 1000)
y = A + B*x

sigma_y = 0
for n in range(0, length):
    sigma_y += (log_cases[n] - A - B * dates[n])**2

sigma_y /= (length - 2)
sigma_y = np.sqrt(sigma_y)

sigma_A = sigma_y * np.sqrt(sum_x_squared / delta)
sigma_B = sigma_y * np.sqrt(length / delta)

print("A = ", A)
print("B = ", B)
print("sigma_A = ", sigma_A)
print("sigma_B = ", sigma_B)

plt.plot(dates, log_cases, "o")
plt.plot(x, y, 'r')
plt.xlabel("days")
plt.ylabel("log(cases)")


































