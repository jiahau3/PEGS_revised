import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def linear(x, y):
    return np.polyfit(x, y, 1)

def plotting(x1, y1, x2, y2, fitted_formula, Nameoutput):  
    plt.figure(figsize=(5,5))
#    axes = fig.add_axes([0.1,0.1,0.8,0.8])
    plt.plot(x1, y1, 'o')
    plt.plot(x2, y2, 'r-')
    plt.xlabel('$N_{fringe}$', fontsize='14')
    plt.ylabel('$Force(N)$', fontsize='14')
    plt.legend(['original data','$1^{st}$ order fitting'], loc='best')
    plt.text(1, 0, fitted_formula)
    plt.savefig('..\\force_measurement\\'+ Nameoutput)
    plt.show()

def loadfile(filename):      
    return np.array(pd.read_excel(filename))

def calculate_Fsigma(slope, radius):    
    return 4 * slope / (np.pi * radius)

data = loadfile('..\\force_measurement\\opticstress_coefficient_calculation.xlsx')
radius = 4 * 10**(-3)
g_to_newton = 0.00981
N_fringe = data[:-3,11]
force = data[:-3,12] * g_to_newton
name_output = str(int(radius*2*10**(3))) + 'mm_circular_5.png'

p1 = linear(N_fringe, force)
x_fit = np.linspace(N_fringe[0], N_fringe[-1], 10)
y_fit = np.poly1d(p1)
F_sigma = calculate_Fsigma(p1[0], radius)
formula = 'y = ' + str(round(p1[0], 4)) + 'x' + str(round(p1[1], 4)) + ', ' + \
'$F_{\sigma}$ = ' + str(F_sigma)


plotting(N_fringe, force, x_fit, y_fit(x_fit), formula, name_output)

print('slope = ', p1[0])
print('$F_{\sigma}$:', F_sigma)

