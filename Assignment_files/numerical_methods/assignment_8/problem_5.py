'''
Python implementation of solution to problem #5.
'''

import numpy as np
import scipy.integrate as scp
import matplotlib.pyplot as plt


def odeCRK4(F,a,b,N,yINI):
    stages = 4
    c = np.array([0., 0.5, 0.5, 1.0],dtype=np.float64)
    B = np.array([1./6., 1./3.,1./3.,1./6.],dtype=np.float64)
    A = np.array([[0, 0, 0, 0],
                  [0.5,0,0,0],
                  [0,0.5,0,0],
                  [0,0,1.0,0]],
                 dtype=np.float64)
    x = np.linspace(a, b,N,dtype=np.float64)
    #sys_size = 1 # for now, just solve scalar systems
    y = np.zeros(N,dtype=np.float64)
    y[0] = yINI
    h = x[1]-x[0];
    for t in range(N-1):
        Xi = np.zeros(stages,dtype=np.float64)
        for s in range(stages):
            Xi[s] = y[t]
            for i in range(s-1):
                Xi[s] = Xi[s]+h*A[s,i]*F(x[t]+c[i]*h,Xi[i])
                
                
        y[t+1]=y[t]
        for i in range(stages):
            y[t+1]=y[t+1]+h*B[i]*F(x[t]+c[i]*h,Xi[i])
    return x,y


# Parameters
A_tank = 3.13 # m^2, cross sectional area of the tank
A_pipe = 0.06 # m^2, cross sectional area of drain pipe
C = np.pi/12.
K1 = 300 # kg/s
K2 = 1000; # kg/s
rho = 1000. # kg/m^3
g = 9.81 # m/s^2

# function defining the governing equation
def gov_eqn(t,h):
    dh = (1.0/(rho*A_tank))*(K1+K2*np.sin(5.0*C*t)*np.cos(C*t) - 
                                rho*A_pipe*np.sqrt(2.0*g*h))   
    
    return dh


a = 0; b = 150;
N = 1000;
tspan = [a,b] # time span, in seconds for the solution

h0 = [3] # m, initial height of the tank

method = 'RK45'

sol = scp.solve_ivp(gov_eqn, tspan, h0, method,dense_output=True,
                    rtol=1e-12,atol=1e-12)

plt.plot(sol.t,sol.y[0,:])
plt.title('Solved with solve_ivp')
plt.show()

# try to solve with my user-defined Classical RK4 method
tCRK,yCRK = odeCRK4(gov_eqn,a,b,N,3.0)

plt.plot(tCRK,yCRK)
plt.title('Solved with CRK4')
plt.show()
