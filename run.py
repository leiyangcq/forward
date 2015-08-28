__author__ = 'yangleicq'


import numpy as np

import matplotlib.pyplot as mp

import scipy.integrate as integrate

import matplotlib.pyplot as plt

import sympy


x=sympy.symbols('x')
#This will take some time because we are evaluating oscillatory function integration
an,bn=getcoef(expr=x**3)

#check the validity of Fourier expansion
a=0
b=1
an,bn=coef(case=1)


print "Error = "
print YT(b,an,bn)-YTOrigin(b,case=1)
print integrate.quad(YT, a, b,args=(an,bn))[0]  -  integrate.quad(YTOrigin, a, b,args=(5,1))[0]



np.random.seed(123)

QoI = ([getInt_random(an,bn,E=0.5,p=-2) for i in range(5000)])




#histogram
plt.hist(QoI,bins=30,normed=True)


#CDF
sorted_data = np.sort(QoI)
cumulative = np.cumsum(sorted_data)/np.sum(QoI)
mp.plot(sorted_data, np.linspace(0,1,sorted_data.size))





