import numpy as np
from numpy import *
import scipy.signal as sp
import matplotlib.pyplot as plt
from math import cos,sin,exp

Num1 = poly1d([1,0.5])
Den1 = polymul([1,0,2.25],[1,1,2.50])
H1 = sp.lti(Num1,Den1)
t1,x1 = sp.impulse(H1,None,linspace(0,50,501))

plt.figure(1)
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Time response of a spring for decay = 0.5")
plt.plot(t1,x1)

Num2 = poly1d([1,0.05])
Den2 = polymul([1,0,2.25],[1,0.1,2.2525])
H2 = sp.lti(Num2,Den2)
t2,x2 = sp.impulse(H2,None,linspace(0,50,501))

plt.figure(2)
plt.xlabel("t")
plt.ylabel("x(t)")
plt.title("Time response of spring for decay = 0.05") 
plt.plot(t2,x2)

freq = 1.40
for i in range(5):
	f= []
	F_num = poly1d([1,0.05])
	F_den = polymul([1,0.1,freq**2+0.0025],[1,0,2.25])
	H = sp.lti(F_num,F_den)
	t = linspace(0,100,1001)
	for i in range(len(t)):
		f.append(cos(freq*t[i])*exp(-0.05*t[i]))
	t,x,vec = sp.lsim(H,f,t)
	plt.figure()
	plt.xlabel("t")
	plt.ylabel("x(t)")
	plt.title("Time response of a spring for frequency = {}".format(freq))
	freq += 0.05
	plt.plot(t,x)
	plt.show()

num3 = poly1d([1,0,2])
den3 = poly1d([1,0,3,0])
H3 = sp.lti(num3,den3)
t3,x3 = sp.impulse(H3,None,linspace(0,20,201))

num4 = poly1d([2])
den4 = poly1d([1,0,3,0])
H4 = sp.lti(num4,den4)
t4,y4 = sp.impulse(H4,None,linspace(0,20,201))

plt.figure(8)
plt.xlabel("t")
plt.ylabel("x(t) and y(t)")
plt.title("Coupled spring oscillation")
plt.plot(t3,x3,label="x(t)")
plt.plot(t4,y4,label="y(t)")
plt.legend(loc="upper right")

R = 100
L = 1e-6
C = 1e-6
s_2 = L*C
s_1 = R*C

num5 = poly1d([1])
den5 = poly1d([s_2,s_1,1])
H5 = sp.lti(num5,den5)
w,S,phi = H5.bode()

t_30 = linspace(0,30e-6,10000)

plt.figure(9)
plt.xlabel("Frequency")
plt.ylabel("Magnitude of transfer function")
plt.title("Magnitude Response")
plt.semilogx(w,S)
plt.grid()

plt.figure(10)
plt.xlabel("Frequency")
plt.ylabel("Phase of transfer function")
plt.title("Phase Response")
plt.semilogx(w,phi)
plt.grid()

t_30m = np.linspace(0,30e-3,10000) 

Vi = []
Vi2 = []
for i in range(len(t_30)):
	Vi.append(cos(1e3*t_30[i]) - cos(1e6*t_30[i]))
	
for i in range(len(t_30m)):
	Vi2.append(cos(1e3*t_30m[i]) - cos(1e6*t_30m[i]))

t5,V0,vec2 = sp.lsim(H5,Vi,t_30)

t6,V02,vec3 = sp.lsim(H5,Vi2,t_30m)

plt.figure(11)
plt.xlabel("t")
plt.ylabel("$V_{0}(t)$")
plt.title("Output signal for 0 < t < 30$\mu$s")
plt.plot(t5,V0)
plt.grid()

plt.figure(12)
plt.xlabel("t")
plt.ylabel("$V_{0}(t)$")
plt.title("Output signal for 0 < t < 30ms")
plt.plot(t6,V02)
plt.grid()

plt.show()








