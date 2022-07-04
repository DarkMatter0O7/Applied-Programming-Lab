from sympy import *
import scipy.signal as sp
import numpy as np
import matplotlib.pyplot as plt

s = symbols('s')
PI = np.pi

#Creating a function for highpass filter
def highpass(R1,R3,C1,C2,G,Vi):
	A = Matrix([[0,-1,0,1/G],[s*C2*R3/(s*C2*R3+1),0,-1,0],[0,G,-G,1],[-s*C2-1/R1-s*C1,0,s*C2,1/R1]])
	b = Matrix([0,0,0,-Vi*s*C1])
	V = A.inv()*b
	return A,b,V
	
#Creating a function for lowpass filter
def lowpass(R1,R2,C1,C2,G,Vi):
	A = Matrix([[0,0,1,-1/G],[-1/(1+s*R2*C2),1,0,0],[0,-G,G,1],[-1/R1-1/R2-s*C1,1/R2,0,s*C1]])
	b = Matrix([0,0,0,-Vi/R1])
	V = A.inv()*b
	return A,b,V

#Creating a function to convert sympy function in a form so that it could be understood by scipy.signal
def num_den(expression):
	num,denom = expression.as_numer_denom()
	num = [float(i) for i in Poly(num,s).all_coeffs()]
	denom = [float(i) for i in Poly(denom,s).all_coeffs()]
	return num,denom

A1,b1,V1=lowpass(10000,10000,1e-9,1e-9,1.586,1/s) 
Vo = V1[3]

num1,den1 = num_den(Vo)
H1 = sp.lti(num1,den1)
t1,Vo1 = sp.impulse(H1,None,np.linspace(0,0.001,1000))

A2,b2,V2 = lowpass(10000,10000,1e-9,1e-9,1.586,1)#For getting Low pass transfer function
Vo2 = V2[3]

num2,den2 = num_den(Vo2) 
H2 = sp.lti(num2,den2) ##Low pass transfer function##
t2 = np.linspace(0,1e-2,100000)
Vi_t = np.sin(2000*PI*t2) + np.cos(2e6*PI*t2)
t2,Vo2,svec = sp.lsim(H2,Vi_t,t2)

A3,b3,V3 = highpass(10000,10000,1e-9,1e-9,1.586,1) #For getting High pass transfer function
Vo3 = V3[3]
num3,den3 = num_den(Vo3)
H3 = sp.lti(num3,den3) ##High pass transfer function##
w,S,phi = H3.bode()

t3 = np.linspace(0,1e-5,1000)
Vi_t2 = np.sin(2000*PI*t3) + np.cos(2e6*PI*t3)
t3,Vo3,svec = sp.lsim(H3,Vi_t2,t3)

def damped_sin(freq,alpha,t):
	return (np.sin(freq*t)*np.exp(-alpha*t))


A4,b4,V4=highpass(10000,10000,1e-9,1e-9,1.586,1/s) 
Vo4 = V4[3]

num4,den4 = num_den(Vo4)
H4 = sp.lti(num4,den4)
t4,Vo4 = sp.impulse(H4,None,np.linspace(0,1e-3,1000))

Vi_decay1 = damped_sin(2*PI,0.5,np.linspace(0,10,1000)) #Low frequency decaying sinusoid
Vi_decay2 = damped_sin(2*PI*1e5,0.5,np.linspace(0,1e-4,1000)) # High frequency decaying sinusoid

#low pass response
t5,Vo5,svec = sp.lsim(H2,Vi_decay1,np.linspace(0,10,1000))
t6,Vo6,svec = sp.lsim(H2,Vi_decay2,np.linspace(0,1e-4,1000))

#High pass response
t7,Vo7,svec = sp.lsim(H3,Vi_decay1,np.linspace(0,10,1000))
t8,Vo8,svec = sp.lsim(H3,Vi_decay2,np.linspace(0,1e-4,1000))

			##PLOTS
plt.figure(1)
plt.title("Vo vs t step function u(t)")
plt.xlabel("t")
plt.ylabel("Vo")
plt.plot(t1,Vo1)
plt.grid()

plt.figure(2)
plt.title("Vo vs t for given Vi")
plt.xlabel("t")
plt.ylabel("Vo")
plt.plot(t2,Vi_t,label="Vin")
plt.plot(t2,Vo2,label="Vo")
plt.legend(loc="upper right")
plt.grid()

# Magnitude response of the given high pass filter
plt.figure(3)
plt.title("Magnitude response of a High Pass filter")
plt.xlabel("frequency")
plt.ylabel("Magnitude response")
plt.semilogx(w,S)
plt.grid()

#Phase response of the given high pass filter
plt.figure(4)
plt.title("Phase response of a High Pass filter")
plt.xlabel("frequency")
plt.ylabel("Phase response")
plt.semilogx(w,phi)
plt.grid()

#Output response of the high pass filter to the sum of sinusoids input 
plt.figure(5)
plt.title("Vo vs t for given Vi")
plt.xlabel("t")
plt.ylabel("Vo")
plt.plot(t3,Vo3,label="Vo")
plt.legend(loc="upper right")
plt.grid()

#response to decaying sinusoids for lowpass filter for lower frequency sinusoid
plt.figure(6)
plt.title("Low pass filter response of decaying sinusoids")
plt.xlabel("t")
plt.ylabel("Response")
plt.plot(t5,Vi_decay1,label="1Hz frequency input")
plt.plot(t5,Vo5,label="Vo for 1Hz frequency")
plt.legend()
plt.grid()

#response to decaying sinusoids for lowpass filter for higher frequency sinusoid
plt.figure(7)
plt.title("Low pass filter response of decaying sinusoids")
plt.xlabel("t")
plt.ylabel("Response")
plt.plot(t6,Vi_decay2,label="10kHz frequency input")
plt.plot(t6,Vo6,label="Vo for 10kHz frequency")
plt.legend()
plt.grid()

#response to decaying sinusoids for highpass filter for lower frequency sinusoid
plt.figure(8)
plt.title("High pass filter response of decaying sinusoids")
plt.xlabel("t")
plt.ylabel("Response")
plt.plot(t7,Vi_decay1,label="1Hz frequency input")
plt.plot(t7,Vo7,label="Vo for 1Hz frequency")
plt.legend()
plt.grid()

#response to decaying sinusoids for highpass filter for higher frequency sinusoid
plt.figure(9)
plt.title("High pass filter response of decaying sinusoids")
plt.xlabel("t")
plt.ylabel("Response")
plt.plot(t8,Vi_decay2,label="10kHz frequency input")
plt.plot(t8,Vo8,label="Vo for 10kHz frequency")
plt.legend()
plt.grid()

#Step response of the high pass filter
plt.figure(10)
plt.title("Step response of the highpass filter")
plt.xlabel("t")
plt.ylabel("Vo")
plt.plot(t4,Vo4)
plt.grid()
plt.show()












	

