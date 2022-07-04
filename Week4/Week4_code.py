import numpy as np
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
from scipy.integrate import quad
from math import sin,cos,pi,exp

def f1(x):
	return exp(x)

def f2(x):
	return cos(cos(x))
	
t = np.arange(-2.0*pi,4.0*pi,0.01)

y1 = []
y2 = []

for i in range(len(t)):
	y1.append(f1(t[i]))
	y2.append(f2(t[i]))

a0_f1 = (1/(2.0*pi))*quad(f1,0,2.0*pi)[0]
a0_f2 = (1/(2.0*pi))*quad(f2,0,2.0*pi)[0]

c1_f1 = [] 
c2_f2 = []

def u1(t,k):
	return f1(t)*cos(k*t)
def u2(t,k):
	return f2(t)*cos(k*t)

def v1(t,k):
	return f1(t)*sin(k*t)
def v2(t,k):
	return f2(t)*sin(k*t)  

c1_f1.append(a0_f1)
c2_f2.append(a0_f2)

for k in range(1,26):
	c1_f1.append((1/pi)*quad(u1,0,2.0*pi,args=(k))[0])
	c1_f1.append((1/pi)*quad(v1,0,2.0*pi,args=(k))[0])
	
	c2_f2.append((1/pi)*quad(u2,0,2.0*pi,args=(k))[0])
	c2_f2.append((1/pi)*quad(v2,0,2.0*pi,args=(k))[0])
	
f1_list = []
f2_list = []

for i in range(51):	
	f1_list.append(c1_f1[i])
	f2_list.append(c2_f2[i])

f1_vec = np.reshape(f1_list,(-1,1))
f2_vec = np.reshape(f2_list,(-1,1))

n = [i+1 for i in range(51)]

x=np.linspace(0,2*pi,401)
x=x[:-1] 

b1 = np.zeros((400,1)) 
b2 = np.zeros((400,1))

for i in range(400):
	b1[i][0] += exp(x[i])
	b2[i][0] += cos(cos((x[i])))

A = np.zeros((400,51))
A[:,0] = 1.0 

for k in range(1,26):
	row = 0
	while row<400: 
		A[row][2*k-1] += cos(k*x[row]) 
		A[row][2*k] += sin(k*x[row])
		row += 1 

c1 = lstsq(A,b1)[0] 
c2 = lstsq(A,b2)[0]

		## semilog plot of exp(x) vs x from -2.0*pi to 4.0*pi
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Plot of $e^{x}$ vs x in loglog scale")
plt.semilogy(t,y1)
plt.grid()
plt.figure()
	
		## plot of cos(cos(x)) vs x from -2.0*pi to 4.0*pi
plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Plot of cos(cos(x)) vs x between 2$\pi$ and 4$\pi$")
plot(t,y2)
plt.grid()
plt.figure()
		## semilog plot for exp(x) 
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title("Plot of coefficient$_{n}$ vs n for $e^{x}$ on semilog scale")
plt.semilogy(n,list(map(abs,c1_f1)),'ro',label="True value")
plt.legend(loc="upper right")
plt.grid()
plt.figure()
		## loglog plot for exp(x)
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title("Plot of coefficient$_{n}$ vs n for $e^{x}$ on loglog scale")
plt.loglog(n,list(map(abs,c1_f1)),'ro',label="True value")
plt.legend(loc="upper right")
plt.grid()
plt.figure()

		## semilog plot for cos(cos(x))
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title("Plot of coefficient$_{n}$ vs n for cos(cos(x)) on semilog scale")
plt.semilogy(n,list(map(abs,c2_f2)),'ro',label="True value")
plt.legend(loc="upper right")
plt.grid()
plt.figure()

		## loglog plot for cos(cos(x))
plt.xlabel("n")
plt.ylabel("fourier coefficients")
plt.title("Plot of a$_{n}$ and b$_{n}$ vs n for cos(cos(x)) on loglog scale")
plt.loglog(n,list(map(abs,c2_f2)),'ro',label="True value")
plt.legend(loc="upper right")
plt.grid()
plt.figure()

f1_approx = np.dot(A,c1)
f2_approx = np.dot(A,c2)

y1_new = []
y2_new = []

for i in range(len(x)):
	y1_new.append(f1(x[i]))
	y2_new.append(f2(x[i]))

plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("$e^{x}$ vs x in semilog scale and the plot of function obtained by approximation")
plt.semilogy(x,y1_new,'ro',label="True value")
plt.semilogy(x,f1_approx,'go',label="Function approximation")
plt.grid()
plt.legend(loc="upper right")
plt.figure()

plt.xlabel("x")
plt.ylabel("f(x)")
plt.title("Plot of cos(cos(x)) vs x and the plot of function obtained by approximation")
plot(x,y2_new,'ro',label="True value")
plot(x,f2_approx,'go',label="Function approximation")
plt.grid()
plt.legend(loc="upper right")

f1_dev = []
f2_dev = []

for i in range(51):
	f1_dev.append(abs(c1[i] - c1_f1[i]))
	f2_dev.append(abs(c2[i] - c2_f2[i]))

print("For exp(x) :")
print("	Maximum deviation is: {}".format(max(f1_dev)))
print(" ")
print("For cos(cos(x)) :")
print("	Maximum deviation is: {}".format(max(f2_dev)))

show()
	







