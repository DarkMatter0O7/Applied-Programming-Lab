import numpy as np
import math as mt
from scipy.linalg import lstsq 
from matplotlib.pyplot import *
from sys import argv,exit
from scipy.special import *
import matplotlib.pyplot as plt

PI = mt.pi

# script to generate data files for the least squares assignment
from pylab import *
import scipy.special as sp
N=101                           # no of data points
k=9                             # no of sets of data with varying noise

# generate the data points and add noise
t=linspace(0,10,N)              # t vector
y=1.05*sp.jn(2,t)-0.105*t       # f(t) vector
Y=meshgrid(y,ones(k),indexing='ij')[0] # make k copies
scl=logspace(-1,-3,k)           # noise stdev
n=dot(randn(N,k),diag(scl))     # generate k vectors
yy=Y+n                          # add noise to signal

savetxt("fitting.dat",c_[t,yy]) # write out matrix to file

sigma = logspace(-1,-3,9)

f = open("fitting.dat","r")
lines = f.readlines()
time_arr = [] 


def l(lines,row,col):
	x = lines[row].split()
	return float(x[col])
	
for i in range(101):
	time_arr.append(l(lines,i,0))

data_arr = np.zeros((101,9))



for i in range(101):
	for j in range(1,10):
		data_arr[i][j-1] += l(lines,i,j)

noise = np.zeros((101,9))

for i in range(101):
	for j in range(1,10):
		t = time_arr[i]
		noise[i][j-1] += data_arr[i][j-1] - 1.05*jv(2,t) + 0.105*t

def noise_generate(col,noise):
	gauss = []
	for i in range(101):
		gauss.append((1/sigma*mt.sqrt(2.0*PI))*math.exp((-1.0*(noise[i][col-1]**2))/2.0*sigma*sigma))
	return gauss

def f_data(data_arr,col):
	data_col = []
	for i in range(101):
		data_col.append(data_arr[i][col])
	
	return data_col
	
def noise_data(noise,col):
	noise_col = []
	for i in range(101):
		noise_col.append(noise[i][col])
	
	return noise_col

for i in range(9):
	plot(time_arr,f_data(data_arr,i),label="$\sigma${} = {}".format(i+1,"%.3f" % np.logspace(-1,-3,9)[i]))
	plt.xlabel("t")
	plt.ylabel("f(t)+noise")
	plt.legend(loc = "upper right")

def g_func(t,A,B):			## Function g_func(t;A,B) = A*J2(t) + B*t
	g_arr=[]
	for i in range(101):
		t = time_arr[i]
		g_arr.append(A*jv(2,t)+B*t)
	return g_arr

plt.title("Data to be fitted to theory")	
plot(time_arr,g_func(t,1.05,-0.105),"-k",label="True value")
plt.legend(loc = "upper right")
plt.figure()

plt.xlabel("t")
plt.title("Data points for $\sigma$ = 0.1 along with exact function")
plot(time_arr,g_func(t,1.05,-0.105),"-k",label="f(t)")  ## complete plot of actual data set with std dev = 0.1	
errorbar(time_arr[::5],f_data(data_arr,0)[::5],np.logspace(-1,-3,9)[0],fmt='ro',label="Errorbar") ## Error bar plot of every 5th data item of data set with std dev = 0.1

plt.legend(loc = "upper right")
plt.figure()

A0 = 1.05
B0 = -0.105

J_matrix = np.zeros((101,1))
t_matrix = np.zeros((101,1))
	
for i in range(101):
	t = time_arr[i]
	J_matrix[i][0] += jv(2,t)
	t_matrix[i][0] += t

M = np.array((101,2))
M = c_[J_matrix,t_matrix]


g_matrix = np.zeros((101,1))
for i in range(101):
	t= time_arr[i]
	g_matrix[i][0] += A0*jv(2,t) + B0*t
	
if np.allclose(np.dot(M,np.array([[A0],[B0]])),g_matrix):
	print("The two arrays are equal , Hence Verified!!")
else:
	print("The two arrays are not equal")


error = np.zeros((21,21))
A = [(0.1*i) for i in range(21)]
B = [(-0.2+(0.01*i)) for i in range(21)]


for i in range(21):
	for j in range(21):
		sum_k = 0
		for k in range(101):
			t = time_arr[k]
			sum_k += (1/101)*((data_arr[k][0] - (A[i]*jv(2,t) + B[j]*t))**2)
		
		error[i][j] += sum_k

plt.xlabel("A")
plt.ylabel("B")
plt.title("Contour plot of $\epsilon_{ij}$")			
CS = plt.contour(A,B,error,20)
plt.clabel(CS)

lst_fit = []
for i in range(9):
	lst_fit.append((lstsq(M,data_arr[:,i],rcond = None)[0]))

est_errorA = []
est_errorB = []
for i in range(9):
	est_errorA.append(abs(lst_fit[i][0] - A0))
	est_errorB.append(abs(lst_fit[i][1] - B0))

plt.figure()
plt.xlabel("Noise standard deviation($\sigma$)")
plt.ylabel("Mean square Error")
plt.title("Variation of error with noise")
plot(np.logspace(-1,-3,9),est_errorA,label = "Aerr",marker = 'o',linestyle="-")
plot(np.logspace(-1,-3,9),est_errorB,label = "Berr",marker = 'o',linestyle="-")
plt.legend(loc = "upper left")

plt.figure()
plt.xlabel("$\sigma_{n}$")
plt.ylabel("Mean square error")
plt.title("Variation of error with noise")
plt.loglog(np.logspace(-1,-3,9),est_errorA,label="Aerr",marker = 'o',linestyle='None')
plt.loglog(np.logspace(-1,-3,9),est_errorB,label="Berr",marker = 'o',linestyle='None')
for i in range(9):
	plt.axvline(np.logspace(-1,-3,9)[i])
plt.legend(loc = "upper right")

plt.figure()
plt.xlabel("$\sigma_{n}$")
plt.ylabel("Error returned by lstsq in estimation")
plt.title("Error returned by lstsq vs $\sigma_n$")

error_estimate = []
for i in range(9):
	error_estimate.append((lstsq(M,data_arr[:,i],rcond = None)[1]))
plt.loglog(np.logspace(-1,-3,9),error_estimate,marker='o')

show()	
			
			
			
			
			
				

	


	

 
	
	
	

