'''
AUTHOR: RAKSHIT PANDEY
ROLL NO: EE20B106
DATE: 14 May, 2022
'''

import numpy as np
import scipy.signal as sp
from pylab import *

Pi = np.pi
with open("h.csv") as f: 
	lines = f.readlines()

b = []
for i in range(len(lines)):
	b.append(float(lines[i]))
	
w,h = sp.freqz(lines) # freqz is a scipy.signal attribute to generate the frequency response corresponding to filter, which in this case is FIR filter

# From the curve it is very clear that this frequency reponse is for a low pass filter
figure(1)
xlabel(r"$\omega\rightarrow$")
ylabel("Magnitude")
title("Magnitude response of Low pass filter")
grid()
plot(w,abs(h))

figure(2)
xlabel(r"$\omega\rightarrow$")
ylabel("Phase")
title("Phase response of Low pass filter")
grid()
plot(w,angle(h))

n = arange(1,2**10+1)
x_n = cos(0.2*Pi*n) + cos(0.85*Pi*n)

figure(3)
xlabel("n")
ylabel("x[n]")
title("Input signal")
grid()
plot(n,x_n)

# Output using linear convolution
y_n = np.zeros(len(x_n))
for i in arange(len(x_n)):
    for k in arange(len(lines)):
        y_n[i]+=x_n[i-k]*b[k]

figure(4)
xlabel("n")
ylabel("y[n]")
title("Output signal as a result of linear convolution")
grid()
plot(n,real(y_n))

# Output using Circular convolution 
y1=ifft(fft(x_n)*fft(concatenate((b,zeros(len(x_n)-len(b))))))

figure(5)
xlabel("n")
ylabel("y[n]")
title("Output signal as a result of circular convolution")
grid()
plot(n,real(y1))

#Output of circular convolution using linear convolution

def circular_conv(x,h): # Function defined to calculate the circular convolution using the linear convolution
    P = len(h)
    n_ = int(ceil(log2(P)))
    h_ = np.concatenate((h,np.zeros(int(2**n_)-P)))
    P = len(h_)
    n1 = int(ceil(len(x)/2**n_))
    x_ = np.concatenate((x,np.zeros(n1*(int(2**n_))-len(x))))
    y = np.zeros(len(x_)+len(h_)-1)
    for i in range(n1):
        temp = np.concatenate((x_[i*P:(i+1)*P],np.zeros(P-1)))
        y[i*P:(i+1)*P+P-1] += np.fft.ifft(np.fft.fft(temp) * np.fft.fft( np.concatenate( (h_,np.zeros(len(temp)-len(h_))) ))).real
    return y

y_lc = circular_conv(x_n,b)

figure(6)
xlabel("n")
ylabel("y[n]")
title("Output of circular convolution using linear convolution")
grid()
plot(n,real(y_lc[:1024])) # While plotting we must take care to truncate the output values to take only first 1024 values

# Zadoff-chu sequence
with open("x1.csv") as f1: 
	lines1 = f1.readlines()

b2 = []
for i in range(len(lines1)):
	x = complex((lines1[i].rstrip("\n")).replace("i","j")) # First we rstrip to remove "\n" at the end of the string, then we replace i with j for complex to understand the argument
	b2.append(x)

X = np.fft.fft(b2)
b3 = np.roll(b2,5)
cor = np.fft.ifftshift(np.correlate(b3,b2,'full'))

figure(7)
xlabel("t")
ylabel("Correlation")
title("Auto correlation")
grid()
xlim([0,20])
plot(np.arange(0,len(cor)),abs(cor))
show()




 



