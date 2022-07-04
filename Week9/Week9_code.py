from pylab import *
from mpl_toolkits.mplot3d import Axes3D

t1 = linspace(-pi,pi,65); t1=t1[:-1]
dt = t1[1]-t1[0]
fmax=1/dt
y1 = sin(sqrt(2)*t1)
y1[0] = 0
w1 = linspace(-pi*fmax,pi*fmax,65); w1=w1[:-1]
y1 = fftshift(y1)
Y1 = fftshift(fft(y1))/64

figure(1)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$|Y|\rightarrow$")
title("Spectrum of $\sin(\sqrt{2}t)$")
xlim([-10,10])
plot(w1,abs(Y1),lw=2)
grid()

figure(2)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\angle(Y)\rightarrow$")
xlim([-10,10])
plot(w1,angle(Y1),'ro',lw=2)
grid()

#Original function for which we want DFT
t2 = linspace(-3.0*pi,-1.0*pi,65); t2=t2[:-1]
t3 = linspace(pi,3.0*pi,65); t3=t3[:-1]

figure(3)
xlabel(r"$t\rightarrow$")
ylabel(r"$y\rightarrow$")
title("$sin(\sqrt{2}t)$")
xlim([-10,10])
plot(t1,sin(sqrt(2)*t1),'b',lw=2)
plot(t2,sin(sqrt(2)*t2),'r',lw=2)
plot(t3,sin(sqrt(2)*t3),'r',lw=2)
grid()

#function plot with t wrapping every 2*pi
figure(4)
xlabel(r"$t\rightarrow$")
ylabel(r"$y\rightarrow$")
title("$\sin(\sqrt{2}t)$ with t wrapping every 2$\pi$")
xlim([-10,10])
plot(t1,sin(sqrt(2)*t1),'bo',lw=2)
plot(t2,sin(sqrt(2)*t1),'ro',lw=2)
plot(t3,sin(sqrt(2)*t1),'ro',lw=2)
grid()

#For spectrum of a digital ramp
y2 = t1
y2[0]=0
y2=fftshift(y2)
Y2=fftshift(fft(y2))/64

figure(5)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$|Y|$(dB)$\rightarrow$")
title("Spectrum of a digital ramp")
xlim([1,20])
ylim([-20,0])
semilogx(abs(w1),20*log10(abs(Y2)),lw=2)
xticks([1,2,5,10],["1","2","5","10"])
grid()

#sin(sqrt(2)t)*w(t)
n=arange(64)
wnd=fftshift(0.54+0.46*cos(2*pi*n/63))
y=sin(sqrt(2)*t1)*wnd

figure(6)
xlabel(r"$t\rightarrow$")
ylabel(r"$y\rightarrow$")
title("$\sin(\sqrt{2}t)$*w(t) with t wrapping every 2$\pi$")
xlim([-10,10])
plot(t1,y,'bo',lw=2)
plot(t2,y,'ro',lw=2)
plot(t3,y,'ro',lw=2)
grid()

#Spectrum of sin(sqrt(2)*t)*w(t)
y[0]=0
y=fftshift(y)
Y=fftshift(fft(y))/64
w = linspace(-pi*fmax,pi*fmax,65); w=w[:-1]

figure(7)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$|Y|\rightarrow$")
title("Spectrum of $\sin(\sqrt{2}t)$*w(t)")
xlim([-8,8])
plot(w,abs(Y),'b',lw=2)
grid()

figure(8)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\angle(Y)\rightarrow$")
xlim([-8,8])
plot(w,angle(Y),'ro',lw=2)
grid()

# Spectrum plot with more number of points and narrower window of analysis
t4=linspace(-4*pi,4*pi,257);t4=t4[:-1]
dt2=t4[1]-t4[0];fmax2=1/dt2
n2=arange(256)
wnd2=fftshift(0.54+0.46*cos(2*pi*n2/256))
y4=sin(sqrt(2)*t4)
y4=y4*wnd2
y4[0]=0
y4=fftshift(y4) 
Y4=fftshift(fft(y4))/256.0
w4=linspace(-pi*fmax2,pi*fmax2,257);w4=w4[:-1]

figure(9)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$|Y|\rightarrow$")
title("Spectrum of $sin(sqrt{2}t)*w(t)$")
xlim([-4,4])
plot(w4,abs(Y4),'b',lw=2)
grid()

figure(10)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\angle(Y)\rightarrow$")
xlim([-4,4])
plot(w4,angle(Y4),'ro',lw=2)
grid()

#Spectrum of cos(0.86*t)**3 with and without hamming window 
y5=cos(0.86*t4)**3
y5[0]=0
y5=fftshift(y5)
Y5=fftshift(fft(y5))/256

## Plots without windowing
figure(11)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$|Y|\rightarrow$")
title("Spectrum of $cos^{3}(0.86t)}*w(t)$ with hamming window")
xlim([-6,6])
plot(w4,abs(Y5),marker='o',color='b',lw=2)
grid()

figure(12)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\angle(Y)\rightarrow$")
xlim([-6,6])
plot(w4,angle(Y5),'ro',lw=2)
grid()

y5=y5*wnd2
y5[0]=0
y5=fftshift(y5)
Y5=fftshift(fft(y5))/256

## Plots with windowing
figure(13)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$|Y|\rightarrow$")
title("Spectrum of $cos^{3}(0.86t)}*w(t)$ with hamming window")
xlim([-6,6])
plot(w4,abs(Y5),marker='o',color='b',lw=2)
grid()

figure(14)
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\angle(Y)\rightarrow$")
xlim([-6,6])
plot(w4,angle(Y5),'ro',lw=2)
grid()
show()
























