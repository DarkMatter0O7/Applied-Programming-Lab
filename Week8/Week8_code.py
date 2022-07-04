from pylab import *

#Generating spectrum of sin(5*t)
N1 = 128
x1 = linspace(0,2.0*pi,N1+1)
x1 = x1[:-1]
y1=sin(5*x1)
Y1=fftshift(fft(y1))/N1
w1=linspace(-64,63,N1)

#Magnitude spectrum plot for sin(5t)
figure(1)
title("Spectrum of sin(5t)")
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\|Y|\rightarrow$")
xlim([-10,10])
plot(w1,abs(Y1))
grid()

#Phase spectrum plot for sin(5t)
figure(2)
title("Phase plot of sin(5t)")
xlabel(r"k$\rightarrow$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlim([-10,10])
plot(w1,angle(Y1),'ro')
ii=where(abs(Y1)>1e-3)
plot(w1[ii],angle(Y1[ii]),'go')
grid()

N2=512
x2=linspace(-4.0*pi,4.0*pi,N2+1)
x2=x2[:-1]
y2 = (1+0.1*cos(x2))*cos(10*x2)
Y2=fftshift(fft(y2))/N2
w2=linspace(-64,64,N2+1)
w2=w2[:-1]

#Magnitude spectrum plot for (1+0.1cos(t))cos(10t)
figure(3)
title("Spectrum of (1+0.1cos(t))cos(10t)")
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\|Y|\rightarrow$")
xlim([-15,15])
plot(w2,abs(Y2))
grid()

#Phase spectrum plot for (1+0.1cos(t))cos(10t)
figure(4)
title("Phase plot of (1+0.1cos(t))cos(10t)")
xlabel(r"k$\rightarrow$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlim([-15,15])
plot(w2,angle(Y2),'ro')
grid()

N3=512
x3=linspace(-4.0*pi,4.0*pi,N3+1)
x3=x3[:-1]
y3=(sin(x3))**3
Y3=fftshift(fft(y3))/N3
w3=linspace(-64,64,N3+1)
w3=w3[:-1]

#Magnitude spectrum plot for (sin(t))^3
figure(5)
title("Spectrum of $sin^{3}(t)$")
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\|Y|\rightarrow$")
xlim([-10,10])
plot(w3,abs(Y3))
grid()

#Phase spectrum plot for (sin(t))^3
figure(6)
title("Phase plot of $sin^{3}(t)$")
xlabel(r"k$\rightarrow$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlim([-10,10])
plot(w3,angle(Y3),'ro')
grid()

x3=linspace(-4.0*pi,4.0*pi,N3+1)
x3=x3[:-1]
y4=(cos(x3))**3
Y4=fftshift(fft(y4))/N3
w3=linspace(-64,64,N3+1)
w3=w3[:-1]

#Magnitude spectrum plot for (cos(t))^3
figure(7)
title("Spectrum of $cos^{3}(t)$")
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\|Y|\rightarrow$")
xlim([-10,10])
plot(w3,abs(Y4))
grid()

#Phase spectrum plot for (cos(t))^3
figure(8)
title("Phase plot of $cos^{3}(t)$")
xlabel(r"k$\rightarrow$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlim([-10,10])
plot(w3,angle(Y4),'ro')
grid()

N5=512
x5=linspace(-4.0*pi,4.0*pi,N5+1)
x5=x5[:-1]
y5=cos(20*x5+5*cos(x5))
Y5=fftshift(fft(y5))/N5
w5=linspace(-64,64,N3+1)
w5=w5[:-1]

#Magnitude spectrum plot for cos(20t+5cos(t))
figure(9)
title("Spectrum of $cos(20t+5cos(t))$")
xlabel(r"$\omega\rightarrow$")
ylabel(r"$\|Y|\rightarrow$")
xlim([-40,40])
plot(w5,abs(Y5))
grid()

#Phase spectrum plot for cos(20t+5cos(t))
figure(10)
title("Phase plot of $cos(20t+5cos(t))$")
xlabel(r"k$\rightarrow$")
ylabel(r"$\angle Y(\omega)\rightarrow$")
xlim([-40,40])
ii=where(abs(Y5)>1e-3)
plot(w5[ii],angle(Y5[ii]),'go')
grid()

def err(N,range_t):
	t = linspace(-range_t,range_t,N+1)
	t=t[:-1]
	y = exp(-0.5*t**2)
	Y = fftshift(abs(fft(y)))/N
	Y=Y*sqrt(2.0*pi)/max(Y)
	w=linspace(-32,32,N+1)
	w=w[:-1]
	Y_ctft = exp(-0.5*w**2)*sqrt(2.0*pi)
	
	figure()
	title(r"Spectrum of $e^{-t^{2}/2}$")
	xlabel(r"$\omega\rightarrow$")
	ylabel(r"$\|Y|\rightarrow$")
	xlim([-10,10])
	plot(w,abs(Y))
	grid()
	
	figure()
	title(r"Spectrum of $e^{-t^{2}/2}$")
	xlabel(r"k$\rightarrow$")
	ylabel(r"$\angle Y(\omega)\rightarrow$")
	xlim([-10,10])
	plot(w,angle(Y),'ro')
	ii=where(abs(Y)>1e-3)
	plot(w[ii],angle(Y[ii]),'go')
	grid()
	show()
	
	print("For N={} and time window range={}:".format(N,range_t))
	print("The value of maximum error is: {}".format(abs(Y-Y_ctft).max()))
	print("")

err(512,4.0*pi)
err(512,8.0*pi)
err(512,12.0*pi)
err(256,8.0*pi)
err(1024,8.0*pi)







