import sys
import math
import numpy
import matplotlib.pyplot as plt
from scipy import integrate as integrate

# gravitational constant (aka 9.81)
g = 9.80665
# max speed of the train in m/s
Vmax = 20.11111

# user inputs
# distance between stations
c = 0
while c == 0:
    try:
        c = float(input('Horizontal distance between stations (m): '))
    except ValueError:
        print("Not a valid number! Please try again")

# revision
z = 32
# height difference between stations
h = 0
while h == 0:
    try:
        h = float(input('Vertical distance between stations (m): '))
    except ValueError:
        print("Not a valid number! Please try again")

h = abs(h)
h = -1 * h

# time between stations

# mass of cart
m = 327000
# added energy


a = (5*z)/(c**4)
b = (a*c**4-12*h)/(2*a*c**3)


def f(x):
    y = -1*a*x*(x-b)*(x-c)
    return y


k = numpy.ceil((12*h*c)/(-5*c))+20
z = k

a = (5*z)/(c**4)
b = (a*c**4-12*h)/(2*a*c**3)

j = abs(integrate.quad(f, 0, b)[0])


while numpy.sqrt(2*g*j) < Vmax:
    a = (5*z)/(c**4)
    b = (a*c**4-12*h)/(2*a*c**3)

    j = abs(integrate.quad(f, 0, b)[0])

    z = z+1

z = z-2


a = (5*z)/(c**4)
b = (a*c**4-12*h)/(2*a*c**3)

j = abs(integrate.quad(f, 0, b)[0])


D1 = numpy.hypot(b, j)
E1 = m*g*j
V1 = numpy.sqrt(2*g*j)
T1 = (2*D1)/V1


D2 = numpy.hypot(c-b, j+h)
E2 = m*g*(j+h)
#V2 = numpy.sqrt(2*((E2)/m))
T2 = (2*D2)/V1

Et = E1-E2
Tmin = T1+T2

Tt = Tmin

t = 0
while t == 0 or t < Tmin:
    print("Time to travel between stations (s)\nMinimum is ", end='')
    print(Tmin, end='')
    try:
        t = float(input(': '))
    except ValueError:
        print("Not a valid number! Please try again")
    if t < Tmin:
        print("Time must be greater than the minimum!\n")


D1 = numpy.hypot(c-b, j+h)
E1 = m*g*(j+h)
V1 = numpy.sqrt(2*((E1+Et)/m))
T1 = (2*D1)/V1


D2 = numpy.hypot(b, j)
E2 = m*g*(j)
#V2 = numpy.sqrt(2*((E2)/m))
T2 = (2*D2)/V1


while z >= k and Tt < t:

    a = (5*z)/(c**4)
    b = (a*c**4-12*h)/(2*a*c**3)

    j = abs(integrate.quad(f, 0, b)[0])

    D1 = numpy.hypot(b, j)
    E1 = m*g*j
    V1 = numpy.sqrt(2*g*j)
    T1 = (2*D1)/V1

    D2 = numpy.hypot(c-b, j+h)
    E2 = m*g*(j+h)
    #V2 = numpy.sqrt(2*((E2)/m))
    T2 = (2*D2)/V1

    Et = E1-E2
    Tt = T1+T2

    z = z-1


z = z+2

a = (5*z)/(c**4)
b = (a*c**4-12*h)/(2*a*c**3)

j = abs(integrate.quad(f, 0, b)[0])

D1 = numpy.hypot(b, j)
E1 = m*g*j
V1 = numpy.sqrt(2*g*j)
T1 = (2*D1)/V1

D2 = numpy.hypot(c-b, j+h)
E2 = m*g*(j+h)
#V2 = numpy.sqrt(2*((E2)/m))
T2 = (2*D2)/V1

Et = E1-E2
Tt = T1+T2


# energy use of straight line
D3 = numpy.hypot(c/2, abs(h)/2)
a = (2*D3)/((t/2)**2)
D4 = numpy.hypot(c, abs(h))/2
E3 = a*D4*m
E4 = m*g*abs(h)
Eline = E3+E4

X = numpy.arange(0, c, 0.01)


def F(x):
    res = numpy.zeros_like(x)
    for i, val in enumerate(x):
        y, err = integrate.quad(f, 0, val)
        res[i] = y
    return res


plt.plot(X, F(X))

print("Ideal path is represented by one where z = ", int(z))
print("It takes exactly", int(Tt), " seconds to go between the stations")
print("Energy use of this path is", int(Et/1000), "kilojoules")
print("Energy use of a straight line is ", int(Eline/1000), "kilojoules")
print("Here is the plot:")


plt.show()
