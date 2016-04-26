# glass pla!te changes pathlength
n = 1.4
l = 800
a= linspace(1, 40, 100) * pi / 180
b = arcsin(sin(a)/n)
r = a - b
d = 3000000             # nm
p = (( n * d / cos(b)  + d ) - (n * d + d * cos(r) / cos(b)))/l * 4 * 180 / pi
figure;plot(p-p[1])
dd = (( n * d / cos(b)  + d ) - (n * d + d * cos(r) / cos(b))) / 1000
plot(diff(dd))
