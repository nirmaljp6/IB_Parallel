clc
clear all

offsetx = 0.5; offsety = 1.5;
len = 3; m = 860; n = m;

delta = len/m;
x = -offsetx:delta:-offsetx + len;
y = -offsety:delta:-offsety + len*n/m;

deltasd = delta*1

[~,i] = min(abs(x-(-0.001-4*deltasd)))
sdxi = x(i)

[~,i] = min(abs(x-(1+4*deltasd)))
sdxe = x(i)

[~,i] = min(abs(y-(-0.34-4*deltasd)))
sdyi = y(i)

[~,i] = min(abs(y-(0.22+4*deltasd)))
sdye = y(i)
