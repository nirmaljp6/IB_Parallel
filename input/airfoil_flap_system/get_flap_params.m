%gets data points for 1 rigid airfoil and multiple rigid or torsional flaps on it
clc
clear all

loc = [0.25,0.5,0.75]; %query location
addflush = 5; %how much flush to add
flap_len = 0.2;
fam = 2; %family for the flap
alpha = 20;
addtranslatex = 0;%-0.14; %how much translation to add
addtranslatey = 0*1.18; %how much translation to add
ii = 1; %rigid body number 
rigid_bodies = 1; %total rigid bodies to create

len = 3;
M = 860;
ds = len/M*2;


%build baseline airfoil
[xairfoilbase,yairfoilbase,ds]= naca00_ibfm(12,ds,1,0,'body.001.inp',1);
%apply final rotation
[xairfoil,yairfoil,ds]= naca00_ibfm(12,ds,1,alpha,'body.001.inp',1);
figure(100)
hold on

for i=1:length(loc)
    clear xhat yhat
%nearest point
[xmin, idx] = min(abs(xairfoilbase(1:round(length(xairfoilbase)/2))-loc(i)));  %checking only on the leading edge

%slope at given point
flush = atand((yairfoilbase(idx+1)-yairfoilbase(idx-1))/(xairfoilbase(idx+1)-xairfoilbase(idx-1)));
flush = flush + addflush; 

%get parameters for the flap
xpivot = xairfoil(idx);
ypivot = yairfoil(idx);
flush = flush-alpha; %provide inclination to the flap
[ds, xhat, yhat] = build_plate(flap_len, flush, ds, xpivot, ypivot, 2);
xhat = xhat + addtranslatex;
yhat = yhat + addtranslatey;
plot(xhat,yhat,'b')

%Write the flap body points to a data file:
fileID = fopen(['body.' num2str(rigid_bodies+length(loc)*(ii-1)+i,'%3.3i') '.inp'],'w');
fprintf(fileID,'%-6d \n',length(xhat));
fprintf(fileID,'%-6d \n',fam);

if fam==2
    aa=0.0; %an extra angle
    itheta=0.001;
    ktheta=0.001;
    ctheta= 0;
    fprintf(fileID,'%-20.16f \n',flush);
    fprintf(fileID,'%-20.16f \n',aa);
    fprintf(fileID,'%-20.16f \n',itheta);
    fprintf(fileID,'%-20.16f \n',ktheta);
    fprintf(fileID,'%-20.16f \n',ctheta);
end

for j = 1 : length(xhat)
fprintf(fileID,'%-20.16f %-20.16f\n',xhat(j), yhat(j));
end
fclose(fileID);

idxs(i) = idx;
end

xairfoil(idxs) = [];
yairfoil(idxs) = [];
xairfoil = xairfoil + addtranslatex;
yairfoil = yairfoil + addtranslatey;
plot(xairfoil,yairfoil,'r')
axis equal

%save airfoil
fid = fopen(['body.' num2str(ii,'%3.3i') '.inp'],'w');
fprintf(fid,'%-6d \n',length(xairfoil));
fprintf(fid,'%-6d \n',1);

for k=1:length(xairfoil)
fprintf(fid,'%-20.16f %-20.16f\n', xairfoil(k),yairfoil(k));
end
fclose(fid);