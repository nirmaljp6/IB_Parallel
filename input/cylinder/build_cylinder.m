function [ds, xhat, yhat] = build_cylinder(r, h, fam)

  %build cylinder of radius r using flow grid spacing of h
  %(cylinder will have a spacing of 2h)

  circum = 2 * pi * r; %Circumference of the circle


  %Get # of points such that ds = 2h
  n = floor( circum / h / 2 );

  int =  2*pi / n ;
  spt = 0 : int : (n-1)*int;
  xfun = @(z) r*cos(z);
  yfun = @(z) r*sin(z);
  xhat = xfun(spt);
  yhat = yfun(spt);

plot(xhat,yhat,'b')

%sanity check: make sure ds is equal to 2 * h
ds = sqrt( (xhat(2) - xhat(1))^2 + (yhat(2) - yhat(1))^2 ) ;

%Write the points to a data file:

fileID = fopen('body.001.inp','w');
fprintf(fileID,'%-6d \n',n);
fprintf(fileID,'%-1s \n','F');
fprintf(fileID,'%-6d \n',fam);
for j = 1 : n
fprintf(fileID,'%-20.16f %-20.16f\n',xhat(j), yhat(j));
end
fclose(fileID);

%plot(xhat,yhat)
