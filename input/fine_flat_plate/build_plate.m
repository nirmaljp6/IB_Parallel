function [ds, xhat, yhat] = build_plate(l, gamma, h, x0, y0, fam)

  %build plate of length l inclided at angle gamma using flow grid spacing of h
  %(cylinder will have a spacing of 2h)
  
  gamma=deg2rad(gamma);
  
  spt=0:round(h,10):l;
  n=length(spt);
  xfun=@(z) x0+z*cos(gamma);
  yfun=@(z) y0+z*sin(gamma);
  xhat=xfun(spt);
  yhat=yfun(spt);



%sanity check: make sure ds is equal to 2 * h
ds = sqrt( (xhat(2) - xhat(1))^2 + (yhat(2) - yhat(1))^2 ) ;



