function [x,y,ds] = build_ellipse(dist,chord_length,aspect_ratio,x0,y0,alpha,filename,fam)

% spacing
ds = max(min(abs(dist),.1),10*eps);

% chord length
c = max(chord_length, ds);

%iteration number
it=0; 
err=1;

while err>1e-4
    
    a = [c/2, c/2/aspect_ratio]; 

    % function for computing the upper boundary as function of x (centerline)
    ubd = @(x)( sqrt(1-(x/a(1)).^2)*a(2) );    

    % half arc length
    arcTot = quadl( ...
      @(x)(sqrt(1+(a(2)/a(1)^2*(1-x.^2/a(1)^2).^(-0.5).*(-x)).^2 )),-c/2,0 );

    % upper bound on equally spaces
    ds = ds+ds/10000;
    N = ceil(arcTot/ds);

    % allocate temporary array for points
    XY = zeros(2,N);

    % find the first point
    xf=-c/2;

    % first point on upper boundary
    XY(1,1) = xf;
    XY(2,1) = ubd(xf);

    % the rest of the points on the upper boundary
    for k=2:N
        % function used to determine the next x value
        minFnc = @(x)( (x-XY(1,k-1)).^2 + ( ubd(x) - XY(2,k-1)).^2 - ds.^2);

        % choose next value of x such that new point is a disntace 
        % DS from previous
        [xf,fval,exitFlag] = fzero(minFnc,XY(1,k-1)+0.5*ds);

        % add to array
        XY(1,k) = xf;
        XY(2,k) = ubd(xf);
    end
    err = abs(XY(1,end)); 
    it = it+1;
end

XY(1,N+1:2*N-1) = -flip(XY(1,1:N-1));
XY(2,N+1:2*N-1) = flip(XY(2,1:N-1));
N=size(XY,2);

XY = [XY(1,:), XY(1,N-1:-1:2); XY(2,:), -XY(2,N-1:-1:2)]; %[XY(1,2*N:-1:1),XY(1,1:2*N);XY(2,:-1:1),-XY(2,1:cnt)];

% % rotate by alpha
rotMat = [ cos(alpha/180*pi), sin(alpha/180*pi); -sin(alpha/180*pi), cos(alpha/180*pi) ];
XY = rotMat*XY;

XY(1,:) = XY(1,:) + x0;
XY(2,:) = XY(2,:) + y0;
% 
% prep outputs
x= XY(1,:);
y = XY(2,:);

plot(x,y,'.k')
axis equal


% write surafce boundary output file
if nargin>4
  fid = fopen(filename,'w');
  % number of points
  fprintf(fid,'%-6d \n',size(XY,2));
  fprintf(fid,'%-6d \n',fam);
  
  for k=1:size(XY,2)
    fprintf(fid,'%-20.16f %-20.16f\n', XY(1,k),XY(2,k));
  end
  fclose(fid);
end

end


