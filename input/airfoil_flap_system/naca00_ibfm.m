function [x,y,ds] = naca00_ibfm(thickness,dist,chord_length,alpha,filename,fam)

%NACA00 - Generate NACA00XX airfoil
% Uniform distribution of points define the lower and upper boundary
% of a NACA00XX airfoil.
% 
% Syntax:  [x,y,ang] = naca00(thickness,ds,alpha,filename)
%
% Inputs:
%    thickness - airfoil thickness as percent of chord
%    dist        - nominal spacing between consecutive points
%    chord_length
%    alpha     - angle of attack (degrees)
%    filename  - write output file
%
% Outputs:
%    x,y - arrays with the (x,y) locations of points defining airfoil 
%

% Author: Sebastian Liska
% email: sebastian.liska@gmail.com
% July 2011; Last revision: 17-July-2011

% Modified: Jeesoon Choi
% email: jeesoonchoi@naver.com
% Added some functions to rescale the chord length 


% spacing
ds = max(min(abs(dist),.1),10*eps);

% thickness (% of c)
thck = min(abs(thickness),100)/100;

% chord length
c = max(chord_length, ds);
chord_old=0;
chord_new=0;

% iteration to get chord length scaled to 1 after the cuting the trailing
% edge

%iteration number
it=0; 
it_max=10e1;

while(1)
    
if(abs(chord_length - chord_new) <= abs(chord_length - chord_old) || abs(chord_length - chord_new) >= 0.02)
    
    chord_old=chord_new;
    c= c + 0.001;
    
    if(it > it_max)
        fprintf('exceeded the maximum iteration which is 10e4')
        break;
    end 
    
    % naca coefficients
    % note that stadard coefficients are:
    % [0.2969, -0.1260, -0.3516, 0.2843, -0.1015]
    % they are not used since the airfoil does not close at trailing edge
    a = 5*thck*c*[0.2969, -0.1260, -0.3516, 0.2843, -0.1036]; 

    % function for computing the upper boundary as function of x (centerline)
    ubd = @(x)( a(1)*sqrt(x/c) + a(2)*(x/c) + a(3)*(x/c).^2 + a(4)*(x/c).^3 + a(5)*(x/c).^4 );

    % derivative coefficients
    b = a/c.*[0.5, 1, 2, 3, 4];

    % tot arc length
    arcTot = quadl( ...
      @(x)(sqrt(1+( b(1)./sqrt(x/c)+b(2)+b(3)*(x/c)+b(4)*(x/c).^2+b(5)*(x/c).^3 ).^2 )),0,1 );

    % upper bound on equally spaces
    N = ceil(1.1*arcTot/ds);

    % allocate temporary array for points
    XY = zeros(2,N);
    XY_prime = zeros(2,N);

    % find the first point
    % airfoil symmetric about centerline
    % two points a distance DS apart define leading edge
    fnc =@(x)( ubd(x) - ds/2);
    [xf,fval,exitFlag] = fzero(fnc,[2*eps,ds]);

    % first point on upper boundary
    XY(1,1) = xf;
    XY(2,1) = ubd(xf);

    % the rest of the points on the upper boundary
    for k=2:N
      if ( XY(1,k-1)>0.5 && (XY(2,k-1)<ds/2) )
        cnt = k-1-1;
        chord_new=XY(1,cnt);
        break;
      else

        % function used to determine the next x value
        minFnc = @(x)( (x-XY(1,k-1)).^2 + ( ubd(x) - XY(2,k-1)).^2 - ds.^2);

        % choose next value of x such that new point is a disntace 
        % DS from previous
        [xf,fval,exitFlag] = fzero(minFnc,[XY(1,k-1),XY(1,k-1)+ds]);

        % add to array
        XY(1,k) = xf;
        XY(2,k) = ubd(xf);

      end
    end

    % check the upper boundary points to remove the points that are closer than
    % the distance(=grid space)
    
    it=it+1;

else    
    

% get lower boundary
XY = [XY(1,cnt:-1:1),XY(1,1:cnt);XY(2,cnt:-1:1),-XY(2,1:cnt)];

% rotate by alpha
rotMat = [ cos(alpha/180*pi), sin(alpha/180*pi); -sin(alpha/180*pi), cos(alpha/180*pi) ];
XY = rotMat*XY;

% prep outputs
x= XY(1,:);
y = XY(2,:);



break;

end

end

end


