function [xn, yn, dt, wn, stress, xb, yb, un, vn, u0pn, v0pn, u0rn, v0rn, pn] = getdata(dir,nproc,it,lev,pressure)

% Reads the data file at istep=it, and returns:
%    xb,yb :     The x,y coordinates of each body/actuator point
%    codeb :     A number identifying the body to which belongs each body/actuator point
%    xn,yn :     Coordinates of each grid point (Cartesian grid)
%    un,vn :     Velocity components at each grid point (based on fluxes from code, which are evaluated at cell faces)
%    wn :        Vorticity at grid point (based on circulation from code, which is evaluated at cell vertices)
%    sn :        Streamfunction at each grid point
%    pn :        Pressure at each grid point
%    filename :  Full name and path of the file read
%
% Input parameters:
%    dir :       Directory where the file is located
%    it :        Number of the time step of the file: this identifies which file should be read
%    lev :       Grid level where data is required
%    pressure :  1 if pressure data is required, 0 if not
%
% Note: The "*.var" files read here have been written to Fortran's 
% unformatted format so in order for these (binary) files to be read, one 
% needs to read several Fortran markers, which are of no interest here.
%
% -------------------------------------------------------------------------

for rank=0:nproc-1
% ---- Open file ----------------------------------------------------------
filename = [dir '/data' num2str(rank,'%2.2i') '_' num2str(it,'%7.7i') '.var'];
fid = fopen(filename,'r');

% ---- Read first line of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    m = fread(fid,1,'int');           % cells in x direction
    n = fread(fid,1,'int');           % cells in y direction
    mg = fread(fid,1,'int');          % total number of grid levels
    nb = fread(fid,1,'int');          % total number of body points
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read second line of file -------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    rey = fread(fid,1,'real*8' );     % Reynolds number
    dt = fread(fid,1,'real*8' );      % Size of time step
    len = fread(fid,1,'real*8' );     % Size of the smallest grid in the x direction
% Note: 'len' sets the size in the y direction too since the grid spacing is uniform
    offsetx = fread(fid,1,'real*8' ); % x-distance from lower-right corner smallest grid to origin
    offsety = fread(fid,1,'real*8' ); % y-distance of lower-right corner of smallest grid to origin
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read third line of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    xsi = fread(fid,1,'int');           % starting grid index  for vorticity in x direction
    xse = fread(fid,1,'int');           % ending grid index  for vorticity in x direction
    xsm = fread(fid,1,'int');          % grid points for vorticity in x direction
    ysi = fread(fid,1,'int');           % starting grid index  for vorticity in y direction
    yse = fread(fid,1,'int');           % ending grid index  for vorticity in y direction
    ysm = fread(fid,1,'int');          % grid points for vorticity in y direction
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read third line of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    xui = fread(fid,1,'int');           % starting grid index  for vorticity in x direction
    xue = fread(fid,1,'int');           % ending grid index  for vorticity in x direction
    xum = fread(fid,1,'int');          % grid points for vorticity in x direction
    yui = fread(fid,1,'int');           % starting grid index  for vorticity in y direction
    yue = fread(fid,1,'int');           % ending grid index  for vorticity in y direction
    yum = fread(fid,1,'int');          % grid points for vorticity in y direction
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read third line of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
    xvi = fread(fid,1,'int');           % starting grid index  for vorticity in x direction
    xve = fread(fid,1,'int');           % ending grid index  for vorticity in x direction
    xvm = fread(fid,1,'int');          % grid points for vorticity in x direction
    yvi = fread(fid,1,'int');           % starting grid index  for vorticity in y direction
    yve = fread(fid,1,'int');           % ending grid index  for vorticity in y direction
    yvm = fread(fid,1,'int');          % grid points for vorticity in y direction
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Compute grid related parameters ------------------------------------

% If we are not considering the smallest grid level, we need to compute the
% grid spacing and x and y offsets for the current grid level. 
%
% Note:  the cells in each grid level are twice as large in both directions
% as the ones of the previous grid level

fac = 2^(lev-1);                                        
delta = len ./ m *fac;                                  % Grid spacing in both directions for current grid
offx = 2^(lev-1) * len/2 - len/2 + offsetx;             % Offset in x direction for current grid
offy = 2^(lev-1) * (n*len/m)/2 - (n*len/m)/2 + offsety; % Offset in y direction for current grid

% ---- Read fourth and fifth lines of file --------------------------------------------
% --> Vorticity
for k=1:mg
temp = fread(fid, 1, 'float32');        % Fortran marker
    omega(ysi-1:yse-1,xsi-1:xse-1,k) =  transpose(reshape(fread(fid, xsm*ysm, 'real*8'),xsm,ysm));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
    stream(ysi-1:yse-1,xsi-1:xse-1,k) =  transpose(reshape(fread(fid, xsm*ysm, 'real*8'),xsm,ysm));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
    qx(yui:yue,xui:xue,k) =  transpose(reshape(fread(fid, xum*yum, 'real*8'),xum,yum));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
   qx0p(yui:yue,xui:xue,k) =  transpose(reshape(fread(fid, xum*yum, 'real*8'),xum,yum));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
    qx0r(yui:yue,xui:xue,k) =  transpose(reshape(fread(fid, xum*yum, 'real*8'),xum,yum));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
    qy(yvi:yve,xvi:xve,k) =  transpose(reshape(fread(fid, xvm*yvm, 'real*8'),xvm,yvm));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
   qy0p(yvi:yve,xvi:xve,k) =  transpose(reshape(fread(fid, xvm*yvm, 'real*8'),xvm,yvm));
temp = fread(fid, 1, 'float32');        % Fortran marker

temp = fread(fid, 1, 'float32');        % Fortran marker
    qy0r(yvi:yve,xvi:xve,k) =  transpose(reshape(fread(fid, xvm*yvm, 'real*8'),xvm,yvm));
temp = fread(fid, 1, 'float32');        % Fortran marker
end

% ---- Read sixth lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');        % Fortran marker
    for k = 1:mg
        rhs_old(ysi-1:yse-1,xsi-1:xse-1,k) = transpose(reshape(fread(fid, xsm*ysm, 'real*8'), xsm,ysm));
    end
temp = fread(fid, 1, 'float32');        % Fortran marker

% ---- Read sixth lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
	rot_angle = fread(fid,1,'real*8' );   % rotating angle of the frame
	rox = fread(fid,1,'real*8' );     % x-coordinate of the center of rotation
	roy = fread(fid,1,'real*8' );     % y-coordinate of the center of rotation
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read seventh lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
	xfm = fread(fid,1,'int' );   % rotating angle of the frame
    xfi = fread(fid,1,'int' );   % rotating angle of the frame
    xfe = fread(fid,1,'int' );   % rotating angle of the frame
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read eigth lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
	stress(xfi:xfe,1) = fread(fid,xfm,'real*8');   % rotating angle of the frame
temp = fread(fid, 1, 'float32');      % Fortran marker
% ---- Read ninth lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
	xb(xfi:xfe,1) = fread(fid,xfm,'real*8');   % rotating angle of the frame
temp = fread(fid, 1, 'float32');      % Fortran marker
% ---- Read tenth lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
	yb(xfi:xfe,1) = fread(fid,xfm,'real*8');   % rotating angle of the frame
temp = fread(fid, 1, 'float32');      % Fortran marker

% ---- Read seventh lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
	mt_body = fread(fid,1,'int' );   % rotating angle of the frame
temp = fread(fid, 1, 'float32');      % Fortran marker
% ---- Read seventh lines of file --------------------------------------------
if mt_body>0
temp = fread(fid, 1, 'float32');      % Fortran marker
    theta = fread(fid,mt_body,'real*8' );
    thetad = fread(fid,mt_body,'real*8' );
    thetadd = fread(fid,mt_body,'real*8' );
temp = fread(fid, 1, 'float32');      % Fortran marker
end
% ---- Read eigth lines of file --------------------------------------------
temp = fread(fid, 1, 'float32');      % Fortran marker
if isempty(temp)==0
stress(xfi:xfe,1) = fread(fid,xfm,'real*8');   % rotating angle of the frame
end
temp = fread(fid, 1, 'float32');      % Fortran marker    


% ---- Done reading file --------------------------------------------------
fclose(fid);

% ---- Pressure Computations ----------------------------------------------
if(pressure == 1)
    % ---- Open file ----------------------------------------------------------
    filename = [dir '/pressure' num2str(rank,'%2.2i') '_' num2str(it,'%7.7i') '.var'];
    fid = fopen(filename,'r');
    % ---- Read first line of file --------------------------------------------
    temp = fread(fid, 1, 'float32');      % Fortran marker
    xpi = fread(fid,1,'int');           % starting grid index  for vorticity in x direction
    xpe = fread(fid,1,'int');           % ending grid index  for vorticity in x direction
    xpm = fread(fid,1,'int');          % grid points for vorticity in x direction
    ypi = fread(fid,1,'int');           % starting grid index  for vorticity in y direction
    ype = fread(fid,1,'int');           % ending grid index  for vorticity in y direction
    ypm = fread(fid,1,'int');          % grid points for vorticity in y direction
    temp = fread(fid, 1, 'float32');      % Fortran marker
    
    for k=1:mg
        temp = fread(fid, 1, 'float32');        % Fortran marker
        pres(ypi:ype,xpi:xpe,k) =  transpose(reshape(fread(fid, xpm*ypm, 'real*8'),xpm,ypm));
        temp = fread(fid, 1, 'float32');        % Fortran marker
    end

end

end

% -Vorticity
omega = omega / delta.^2;

% ---- Combine q and q0 and convert to velocity  --------------------------
% Note:
% Total flux is q+q0 (q0 = q0p + q0r)
% Need to divide by cell face length 'delta' to convert to velocity

for i = 1:mg                                    
    u(:,:,i) = qx(:,:,i)./ delta ;               % x-velocity relative to the body
    u0p(:,:,i) = qx0p(:,:,i)./ delta ;             % x-velocity due to potential flow (due to motion of the grid)
    u0r(:,:,i) = qx0r(:,:,i)./ delta ;             % x-velocity due to rotational flow (due to motion of the grid)
    v(:,:,i) = qy(:,:,i)./ delta ;               % y-velocity relative to the body
    v0p(:,:,i) = qy0p(:,:,i)./ delta ;             % y-velocity due to potential flow (due to motion of the grid)
	v0r(:,:,i) = qy0r(:,:,i)./ delta ;             % y-velocity due to rotational flow (due to motion of the grid)
end

u0 = u0p + u0r; v0 = v0p + v0r;

s0(1,1) = ( -v0(1,1,lev)+u0(1,2,lev) )*delta;    % Compute streamfunction component due to potential flow (due to motion of the grid)
for j = 3:n
    s0(j-1,1) = s0(j-2,1) + u0(j-1,2,lev)*delta; % x-direction
end
for i = 3:m
    s0(:,i-1) = s0(:,i-2) - v0(2:n,i-1,lev)*delta; % y-direction
end

% ---- Interpolate all variables to the same grid -------------------------

% Note: All variables are interpolated to the cell vertices

% Create the grid that will be used for all variables
[xn, yn] = meshgrid(delta:delta:(m-1)*delta, delta:delta:(n-1)*delta);
xn = xn - offx;
yn = yn - offy;

% Grid for x-velocities (vertical cell faces)
[xu, yu] = meshgrid(0:delta:m*delta, delta/2:delta:(n-0.5)*delta); 
xu = xu - offx;
yu = yu - offy;

% Grid for y-velocities (horizontal cell faces)
[xv, yv] = meshgrid(delta/2:delta:(m-0.5)*delta, 0:delta:n*delta);    
xv = xv - offx;
yv = yv - offy;

% Grid for vorticity and streamfunction (cell vertices)
[xw, yw] = meshgrid(delta:delta:(m-1)*delta, delta:delta:(n-1)*delta); 
xw = xw - offx;
yw = yw - offy;

% Interpolate all variables accordingly to xn, yn
un = interp2(xu,yu,u(:,:,lev),xn,yn);
vn = interp2(xv,yv,v(:,:,lev),xn,yn);
u0pn = interp2(xu,yu,u0p(:,:,lev),xn,yn);
v0pn = interp2(xv,yv,v0p(:,:,lev),xn,yn);
u0rn = interp2(xu,yu,u0r(:,:,lev),xn,yn);
v0rn = interp2(xv,yv,v0r(:,:,lev),xn,yn);
wn = interp2(xw,yw,omega(:,:,lev),xn,yn);
sn = interp2(xw,yw,s0(:,:)+stream(:,:,lev),xn,yn);

pn = 0;
if (pressure==1)
% Grid for pressure (cell center)
[xp, yp] = meshgrid(delta/2:delta:(m-0.5)*delta, delta/2:delta:(n-0.5)*delta); 
xp = xp - offx;
yp = yp - offy;
pn = interp2(xp,yp,pres(:,:,lev),xn,yn);
end
end
