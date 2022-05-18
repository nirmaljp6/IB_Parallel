clc; clear all;
% Script to plot vorticity and body and forces
% Plots the Vorticity field and velocity fields around the body
% Plots the x and y forces on the body

% -------- General Parameters ---------------------------------------------
dir=pwd;
filename=[dir '/nproc.dat'];
nproc = csvread(filename);
% Number of the time step of the file: this identifies which file should be read
istart = 70476; iend = 74360; inv = 971;
dt = 0.0004375;

for it=istart:inv:iend-inv
% for it=[7100,11700,12700]
    
% Contour maximum values and number of contour levels
% Vorticity
cmax_w = 10;
clev_w = 5;
% Velocity
cmax_v = 1;
clev_v = 21;

% Range for plots
range = [-2 2 -2 2];

% Plot the pressure field?
pressure=1;

% Get data from specified file (see getdata for details) : 1st grid level
[xn, yn, dt, wn, stress, xb, yb, un, vn, u0pn, v0pn, u0rn, v0rn, pn] ...
    = getdata(dir,nproc,it,5,pressure);

u = un + u0pn + u0rn;
v = vn + v0pn + v0rn;
Un = sqrt(u.^2 + v.^2);

figure
pcolor(xn,yn,pn)
colormap(parula)
shading flat
hold on
ittt = 1:10:859;
% quiver(xn(ittt,ittt), yn(ittt,ittt), u(ittt,ittt), v(ittt,ittt), 3,'k');
caxis([-1.2 1.2])
% caxis([-cmax_w cmax_w])
% colorbar
axis equal;

hold on
fill(xb(1:287),yb(1:287),'w');
hold off

hold on
plot(xb(575:603),yb(575:603),'w','linewidth',2)
hold off

xlabel('$x$', 'interpreter', 'latex')
ylabel('$y$', 'interpreter', 'latex')

xlim([-0.1,1.4]); ylim([-0.5,0.25]);
colorbar

drawnow;

% body = [xb(1:287),yb(1:287);xb(575:603),yb(575:603)];
% dlmwrite(['snapshot_' num2str(it) '.csv'],pn,'delimiter', ',', 'precision', 9);
% dlmwrite(['body' num2str(it) '.csv'],body,'delimiter', ',', 'precision', 9);

end
% dlmwrite(['gridx.csv'],xn,'delimiter', ',', 'precision', 9);
% dlmwrite(['gridy.csv'],yn,'delimiter', ',', 'precision', 9);

% for rank=0:nproc-1
% % ---- Open file ----------------------------------------------------------
% forcefile = ['/force' num2str(nproc,'%2.2i') '_' num2str(rank,'%2.2i') '.dat'];   % Filename for Force Data
% force_temp=importdata(strcat(pwd,forcefile));   % Import all data from file
% % force_temp=force_temp(1:23000,:);
% if (rank==0)
%     forces=force_temp;
% else
%     forces(:,2:end)=forces(:,2:end)+force_temp(:,2:end);
% end
% 
% % forcefile = ['/forcerdst' num2str(rank,'%2.2i') '.dat'];   % Filename for Force Data
% % force_temp=importdata(strcat(pwd,forcefile));   % Import all data from file
% % % force_temp=force_temp(1:23000,:);
% % if (rank==0)
% %     forces_rdst=force_temp;
% % else
% %     forces_rdst(:,2:end)=forces_rdst(:,2:end)+force_temp(:,2:end);
% % end
% end
% 
% figure
% subplot(2,1,1)
% plot(forces(:,1)*dt, forces(:,4)+0*forces(:,5))
% ylim([-0.2,2])
% ylabel('C_L')
% % legend('Airfoil','Flap','Location','northeast')
% 
% % subplot(2,1,2)
% % plot(forces(:,1)*dt, forces(:,3))
% % ylim([-0.2,2])
% % ylabel('C_L')
% % % legend('Airfoil','Flap','Location','northeast')
% % 
% % figure
% % subplot(2,1,1)
% % plot(forces(:,1)*dt, forces(:,2))
% % ylim([-0.2,1.2])
% % ylabel('C_D')
% % % legend('Airfoil','Flap','Location','northeast')
% 
% subplot(2,1,2)
% plot(forces(:,1)*dt, forces(:,2)+forces(:,3))
% ylim([-0.2,1.2])
% ylabel('C_D')
% 
% figure
% hold
% for i_bdy=2:2
%     thetfile = ['/theta' num2str(nproc,'%2.2i') '_' num2str(i_bdy,'%2.2i') '.dat'];   % Filename for Force Data
%     thet{i_bdy}=importdata(strcat(pwd,thetfile));   % Import all data from file
%     plot(thet{i_bdy}(:,1)*dt, -thet{i_bdy}(:,2)*180/pi)
% end
% 
% 
% %plot(thet{1}(:,1), thet{1}(:,2))
% 
% tlim=15; tlim2=1000;
% ii=1;
% time = forces(:,1)*dt;
% itime = forces(:,1);
% 
%         cl = forces(:,4);% + forces{ii}(:,5);
%     
%     cl_mean = mean(cl(time>tlim & time<tlim2));
%     cl = max(cl,cl_mean);
% 
% cl_smooth = smooth(cl(itime<=900000),50,'moving');
% % cl_smooth = smooth(cl{ii},1000,'moving');
% [peaks, locs] = findpeaks(cl_smooth);
% 
% figure
% hold
% plot(forces(:,1)*dt,cl)
% plot(locs*dt,peaks,'.')
% ylim([-0.1 2])
% 
% m=5; %number of snapshots to save
% di = locs(end) - locs(end-1);
% irestart = locs(end) + di
% isave = round(di/(m-1))
% istop = irestart + isave*(m-1)
% 
% irestart:isave:istop
