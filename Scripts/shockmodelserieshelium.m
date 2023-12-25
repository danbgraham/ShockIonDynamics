Vsw = -320e3; 
n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;

l = 10e3; 

xvec = [-10e2:1:5e2]*1e3;
[xx,yy] = meshgrid(xvec,xvec);
Byxx = -B0*tanh(xx/l)+B1;

By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;

Units = irf_units;
mu0 = Units.mu0;
e = Units.e;
mi = Units.mp;

Jz = -B0*sech(xvec/l).^2/(mu0*l);
Vx = Vsw*nsw./ni;

Emultiplier = 1;

Ex = -Jz.*By./(e*ni) + P0*sech(xvec/l).^2./(e*ni*l);
Ez = -Vx.*By;

dx = 1e3;
phi = cumsum(Ex)*dx;

opts = odeset('RelTol',1e-8,'AbsTol',1e-10);
%%
if 1
%plot(d1(:,1),d1(:,2))
%hold on
%plot(d1b(:,1),d1b(:,2))
dv = 10e3;
vxvec = -1300e3:dv:700e3;
vzvec = -1000e3:dv:1000e3;
xpos = (-2000:10:200)*1e3;
%xpos = (-50:10:100)*1e3;

plotpos1 = -400e3;
plotpos2 = 50e3;
[~,idx1] = min(abs(xpos-plotpos1));
[~,idx2] = min(abs(xpos-plotpos2));

[VX,VZ] = meshgrid(vxvec, vzvec);
fxz = zeros(size(VX));
fxvxvz = repmat(fxz,1,1,length(xpos));
fxvxvz = squeeze(permute(fxvxvz,[3 1 2]));

fxvxvztemp = squeeze(zeros(size(fxvxvz(1,:,:))));

Ti = 12;
vi = sqrt(2*e*Ti/mi/4);
nsw = 9e6;
vsw = -320e3;
spos = 200e3;
y0 = 0;
z0 = 0;
vy0 = 0;

%%

for kk = 1:length(xpos)

x0 = xpos(kk);  
tic
for ii = 1:length(vxvec)
	parfor jj = 1:length(vzvec)
    [t1b,d1b] = ode45(@temp,[20 0],[x0; vxvec(ii); y0; vy0; z0; vzvec(jj)],opts);
    if d1b(end,1) > spos
      [a,~] = find(d1b(:,1) > spos);
      idx = a(1)-1;
      fxvxvztemp(jj,ii) = nsw/(pi*vi^2)*exp((-(d1b(idx,2) - vsw)^2-d1b(idx,6)^2)/vi^2);
    else
      fxvxvztemp(jj,ii) = 0;
    end
  end

fxvxvz(kk,:,:) = fxvxvztemp;
end
toc

kk
end

fdisthelium = struct('xpos',xpos,'vxvec',vxvec,'vzvec',vzvec,'fxvxvz',fxvxvz,'Emultiplier',Emultiplier,'dv',dv);
save('fdisthelium.mat','fdisthelium')

%%
fn=figure;
set(fn,'Position',[10 10 800 750])
h(1)=axes('position',[0.065 0.89 0.85 0.10]);
h(2)=axes('position',[0.065 0.78 0.85 0.10]);
h(3)=axes('position',[0.065 0.63 0.905 0.14]);
h(4)=axes('position',[0.065 0.48 0.917 0.14]);
h(5)=axes('position',[0.08 0.06 0.40 0.35]); % [x y dx dy]
h(6)=axes('position',[0.58 0.06 0.40 0.35]);
set(fn,'defaultLineLineWidth',2);
set(fn,'defaultAxesFontSize',12)

plot(h(1),xvec/1e3,By*1e9,'k');
hold(h(1),'on')
plot(h(1),[xpos(idx1) xpos(idx1)]/1e3,[0 65],'m--');
plot(h(1),[xpos(idx2) xpos(idx2)]/1e3,[0 5],'m--');
hold(h(1),'off')
set(h(1),'Ylim',[30 65]);
set(h(1),'XTickLabel',[])
ylabel(h(1),'B_y (nT)','FontSize',12)
ax2 = axes('Position',get(h(1),'Position'),...
		'XAxisLocation','top',...
		'YAxisLocation','right',...
		'Color','none',...
		'XColor','k','YColor','r',...
		'XTickLabel',[],'Xtick',[]);
axes(ax2)
line(xvec/1e3,ni/1e6,'Color','r','Parent',ax2);
set(ax2,'Ylim',[0 35]);
ylabel('n_e (cm^{-3})','FontSize',12)

plot(h(2),xvec/1e3,Ex*1e3,'k');
ylabel(h(2),'E_x (mV m^{-1})','FontSize',12)
set(h(2),'XTickLabel',[])
hold(h(2),'on')
plot(h(2),[xpos(idx1) xpos(idx1)]/1e3,[0 35],'m--');
plot(h(2),[xpos(idx2) xpos(idx2)]/1e3,[0 35],'m--');
hold(h(2),'off')
set(h(2),'Ylim',[0 30]);

fvx = squeeze(sum(fxvxvz,2))*dv;
fvz = squeeze(sum(fxvxvz,3))*dv;

pcolor(h(3),xpos/1e3,vxvec/1e3,log10(fvx'))
shading(h(3),'flat')
axis(h(3),[-2000 500 -1300 700])
ylabel(h(3),'v_x (km s^{-1})','FontSize',12)
caxis(h(3),[-4 2]);
hold(h(3),'on')
plot(h(3),[xpos(idx1) xpos(idx1)]/1e3,[-800 800],'m--');
plot(h(3),[xpos(idx2) xpos(idx2)]/1e3,[-800 800],'m--');
hold(h(3),'off')
c = colorbar('peer',h(3));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(3),'XTickLabel',[])
colormap(h(3),'jet')

pcolor(h(4),xpos/1e3,vzvec/1e3,log10(fvz'))
shading(h(4),'flat')
axis(h(4),[-2000 500 -1000 1000])
ylabel(h(4),'v_z (km s^{-1})','FontSize',12)
caxis(h(4),[-4 2]);
hold(h(4),'on')
plot(h(4),[xpos(idx1) xpos(idx1)]/1e3,[-800 800],'m--');
plot(h(4),[xpos(idx2) xpos(idx2)]/1e3,[-800 800],'m--');
hold(h(4),'off')
c = colorbar('peer',h(4));
c.Label.String = 'log_{10}f_i (s m^{-4})';
xlabel(h(4),'x (km)')
colormap(h(4),'jet')


xcirc = -abs(Vsw)/1e3 + 2*abs(Vsw)*sind(1:361)/1e3;
ycirc = 2*abs(Vsw)*cosd(1:361)/1e3;

pcolor(h(5),vxvec/1e3,vzvec/1e3,log10(squeeze(fxvxvz(idx1,:,:))));
shading(h(5),'flat')
hold(h(5),'on')
plot(h(5),xcirc,ycirc,'r')
hold(h(5),'off')
xlabel(h(5),'v_x (km s^{-1})','FontSize',14)
ylabel(h(5),'v_z (km s^{-1})','FontSize',14)
c = colorbar('peer',h(5));
c.Label.String = 'f_i (s^2 m^{-5})';
colormap(h(5),'jet')
caxis(h(5),[-13 -2])
%caxis([1e-5 log10(max(max(fxz)))])

pcolor(h(6),vxvec/1e3,vzvec/1e3,log10(squeeze(fxvxvz(idx2,:,:))));
shading(h(6),'flat')
hold(h(6),'on')
plot(h(6),xcirc,ycirc,'r')
hold(h(6),'off')
xlabel(h(6),'v_x (km s^{-1})','FontSize',14)
ylabel(h(6),'v_z (km s^{-1})','FontSize',14)
c = colorbar('peer',h(6));
c.Label.String = 'f_i (s^2 m^{-5})';
colormap(h(6),'jet')
caxis(h(6),[-13 -2])
set(gcf,'color','w')
set(h(1:6),'fontsize',12)

end

%%
function dydt = temp(t,q)

l = 10e3;
vsw = -320e3;

n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 

P0 = 0.0425e-9;
P1 = 0.0575e-9;

mu0 = 1.2566e-06;
eomi = 9.5788e+07/4;
e = 1.6022e-19;

By = -B0*tanh(q(1)/l)+B1;
ni = -n0*tanh(q(1)/l)+n1;
Jz = -B0*sech(q(1)/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni) + P0*sech(q(1)/l).^2./(e*ni*l);
Vx = vsw*(n1-n0)./ni;
Ey = 0;
Ez = -Vx*By;

dydt = [q(2) ; eomi*Ex - eomi*q(6)*By; ...
        q(4) ; eomi*Ey; ...
        q(6) ; eomi*Ez + eomi*q(2)*By];
end