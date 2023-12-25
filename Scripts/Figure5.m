%% 
ic = 1;
Tint = irf.tint('2023-04-24T03:49:50.00Z/2023-04-24T03:50:20.00Z');
Tintshock = irf_time('2023-04-24T03:50:12.32Z','iso>epochtt');

ndsl = [0.9441 0.0997 -0.3142];
t1dsl = [-0.1524 -0.7125 -0.6849];
t2dsl = [-0.2920 0.6949 -0.6571];
ndsl = ndsl/norm(ndsl);
t1dsl = t1dsl/norm(t1dsl);
t2dsl = t2dsl/norm(t2dsl);

iPDist = mms.get_data('PDi_fpi_brst_l2',Tint,ic);
iPDist.data(:,1:15,:,:) = 0;

vg2D = linspace(-1400,1400,100);

vlimn = [-1400,800];
vlimt1 = [-600,400];
vlimt2 = [-800,1200];

%%
Tints1 = irf_time('2023-04-24T03:50:06.80Z','iso>epochtt');
Tints2 = irf_time('2023-04-24T03:50:08.00Z','iso>epochtt');
Tints3 = irf_time('2023-04-24T03:50:09.70Z','iso>epochtt');
Tints4 = irf_time('2023-04-24T03:50:11.00Z','iso>epochtt');
Tints5 = irf_time('2023-04-24T03:50:12.5Z','iso>epochtt');

[~,indT1] = min(abs(iPDist.time-Tints1));
[~,indT2] = min(abs(iPDist.time-Tints2));
[~,indT3] = min(abs(iPDist.time-Tints3));
[~,indT4] = min(abs(iPDist.time-Tints4));
[~,indT5] = min(abs(iPDist.time-Tints5));

nMC = 5e3;
f2Dnt21 = iPDist(indT1).reduce('2D',ndsl,t2dsl,'base','cart','vg',vg2D,'nMC',nMC);
f2Dnt22 = iPDist(indT2).reduce('2D',ndsl,t2dsl,'base','cart','vg',vg2D,'nMC',nMC);
f2Dnt23 = iPDist(indT3).reduce('2D',ndsl,t2dsl,'base','cart','vg',vg2D,'nMC',nMC);
f2Dnt24 = iPDist(indT4).reduce('2D',ndsl,t2dsl,'base','cart','vg',vg2D,'nMC',nMC);
f2Dnt25 = iPDist(indT5).reduce('2D',ndsl,t2dsl,'base','cart','vg',vg2D,'nMC',nMC);

%% Model ion distributions
load('fdistproton.mat')
load('fdistalpha.mat')
load('fdisthelium.mat')

xposp = fdistproton.xpos;
fxvxvzp = fdistproton.fxvxvz*8/9;
vxvec = fdistproton.vxvec;
vzvec = fdistproton.vzvec;
Emultiplier = fdistproton.Emultiplier;
dv = fdistproton.dv;
vxmat = ones(size(xposp'))*vxvec;
vzmat = ones(size(xposp'))*vzvec;

fvxp = squeeze(sum(fxvxvzp,2))*dv;
fvzp = squeeze(sum(fxvxvzp,3))*dv;

% Alphas
fxvxvza = fdistalpha.fxvxvz*0.1*8/9;
fvxa = squeeze(sum(fxvxvza,2))*dv;
fvza = squeeze(sum(fxvxvza,3))*dv;

% Helium
fxvxvzhe = fdisthelium.fxvxvz*0.01*8/9;
fvxhe = squeeze(sum(fxvxvzhe,2))*dv;
fvzhe = squeeze(sum(fxvxvzhe,3))*dv;

%Convert to spacecraft frame
vxmatp = vxmat-130e3;
vzmatp = vzmat+240e3;

vxmata = (vxmat-130e3)*sqrt(2);
vzmata = (vzmat+240e3)*sqrt(2);

vxmathe = (vxmat-130e3)*2;
vzmathe = (vzmat+240e3)*2;

fvxhesc = interp2(xposp,vxmathe(1,:)',fvxhe',xposp,vxmatp(1,:)');
fvxasc = interp2(xposp,vxmata(1,:)',fvxa',xposp,vxmatp(1,:)');

fallvx = fvxp'+fvxasc+fvxhesc;

fvzhesc = interp2(xposp,vzmathe(1,:)',fvzhe',xposp,vzmatp(1,:)');
fvzasc = interp2(xposp,vzmata(1,:)',fvza',xposp,vzmatp(1,:)');

fallvz = fvzp'+fvzasc+fvzhesc;


post1 = 30e3;
post2 = -190e3;
post3 = -390e3;
post4 = -560e3;
post5 = -670e3;
[~,idx1] = min(abs(post1 - xposp));
[~,idx2] = min(abs(post2 - xposp));
[~,idx3] = min(abs(post3 - xposp));
[~,idx4] = min(abs(post4 - xposp));
[~,idx5] = min(abs(post5 - xposp));

distp1 = squeeze(fxvxvzp(idx1,:,:));
distp2 = squeeze(fxvxvzp(idx2,:,:));
distp3 = squeeze(fxvxvzp(idx3,:,:));
distp4 = squeeze(fxvxvzp(idx4,:,:));
distp5 = squeeze(fxvxvzp(idx5,:,:));

disthe1 = squeeze(fxvxvzhe(idx1,:,:));
disthe2 = squeeze(fxvxvzhe(idx2,:,:));
disthe3 = squeeze(fxvxvzhe(idx3,:,:));
disthe4 = squeeze(fxvxvzhe(idx4,:,:));
disthe5 = squeeze(fxvxvzhe(idx5,:,:));

dista1 = squeeze(fxvxvza(idx1,:,:));
dista2 = squeeze(fxvxvza(idx2,:,:));
dista3 = squeeze(fxvxvza(idx3,:,:));
dista4 = squeeze(fxvxvza(idx4,:,:));
dista5 = squeeze(fxvxvza(idx5,:,:));

disthe1 = interp2((vxvec-130e3)*2,(vzvec'+240e3)*2,disthe1,(vxvec-130e3),(vzvec'+240e3));
disthe2 = interp2((vxvec-130e3)*2,(vzvec'+240e3)*2,disthe2,(vxvec-130e3),(vzvec'+240e3));
disthe3 = interp2((vxvec-130e3)*2,(vzvec'+240e3)*2,disthe3,(vxvec-130e3),(vzvec'+240e3));
disthe4 = interp2((vxvec-130e3)*2,(vzvec'+240e3)*2,disthe4,(vxvec-130e3),(vzvec'+240e3));
disthe5 = interp2((vxvec-130e3)*2,(vzvec'+240e3)*2,disthe5,(vxvec-130e3),(vzvec'+240e3));

dista1 = interp2((vxvec-130e3)*sqrt(2),(vzvec'+240e3)*sqrt(2),dista1,(vxvec-130e3),(vzvec'+240e3));
dista2 = interp2((vxvec-130e3)*sqrt(2),(vzvec'+240e3)*sqrt(2),dista2,(vxvec-130e3),(vzvec'+240e3));
dista3 = interp2((vxvec-130e3)*sqrt(2),(vzvec'+240e3)*sqrt(2),dista3,(vxvec-130e3),(vzvec'+240e3));
dista4 = interp2((vxvec-130e3)*sqrt(2),(vzvec'+240e3)*sqrt(2),dista4,(vxvec-130e3),(vzvec'+240e3));
dista5 = interp2((vxvec-130e3)*sqrt(2),(vzvec'+240e3)*sqrt(2),dista5,(vxvec-130e3),(vzvec'+240e3));


vxp = vxvec - 130e3;
vzp = vzvec + 240e3;


distall1 = distp1+dista1+disthe1;
distall2 = distp2+dista2+disthe2;
distall3 = distp3+dista3+disthe3;
distall4 = distp4+dista4+disthe4;
distall5 = distp5+dista5+disthe5;

%%
h=irf_plot(10,'newfigure');
%h=irf_figure(540+ic,8);
xSize=1000; ySize=450;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.89;
ywidth = 0.15;
ywidth2 = 0.1;

set(h(1),'position',[0.015 0.59 0.27 0.37]);
set(h(2),'position',[0.185 0.59 0.27 0.37]);
set(h(3),'position',[0.355 0.59 0.27 0.37]);
set(h(4),'position',[0.525 0.59 0.27 0.37]);
set(h(5),'position',[0.718 0.59 0.27 0.37]);

set(h(6),'position',[0.015 0.11 0.27 0.37]);
set(h(7),'position',[0.185 0.11 0.27 0.37]);
set(h(8),'position',[0.355 0.11 0.27 0.37]);
set(h(9),'position',[0.525 0.11 0.27 0.37]);
set(h(10),'position',[0.718 0.11 0.27 0.37]);

c = [55,137,187;...
  106,193,165;...
  172,220,166;...
  230,244,157;...
  255,254,194;...
  253,223,144;...
  251,173,104;...
  242,109,074;...
  211,064,082]/255;
cmap = interp1(linspace(1,64,size(c,1)),c,1:64);

Vswscn = -450;
Vswsct2 = 240;
Vswsh = 330;

xcirc = -abs(Vswscn) + 2*abs(Vswsh)*sind(1:361);
ycirc = Vswsct2 + 2*abs(Vswsh)*cosd(1:361);

vlimt2 = [-900,1300];

f2Dnt21.plot_plane(h(1),'colorbar',1,'contour',6)
colormap(h(1),'jet');
xlabel(h(1),'v_{n} (km s^{-1})','fontsize',14)
ylabel(h(1),'v_{t2} (km s^{-1})','fontsize',14)
hold(h(1),'on')
plot(h(1),xcirc,ycirc,'b','linewidth',1.5)
plot(h(1),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(1),'off')
c=colorbar('peer',h(1));
delete(c);
set(h(1),'xgrid','off','ygrid','off');
caxis(h(1),[-9 -3]);
axis(h(1),'square')
irf_legend(h(1),'(a)',[0.99 0.99],'color','r','FontSize',14);
colormap(h(1),cmap)
axis(h(1),[vlimn vlimt2])
tinttemp = irf_time(f2Dnt21.time,'epochtt>utc');
title(h(1),[tinttemp(12:21) ' UTC']);

f2Dnt22.plot_plane(h(2),'colorbar',1,'contour',6)
colormap(h(2),'jet');
xlabel(h(2),'v_{n} (km s^{-1})','fontsize',14)
ylabel(h(2),' ','fontsize',12)
c=colorbar('peer',h(2));
delete(c);
hold(h(2),'on')
plot(h(2),xcirc,ycirc,'b','linewidth',1.5)
plot(h(2),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(2),'off')
set(h(2),'xgrid','off','ygrid','off');
set(h(2),'YTickLabel',[])
caxis(h(2),[-9 -3]);
axis(h(2),'square')
irf_legend(h(2),'(b)',[0.99 0.99],'color','r','FontSize',14);
colormap(h(2),cmap)
axis(h(2),[vlimn vlimt2])
tinttemp = irf_time(f2Dnt22.time,'epochtt>utc');
title(h(2),[tinttemp(12:21) ' UTC']);

f2Dnt23.plot_plane(h(3),'colorbar',1,'contour',6)
colormap(h(3),'jet');
xlabel(h(3),'v_{n} (km s^{-1})','fontsize',14)
ylabel(h(3),' ','fontsize',12)
c=colorbar('peer',h(3));
delete(c);
hold(h(3),'on')
plot(h(3),xcirc,ycirc,'b','linewidth',1.5)
plot(h(3),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(3),'off')
set(h(3),'xgrid','off','ygrid','off');
set(h(3),'YTickLabel',[])
caxis(h(3),[-9 -3]);
axis(h(3),'square')
irf_legend(h(3),'(c)',[0.99 0.99],'color','r','FontSize',14);
colormap(h(3),cmap)
axis(h(3),[vlimn vlimt2])
tinttemp = irf_time(f2Dnt23.time,'epochtt>utc');
title(h(3),[tinttemp(12:21) ' UTC']);

f2Dnt24.plot_plane(h(4),'colorbar',1,'contour',6)
colormap(h(4),'jet');
xlabel(h(4),'v_{n} (km s^{-1})','fontsize',14)
ylabel(h(4),' ','fontsize',12)
c=colorbar('peer',h(4));
delete(c);
hold(h(4),'on')
plot(h(4),xcirc,ycirc,'b','linewidth',1.5)
plot(h(4),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(4),'off')
set(h(4),'xgrid','off','ygrid','off');
set(h(4),'YTickLabel',[])
caxis(h(4),[-9 -3]);
axis(h(4),'square')
irf_legend(h(4),'(d)',[0.99 0.99],'color','r','FontSize',14);
colormap(h(4),cmap)
axis(h(4),[vlimn vlimt2])
tinttemp = irf_time(f2Dnt24.time,'epochtt>utc');
title(h(4),[tinttemp(12:21) ' UTC']);

f2Dnt25.units = 's^2 m^{-5}';
f2Dnt25.plot_plane(h(5),'colorbar',1,'contour',6)
colormap(h(5),'jet');
xlabel(h(5),'v_{n} (km s^{-1})','fontsize',14)
ylabel(h(5),' ','fontsize',12)
c=colorbar('peer',h(5));
hold(h(5),'on')
plot(h(5),xcirc,ycirc,'b','linewidth',1.5)
plot(h(5),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(5),'off')
set(h(5),'xgrid','off','ygrid','off');
set(h(5),'YTickLabel',[])
caxis(h(5),[-9 -3]);
axis(h(5),'square')
irf_legend(h(5),'(e)',[0.99 0.99],'color','r','FontSize',14);
ylabel(c,'f_i (s^2 m^{-5})','FontSize',14);
colormap(h(5),cmap)
axis(h(5),[vlimn vlimt2])
irf_legend(h(5),'H^{+}',[0.5 0.35],'fontsize',12,'color','r')
irf_legend(h(5),'He^{2+}',[0.19 0.36],'fontsize',12,'color','r')
irf_legend(h(5),'He^{+}',[0.15 0.85],'fontsize',12,'color','r')
irf_legend(h(5),'H^{+}',[0.75 0.85],'fontsize',12,'color','r')
tinttemp = irf_time(f2Dnt25.time,'epochtt>utc');
title(h(5),[tinttemp(12:21) ' UTC']);

% plot model distributions

Vswscn = -450;
Vswsct2 = 240;
Vswsh = 330;

xcirc = -abs(Vswscn) + 2*abs(Vswsh)*sind(1:361);
ycirc = Vswsct2 + 2*abs(Vswsh)*cosd(1:361);

vlimt2 = [-900,1300];
vlimn = [-1400,800];

ftemp = distall5;
ftemp(ftemp < 1e-9) = NaN;
pcolor(h(6),vxp/1e3,vzp/1e3,log10(ftemp));
shading(h(6),'flat')
hold(h(6),'on')
plot(h(6),xcirc,ycirc,'b','linewidth',1.5)
plot(h(6),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(6),'off')
xlabel(h(6),'v_{n} (km s^{-1})','FontSize',14)
ylabel(h(6),'v_{t2} (km s^{-1})','FontSize',14)
colormap(h(6),cmap)
caxis(h(6),[-9 -3])
axis(h(6),'square')
axis(h(6),[vlimn vlimt2])
irf_legend(h(6),'(f)',[0.98 0.98],'fontsize',14,'color','r')
irf_legend(h(6),'n = -670 km',[0.04 0.03],'fontsize',12,'color','k')

ftemp = distall4;
ftemp(ftemp < 1e-9) = NaN;
pcolor(h(7),vxp/1e3,vzp/1e3,log10(ftemp));
shading(h(7),'flat')
hold(h(7),'on')
plot(h(7),xcirc,ycirc,'b','linewidth',1.5)
plot(h(7),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(7),'off')
xlabel(h(7),'v_{n} (km s^{-1})','FontSize',14)
ylabel(h(7),' ','FontSize',14)
set(h(7),'YTickLabel',[])
colormap(h(7),cmap)
caxis(h(7),[-9 -3])
axis(h(7),'square')
axis(h(7),[vlimn vlimt2])
irf_legend(h(7),'(g)',[0.98 0.98],'fontsize',14,'color','r')
irf_legend(h(7),'n = -560 km',[0.04 0.03],'fontsize',12,'color','k')

ftemp = distall3;
ftemp(ftemp < 1e-9) = NaN;
pcolor(h(8),vxp/1e3,vzp/1e3,log10(ftemp));
shading(h(8),'flat')
hold(h(8),'on')
plot(h(8),xcirc,ycirc,'b','linewidth',1.5)
plot(h(8),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(8),'off')
xlabel(h(8),'v_{n} (km s^{-1})','FontSize',14)
ylabel(h(8),' ','FontSize',14)
set(h(8),'YTickLabel',[])
colormap(h(8),cmap)
caxis(h(8),[-9 -3])
axis(h(8),'square')
axis(h(8),[vlimn vlimt2])
irf_legend(h(8),'(h)',[0.98 0.98],'fontsize',14,'color','r')
irf_legend(h(8),'n = -390 km',[0.04 0.03],'fontsize',12,'color','k')

ftemp = distall2;
ftemp(ftemp < 1e-9) = NaN;
pcolor(h(9),vxp/1e3,vzp/1e3,log10(ftemp));
shading(h(9),'flat')
hold(h(9),'on')
plot(h(9),xcirc,ycirc,'b','linewidth',1.5)
plot(h(9),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(9),'off')
xlabel(h(9),'v_{n} (km s^{-1})','FontSize',14)
ylabel(h(9),' ','FontSize',14)
set(h(9),'YTickLabel',[])
colormap(h(9),cmap)
caxis(h(9),[-9 -3])
axis(h(9),'square')
axis(h(9),[vlimn vlimt2])
irf_legend(h(9),'(i)',[0.98 0.98],'fontsize',14,'color','r')
irf_legend(h(9),'n = -190 km',[0.04 0.03],'fontsize',12,'color','k')

ftemp = distall1;
ftemp(ftemp < 1e-9) = NaN;
pcolor(h(10),vxp/1e3,vzp/1e3,log10(ftemp));
shading(h(10),'flat')
hold(h(10),'on')
plot(h(10),xcirc,ycirc,'b','linewidth',1.5)
plot(h(10),Vswscn,Vswsct2,'b.','linewidth',2,'MarkerSize',14)
hold(h(10),'off')
xlabel(h(10),'v_{n} (km s^{-1})','FontSize',14)
ylabel(h(10),' ','FontSize',14)
set(h(10),'YTickLabel',[])
colormap(h(10),cmap)
caxis(h(10),[-9 -3])
axis(h(10),'square')
axis(h(10),[vlimn vlimt2])
c = colorbar('peer',h(10));
irf_legend(h(10),'H^{+}',[0.5 0.40],'fontsize',12,'color','r')
irf_legend(h(10),'He^{2+}',[0.23 0.53],'fontsize',12,'color','r')
irf_legend(h(10),'He^{+}',[0.18 0.80],'fontsize',12,'color','r')
irf_legend(h(10),'H^{+}',[0.75 0.80],'fontsize',12,'color','r')
c.Label.String = 'f_i (s^2 m^{-5})';
irf_legend(h(10),'(j)',[0.98 0.98],'fontsize',14,'color','r')
irf_legend(h(10),'n = 30 km',[0.04 0.03],'fontsize',12,'color','k')