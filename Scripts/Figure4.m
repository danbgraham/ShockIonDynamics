%% Bowshock model
Vsw = -320e3; 
n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;

l = 10e3; 

xvec = [-20e2:1:2e2]*1e3;
dx = median(diff(xvec));

By = -B0*tanh(xvec/l)+B1;
ni = -n0*tanh(xvec/l)+n1;

Units = irf_units;
mu0 = Units.mu0;
e = Units.e;
mi = Units.mp;

Jz = -B0*sech(xvec/l).^2/(mu0*l);
Vx = Vsw*nsw./ni;

Ex = -Jz.*By./(e*ni)+P0*sech(xvec/l).^2./(e*ni*l);
Ez = -Vx.*By;

phi = -cumsum(Ex)*dx;

%% 
load('fdistproton.mat');
load('fdisthelium.mat');
load('fdistalpha.mat');

% Protons
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
%% Compute moments of all dists in SC frame
npsc = squeeze(sum(fallvx',2))*dv;

vxsc = sum(fallvx'.*vxmatp*dv,2)./npsc;
vzsc = sum(fallvz'.*vzmatp*dv,2)./npsc;

vxscmat = vxsc*ones(size(vxvec));
vzscmat = vzsc*ones(size(vzvec));

Txxsc = sum(fallvx'.*(vxmatp - vxscmat).^2*dv,2)./npsc/e*mi;
Tzzsc = sum(fallvz'.*(vzmatp - vzscmat).^2*dv,2)./npsc/e*mi;

Tscs = (Txxsc + Tzzsc + 6)/3;

%% Compute moments
% Protons
np = squeeze(sum(sum(fxvxvzp,3),2))*dv^2;
vxp = sum(fvxp.*vxmat*dv,2)./np;
vzp = sum(fvzp.*vzmat*dv,2)./np;

vxpmat = vxp*ones(size(vxvec));
vzpmat = vzp*ones(size(vzvec));

Txxp = sum(fvxp.*(vxmat - vxpmat).^2*dv,2)./np/e*mi;
Tzzp = sum(fvzp.*(vzmat - vzpmat).^2*dv,2)./np/e*mi;

np = np*1e-6;
vxp = vxp*1e-3;
vzp = vzp*1e-3;

Tps = (Txxp + Tzzp + 3)/3;

% Alpha particles
na = squeeze(sum(sum(fxvxvza,3),2))*dv^2;
vxa = sum(fvxa.*vxmat*dv,2)./na;
vza = sum(fvza.*vzmat*dv,2)./na;

vxamat = vxa*ones(size(vxvec));
vzamat = vza*ones(size(vzvec));

Txxa = sum(fvxa.*(vxmat - vxamat).^2*dv,2)./na/e*4*mi;
Tzza = sum(fvza.*(vzmat - vzamat).^2*dv,2)./na/e*4*mi;

na = na*1e-6;
vxa = vxa*1e-3;
vza = vza*1e-3;

Tas = (Txxa + Tzza + 12)/3;

% Helium ions
nhe = squeeze(sum(sum(fxvxvzhe,3),2))*dv^2;
vxhe = sum(fvxhe.*vxmat*dv,2)./nhe;
vzhe = sum(fvzhe.*vzmat*dv,2)./nhe;

vxhemat = vxhe*ones(size(vxvec));
vzhemat = vzhe*ones(size(vzvec));

Txxhe = sum(fvxhe.*(vxmat - vxhemat).^2*dv,2)./nhe/e*4*mi;
Tzzhe = sum(fvzhe.*(vzmat - vzhemat).^2*dv,2)./nhe/e*4*mi;

nhe = nhe*1e-6;
vxhe = vxhe*1e-3;
vzhe = vzhe*1e-3;

Thes = (Txxhe + Tzzhe + 12)/3;


%% Load data
Tint = irf.tint('2023-04-24T03:49:35.00Z/2023-04-24T03:50:30.00Z');

Bxyz=mms.get_data('B_dmpa_brst_l2',Tint,1);
iPDist = mms.get_data('PDi_fpi_brst_l2',Tint,1);
c_eval('SCpot=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',Tint);',1);

ndsl = [0.9441 0.0997 -0.3142];
t1dsl = [-0.1524 -0.7125 -0.6849];
t2dsl = [-0.2920 0.6949 -0.6571];
ndsl = ndsl/norm(ndsl);
t1dsl = t1dsl/norm(t1dsl);
t2dsl = t2dsl/norm(t2dsl);

Bntt = irf_newxyz(Bxyz,ndsl,t1dsl,t2dsl);

iPDist.data(:,1:15,:,:) = 0;

vlimn = [-1400,800];
vlimt1 = [-600,800];
vlimt2 = [-800,1200];

vg1Dn = linspace(vlimn(1),vlimn(2),150);
vg1Dt1 = linspace(vlimt1(1),vlimt1(2),150);
vg1Dt2 = linspace(vlimt2(1),vlimt2(2),150);
nMC = 2e2;

nvec = irf.ts_vec_xyz(iPDist.time,repmat(ndsl,length(iPDist.time),1));
t1vec = irf.ts_vec_xyz(iPDist.time,repmat(t1dsl,length(iPDist.time),1));
t2vec = irf.ts_vec_xyz(iPDist.time,repmat(t2dsl,length(iPDist.time),1));

f1Dn = iPDist.reduce('1D',nvec,'vg',vg1Dn,'nMC',nMC);
f1Dt1 = iPDist.reduce('1D',t1vec,'vg',vg1Dt1,'nMC',nMC);
f1Dt2 = iPDist.reduce('1D',t2vec,'vg',vg1Dt2,'nMC',nMC);

Tintshock = irf_time('2023-04-24T03:50:12.32Z','iso>epochtt');
Vshock = 100e3;
Tints = Tintshock + [-2000 200]/100;

%% Plot Figure

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

xwidth = 0.41;
ywidth = 0.13;

xwidth2 = 0.41;
ywidth2 = 0.17;

fn=figure;
set(fn,'Position',[10 10 1000 700])
h(1)=axes('position',[0.57 0.99-ywidth2 xwidth2+0.003 ywidth2]);
h(2)=axes('position',[0.57 0.99-ywidth2*2 xwidth2+0.003 ywidth2]);
h(3)=axes('position',[0.57 0.99-ywidth2*3 xwidth2+0.003 ywidth2]);

h(4)=axes('position',[0.065 0.99-ywidth xwidth-0.044 ywidth]);
h(5)=axes('position',[0.065 0.99-2*ywidth-0.002 xwidth ywidth]);
h(6)=axes('position',[0.065 0.99-3*ywidth-0.002*2 xwidth ywidth]);
h(7)=axes('position',[0.065 0.99-4*ywidth-0.002*3 xwidth ywidth]);
h(8)=axes('position',[0.065 0.99-5*ywidth-0.002*4 xwidth ywidth]); % [x y dx dy]
h(9)=axes('position',[0.065 0.99-6*ywidth-0.002*5 xwidth ywidth]);
h(10)=axes('position',[0.065 0.99-7*ywidth-0.002*6 xwidth ywidth]);

h(11)=axes('position',[0.57 0.92-ywidth2*4 xwidth2 ywidth2]);
h(12)=axes('position',[0.57 0.92-ywidth2*5-0.002 xwidth2 ywidth2]);

set(fn,'defaultLineLineWidth',2);
set(fn,'defaultAxesFontSize',14);

%h(1)=irf_panel('Bntt');
irf_plot(h(1),Bntt.x,'k');
hold(h(1),'on')
irf_plot(h(1),Bntt.y,'b');
irf_plot(h(1),Bntt.z,'r');
hold(h(1),'off')
ylabel(h(1),{'B (nT)'},'Interpreter','tex');
irf_zoom(h(1),'y',[-10 65])
irf_legend(h(1),{'B_{n}'},[0.03 0.3],'fontsize',14,'color','k')
irf_legend(h(1),{'B_{t1}'},[0.1 0.3],'fontsize',14,'color','b')
irf_legend(h(1),{'B_{t2}'},[0.17 0.3],'fontsize',14,'color','r')
irf_legend(h(1),'(h)',[0.99 0.94],'color','k','fontsize',14)


tion1 = irf_time([2023 04 24 03 50 06.801])+[-0.075 0.075];
tion2 = irf_time([2023 04 24 03 50 08.001])+[-0.075 0.075];
tion3 = irf_time([2023 04 24 03 50 09.651])+[-0.075 0.075];
tion4 = irf_time([2023 04 24 03 50 11.001])+[-0.075 0.075];
tion5 = irf_time([2023 04 24 03 50 12.501])+[-0.075 0.075];
yi = [650 650];

specrecn = f1Dn.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i} (s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
irf_spectrogram(h(2),specrecn);
hold(h(2),'on')
irf_plot(h(2),[tion1' yi'],'color','m','linewidth',12)
irf_plot(h(2),[tion2' yi'],'color','m','linewidth',12)
irf_plot(h(2),[tion3' yi'],'color','m','linewidth',12)
irf_plot(h(2),[tion4' yi'],'color','m','linewidth',12)
irf_plot(h(2),[tion5' yi'],'color','m','linewidth',12)
hold(h(2),'off')
ylabel(h(2),'v_n (km s^{-1})','interpreter','tex')
colormap(h(2),cmap)
irf_zoom(h(2),'y',vlimn)
caxis(h(2),[-3 2])
irf_legend(h(2),'(i)',[0.99 0.94],'color','k','fontsize',14)

specrecn = f1Dt2.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i} (s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
irf_spectrogram(h(3),specrecn);
ylabel(h(3),'v_{t2} (km s^{-1})','interpreter','tex')
colormap(h(3),cmap)
irf_zoom(h(3),'y',vlimt2)
caxis(h(3),[-3 2])
irf_legend(h(3),'(j)',[0.99 0.94],'color','k','fontsize',14)


set(h(2:3),'xgrid','off','ygrid','off')
set(h(2:3),'Color',0.85*[1 1 1]);

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tints);

plot(h(4),xvec/1e3,By*1e9,'k')
ylabel(h(4),{'B (nT)'},'fontsize',14)
xlabel(h(4),'n (10^3 km)','fontsize',14)
axis(h(4),[xvec(1)/1e3 xvec(end)/1e3 0 65])
set(h(4),'XTickLabel',[])
irf_legend(h(4),'(a)',[0.98 0.96],'fontsize',14)
ax2 = axes('Position',get(h(4),'Position'),...
		'XAxisLocation','top',...
		'YAxisLocation','right',...
		'Color','none',...
		'XColor','k','YColor','r',...
		'XTickLabel',[],'Xtick',[]);
axes(ax2)
line(xvec/1e3,phi,'Color','r','Parent',ax2);
set(ax2,'Ylim',[-550 100]);
set(ax2,'Xlim',[xvec(1) xvec(end)]/1e3);
ylabel('\phi (V)','FontSize',14)

ftemp = fvxp;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(5),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(5),'flat')
axis(h(5),[-2000 200 -900 700])
ylabel(h(5),'v_n (km s^{-1})','FontSize',14)
caxis(h(5),[-2.99 2]);
c = colorbar('peer',h(5));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(5),'XTickLabel',[])
colormap(h(5),cmap)
irf_legend(h(5),'(b)',[0.98 0.96],'fontsize',14)
irf_legend(h(5),'H^{+}',[0.02 0.96],'fontsize',14)

ftemp = fvxa;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(6),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(6),'flat')
axis(h(6),[-2000 200 -900 700])
ylabel(h(6),'v_n (km s^{-1})','FontSize',14)
caxis(h(6),[-2.99 2]);
c = colorbar('peer',h(6));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(6),'XTickLabel',[])
colormap(h(6),cmap)
irf_legend(h(6),'(c)',[0.98 0.96],'fontsize',14)
irf_legend(h(6),'He^{2+}',[0.02 0.96],'fontsize',14)

ftemp = fvxhe;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(7),xposp/1e3,vxvec/1e3,log10(ftemp'))
shading(h(7),'flat')
axis(h(7),[-2000 200 -900 700])
ylabel(h(7),'v_n (km s^{-1})','FontSize',14)
caxis(h(7),[-2.99 2]);
c = colorbar('peer',h(7));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(7),'XTickLabel',[])
colormap(h(7),cmap)
irf_legend(h(7),'(d)',[0.98 0.96],'fontsize',14)
irf_legend(h(7),'He^{+}',[0.02 0.96],'fontsize',14)

ftemp = fvzp;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(8),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(8),'flat')
axis(h(8),[-2000 200 -800 800])
ylabel(h(8),'v_{t2} (km s^{-1})','FontSize',14)
caxis(h(8),[-2.99 2]);
c = colorbar('peer',h(8));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(8),'XTickLabel',[])
colormap(h(8),cmap)
irf_legend(h(8),'(e)',[0.98 0.96],'fontsize',14)
irf_legend(h(8),'H^{+}',[0.02 0.96],'fontsize',14)

ftemp = fvza;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(9),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(9),'flat')
axis(h(9),[-2000 200 -800 800])
ylabel(h(9),'v_{t2} (km s^{-1})','FontSize',14)
caxis(h(9),[-2.99 2]);
c = colorbar('peer',h(9));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(9),'XTickLabel',[])
colormap(h(9),cmap)
irf_legend(h(9),'(f)',[0.98 0.96],'fontsize',14)
irf_legend(h(9),'He^{2+}',[0.02 0.96],'fontsize',14)

ftemp = fvzhe;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(10),xposp/1e3,vzvec/1e3,log10(ftemp'))
shading(h(10),'flat')
axis(h(10),[-2000 200 -800 800])
ylabel(h(10),'v_{t2} (km s^{-1})','FontSize',14)
caxis(h(10),[-2.99 2]);
c = colorbar('peer',h(10));
c.Label.String = 'log_{10}f_i (s m^{-4})';
%set(h(10),'XTickLabel',[])
irf_legend(h(10),'(g)',[0.98 0.96],'fontsize',14)
irf_legend(h(10),'He^{+}',[0.02 0.96],'fontsize',14)
colormap(h(10),cmap)
xlabel(h(10),'n (km)')

set(h(5:10),'xgrid','off','ygrid','off')
set(h(5:10),'Color',0.85*[1 1 1]);

tion1 = -670+[0 10];
tion2 = -560+[0 10];
tion3 = -390+[0 10];
tion4 = -190+[0 10];
tion5 = 30+[0 10];
yi = [650 650];

ftemp = fallvx;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(11),xposp/1e3,(vxvec-130e3)/1e3,log10(ftemp))
shading(h(11),'flat')
axis(h(11),[-2000 200 -1400 800])
hold(h(11),'on')
plot(h(11),tion1',yi','color','m','linewidth',12)
plot(h(11),tion2',yi','color','m','linewidth',12)
plot(h(11),tion3',yi','color','m','linewidth',12)
plot(h(11),tion4',yi','color','m','linewidth',12)
plot(h(11),tion5',yi','color','m','linewidth',12)
hold(h(11),'off')
ylabel(h(11),'v_{n} (km s^{-1})','FontSize',14)
caxis(h(11),[-2.99 2]);
c = colorbar('peer',h(11));
c.Label.String = 'log_{10}f_i (s m^{-4})';
set(h(11),'XTickLabel',[])
irf_legend(h(11),'(k)',[0.98 0.96],'fontsize',14)
irf_legend(h(11),'All ions',[0.02 0.96],'fontsize',14)
title(h(11),'')
colormap(h(11),cmap)
%xlabel(h(11),'n (km)')

ftemp = fallvz;
ftemp(ftemp < 1e-3) = NaN;
pcolor(h(12),xposp/1e3,(vzvec+240e3)/1e3,log10(ftemp))
shading(h(12),'flat')
axis(h(12),[-2000 200 -800 1200])
ylabel(h(12),'v_{t2} (km s^{-1})','FontSize',14)
caxis(h(12),[-2.99 2]);
c = colorbar('peer',h(12));
c.Label.String = 'log_{10}f_i (s m^{-4})';
%set(h(12),'XTickLabel',[])
irf_legend(h(12),'(l)',[0.98 0.96],'fontsize',14)
colormap(h(12),cmap)
xlabel(h(12),'n (km)')

set(h(11:12),'xgrid','off','ygrid','off')
set(h(11:12),'Color',0.85*[1 1 1]);

set(h(1:12),'fontsize',14);

set(gcf,'color','w')