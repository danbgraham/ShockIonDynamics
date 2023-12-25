Vsw = -320e3; 
n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;


l = 10e3; 

xvec = [-1e3:1:1e3]*1e3;
dx = median(diff(xvec));
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

Ex = -Jz.*By./(e*ni)+P0*sech(xvec/l).^2./(e*ni*l);

Ez = -Vx.*By;

phi = cumsum(Ex)*dx;

opts = odeset('RelTol',1e-7,'AbsTol',1e-9);

  
x0 = 5000e3;
y0 = 0;
z0 = 0;
vx0 = -320e3;
vy0 = 0;
vz0 = 0e3;

[t1p,d1p] = ode45(@tempp,[0 50],[x0; vx0; y0; vy0; z0; vz0],opts);
[t1a,d1a] = ode45(@tempa,[0 50],[x0; vx0; y0; vy0; z0; vz0],opts);
[t1he,d1he] = ode45(@temphe,[0 50],[x0; vx0; y0; vy0; z0; vz0],opts);
%[t1a,d1a] = ode45(@tempa,[0 30],[x0; vx0; y0; vy0; z0; vz0],opts);
%[t1he,d1he] = ode45(@temphe,[0 30],[x0; vx0; y0; vy0; z0; vz0],opts);

protonmotion = struct('time',t1p,'posvel',d1p);
alphamotion = struct('time',t1a,'posvel',d1a);
heliummotion = struct('time',t1he,'posvel',d1he);

%% Load MMS data
ic = 1;
Tint = irf.tint('2023-04-24T03:49:50.00Z/2023-04-24T03:50:20.00Z');
Tintshock = irf_time('2023-04-24T03:50:12.32Z','iso>epochtt');
Tints = Tintshock + [-15.385 1.5385];


ndsl = [0.9441 0.0997 -0.3142];
t1dsl = [-0.1524 -0.7125 -0.6849];
t2dsl = [-0.2920 0.6949 -0.6571];
ndsl = ndsl/norm(ndsl);
t1dsl = t1dsl/norm(t1dsl);
t2dsl = t2dsl/norm(t2dsl);

c_eval('Bxyz=mms.get_data(''B_dmpa_brst_l2'',Tint,?);',ic);

Bntt = irf_newxyz(Bxyz,ndsl,t1dsl,t2dsl);
iPDist = mms.get_data('PDi_fpi_brst_l2',Tint,ic);
iPDist.data(:,1:15,:,:) = 0;

vlimn = [-1400,800];
vlimt1 = [-600,400];
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

Vshock = 130e3;

%% Plot example figure; 4 panels

fn=figure;
set(fn,'Position',[10 10 700 700])
h(1)=axes('position',[0.09 0.87 0.90 0.12]);
h(2)=axes('position',[0.08 0.68 0.405 0.12]);
h(3)=axes('position',[0.585 0.68 0.405 0.12]);
h(4)=axes('position',[0.08 0.54 0.405 0.12]);
h(5)=axes('position',[0.585 0.54 0.405 0.12]);
h(6)=axes('position',[0.10 0.33 0.89 0.13]);
h(7)=axes('position',[0.10 0.20 0.89 0.13]);
h(8)=axes('position',[0.10 0.07 0.89 0.13]);
set(fn,'defaultLineLineWidth',2);
set(fn,'defaultAxesFontSize',14)

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

plot(h(1),xvec/1e6,By*1e9)
hold(h(1),'on');
plot(h(1),xvec/1e6,Ex*1e3)
plot(h(1),xvec/1e6,Ez*1e3)
hold(h(1),'off');
ylabel(h(1),{'B (nT)','E (mV m^{-1})'},'fontsize',14)
xlabel(h(1),'n (10^3 km)','fontsize',14)
axis(h(1),[-1 1 0 65])
yticks(h(1),[0 20 40 60])
legend(h(1),{'B_{t1}','E_{n}','E_{t2}'},'Location','northeast','fontsize',14)
irf_legend(h(1),'(a)',[0.88 0.96],'fontsize',14)

plot(h(2),d1p(:,1)/1e6,d1p(:,2)/1e3)
hold(h(2),'on');
plot(h(2),d1a(:,1)/1e6,d1a(:,2)/1e3)
plot(h(2),d1he(:,1)/1e6,d1he(:,2)/1e3)
hold(h(2),'off');
%xlabel(h(2),'n (10^3 km)','fontsize',14)
set(h(2),'XTickLabel',[])
ylabel(h(2),'v_n (km s^{-1})','fontsize',14)
axis(h(2),[-4 2 -400 0])
legend(h(2),{'H^{+}','He^{2+}','He^{+}'},'Location','northeast','fontsize',14)
irf_legend(h(2),'(b)',[0.65 0.06],'fontsize',14)
title(h(2),'Shock frame')

plot(h(3),d1p(:,1)/1e6,(d1p(:,2)-130e3)/1e3)
hold(h(3),'on');
plot(h(3),d1a(:,1)/1e6,(d1a(:,2)-130e3)*sqrt(2)/1e3)
plot(h(3),d1he(:,1)/1e6,(d1he(:,2)-130e3)*2/1e3)
hold(h(3),'off');
%xlabel(h(3),'n (10^3 km)','fontsize',14)
set(h(3),'XTickLabel',[])
ylabel(h(3),'v_n (km s^{-1})','fontsize',14)
axis(h(3),[-4 2 -1000 000])
irf_legend(h(3),'(d)',[0.98 0.96],'fontsize',14)
set(h(3),'XTickLabel',[])
title(h(3),'MMS-FPI frame')

plot(h(4),d1p(:,1)/1e6,d1p(:,6)/1e3)
hold(h(4),'on');
plot(h(4),d1a(:,1)/1e6,d1a(:,6)/1e3)
plot(h(4),d1he(:,1)/1e6,d1he(:,6)/1e3)
hold(h(4),'off');
xlabel(h(4),'n (10^3 km)','fontsize',14)
ylabel(h(4),'v_{t2} (km s^{-1})','fontsize',14)
axis(h(4),[-4 2 -200 200])
irf_legend(h(4),'(c)',[0.98 0.96],'fontsize',14)

plot(h(5),d1p(:,1)/1e6,(d1p(:,6)+240e3)/1e3)
hold(h(5),'on');
plot(h(5),d1a(:,1)/1e6,(d1a(:,6)+240e3)*sqrt(2)/1e3)
plot(h(5),d1he(:,1)/1e6,(d1he(:,6)+240e3)*2/1e3)
hold(h(5),'off');
xlabel(h(5),'n (10^3 km)','fontsize',14)
ylabel(h(5),'v_{t2} (km s^{-1})','fontsize',14)
axis(h(5),[-4 2 100 800])
irf_legend(h(5),'(e)',[0.98 0.96],'fontsize',14)


%h(6)=irf_panel('Bntt');
irf_plot(h(6),Bntt.x,'k');
hold(h(6),'on')
irf_plot(h(6),Bntt.y,'b');
irf_plot(h(6),Bntt.z,'r');
hold(h(6),'off')
ylabel(h(6),{'B (nT)'},'Interpreter','tex');
irf_legend(h(6),{'B_{n}'},[0.85 0.35],'fontsize',14,'color','k')
irf_legend(h(6),{'B_{t1}'},[0.90 0.35],'fontsize',14,'color','b')
irf_legend(h(6),{'B_{t2}'},[0.95 0.35],'fontsize',14,'color','r')
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',14)
irf_zoom(h(6),'y',[-15 65])
title(h(6),['MMS' num2str(ic)])

specrecn = f1Dn.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n} (s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
%h(7)=irf_panel('f1Dn');
%hold(h(7),'on')
%h(7).HandleVisibility = 'off';
irf_spectrogram(h(7),specrecn);
%hold(h(7),'on')
%h(7).HandleVisibility = 'on';
%irf_plot(h(7),vnproton,'linewidth',2,'color',[0 0.4470 0.7410])
%irf_plot(h(7),vnalpha,'linewidth',2,'color',[0.8500 0.3250 0.0980])
%irf_plot(h(7),vnhelium,'linewidth',2,'color',[0.9290 0.6940 0.1250])
%hold(h(7),'off')
%legend(h(7),{'H^{+}','He^{2+}','He^{+}'},'Location','northeast','AutoUpdate','off')
ylabel(h(7),'v_n (km s^{-1})','interpreter','tex')
colormap(h(7),cmap)
irf_zoom(h(7),'y',vlimn)
caxis(h(7),[-3 2])
irf_legend(h(7),'(g)',[0.99 0.94],'color','k','fontsize',14)

specrect2 = f1Dt2.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2} (s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
%h(8)=irf_panel('f1Dt2');
irf_spectrogram(h(8),specrect2);
%hold(h(8),'on')
%irf_plot(h(8),vt2proton,'linewidth',2,'color',[0 0.4470 0.7410])
%irf_plot(h(8),vt2alpha,'linewidth',2,'color',[0.8500 0.3250 0.0980])
%irf_plot(h(8),vt2helium,'linewidth',2,'color',[0.9290 0.6940 0.1250])
%hold(h(8),'off')
ylabel(h(8),'v_{t2} (km s^{-1})','interpreter','tex')
colormap(h(8),cmap)
irf_zoom(h(8),'y',vlimt2)
caxis(h(8),[-3 2])
irf_legend(h(8),'(h)',[0.99 0.94],'color','k','fontsize',14)

set(h(7:8),'xgrid','off','ygrid','off')
set(h(7:8),'Color',0.7*[1 1 1]);

irf_plot_axis_align(h(6:8));
irf_zoom(h(6:8),'x',Tint);
set(h(1:8),'fontsize',14);

set(gcf,'color','w')

%%
function dydt = tempp(t,q)

l = 10e3;
vsw = -320e3;

n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;

mu0 = 1.2566e-06;
eomi = 9.5788e+07;
e = 1.6022e-19;

By = -B0*tanh(q(1)/l)+B1;
ni = -n0*tanh(q(1)/l)+n1;
Jz = -B0*sech(q(1)/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(q(1)/l).^2./(e*ni*l);
Vx = vsw*nsw./ni;
Ey = 0;
Ez = -Vx*By;

dydt = [q(2) ; eomi*Ex - eomi*q(6)*By; ...
        q(4) ; eomi*Ey; ...
        q(6) ; eomi*Ez + eomi*q(2)*By];
end


function dydt = tempa(t,q)

l = 10e3;
vsw = -320e3;

n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;

mu0 = 1.2566e-06;
eomi = 9.5788e+07/2;
e = 1.6022e-19;

By = -B0*tanh(q(1)/l)+B1;
ni = -n0*tanh(q(1)/l)+n1;
Jz = -B0*sech(q(1)/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(q(1)/l).^2./(e*ni*l);
Ey = 0;
Vx = vsw*nsw./ni;
Ez = -Vx.*By;

dydt = [q(2) ; eomi*Ex - eomi*q(6)*By; ...
        q(4) ; eomi*Ey; ...
        q(6) ; eomi*Ez + eomi*q(2)*By];
end

function dydt = temphe(t,q)

l = 10e3;
vsw = -320e3;

n0 = 4*1e6; 
n1 = 13*1e6; 
B0 = 13.5*1e-9; 
B1 = 45.5*1e-9; 
nsw = n1-n0;

P0 = 0.0425e-9;
P1 = 0.0575e-9;

mu0 = 1.2566e-06;
eomi = 9.5788e+07/4;
e = 1.6022e-19;

By = -B0*tanh(q(1)/l)+B1;
ni = -n0*tanh(q(1)/l)+n1;
Jz = -B0*sech(q(1)/l).^2/(mu0*l);
Ex = -Jz.*By./(e*ni)+P0*sech(q(1)/l).^2./(e*ni*l);
Ey = 0;
Vx = vsw*nsw./ni;
Ez = -Vx.*By;

dydt = [q(2) ; eomi*Ex - eomi*q(6)*By; ...
        q(4) ; eomi*Ey; ...
        q(6) ; eomi*Ez + eomi*q(2)*By];
end