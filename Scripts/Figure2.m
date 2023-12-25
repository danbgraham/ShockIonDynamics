Tint = irf.tint('2023-04-24T03:49:35.00Z/2023-04-24T03:50:30.00Z');

ic = 1;
Bxyz=mms.get_data('B_dmpa_brst_l2',Tint,ic);
ne = mms.get_data('Ne_fpi_brst_l2',Tint,ic);
ni = mms.get_data('Ni_fpi_brst_l2',Tint,ic);

Te = mms.get_data('Te_dbcs_fpi_brst_l2',Tint,ic);
Pe = mms.get_data('Pe_dbcs_fpi_brst_l2',Tint,ic);

Vi = mms.get_data('Vi_dbcs_fpi_brst_l2',Tint,ic);
Ti = mms.get_data('Ti_dbcs_fpi_brst_l2',Tint,ic);
Tis = Ti.trace/3;

Tipp = mms.rotate_tensor(Ti,'fac',Bxyz,'qq');
Tipar = irf.ts_scalar(Ti.time,Tipp.xx.data);
Tiperp1 = irf.ts_scalar(Ti.time,Tipp.yy.data);
Tiperp2 = irf.ts_scalar(Ti.time,Tipp.zz.data);

Tippp = irf.ts_scalar(Ti.time,[Tiperp1.data Tiperp2.data Tipar.data]);

iPDist = mms.get_data('PDi_fpi_brst_l2',Tint,ic);
c_eval('SCpot=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',Tint);',ic);
%% 

n = [0.926 0.102 -0.364];
t1 = [-0.189 -0.708 -0.680];
t2 = [-0.327 0.699 -0.636];

ndsl = [0.9441 0.0997 -0.3142];
t1dsl = [-0.1524 -0.7125 -0.6849];
t2dsl = [-0.2920 0.6949 -0.6571];
ndsl = ndsl/norm(ndsl);
t1dsl = t1dsl/norm(t1dsl);
t2dsl = t2dsl/norm(t2dsl);

Bntt = irf_newxyz(Bxyz,ndsl,t1dsl,t2dsl);
Vintt = irf_newxyz(Vi,ndsl,t1dsl,t2dsl);

Tints = irf.tint('2023-04-24T03:50:10.00Z/2023-04-24T03:50:13.50Z');
Exyz = mms.db_get_ts('mms1_edp_brst_l2_hmfe','mms1_edp_hmfe_dsl_brst_l2',Tints);
Tints = irf.tint(Exyz.time);

Tints2 = Tints+[-0.1 0.1];


%% 

Units = irf_units;
Me=Units.me;
Mi=Units.mp;
e=Units.e;
epso=Units.eps0;
Wpe = sqrt(ne.data*1e6*e^2/Me/epso);
FpeTS = irf.ts_scalar(ne.time,Wpe/2/pi);
FpiTS = FpeTS*sqrt(Me/Mi);

nf = 400;
Ewavelet = irf_wavelet(Exyz,'nf',nf,'f',[50 32e3],'wavelet_width',5.36*20);

nc = 10;
idx = [nc/2:nc:length(Ewavelet.t)-nc/2];
Ewavelettimes = Ewavelet.t(idx);
Ewaveletx = zeros(length(idx),nf);
Ewavelety = zeros(length(idx),nf);
Ewaveletz = zeros(length(idx),nf);
for ii = [1:length(idx)];
        Ewaveletx(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,1}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewavelety(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,2}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
        Ewaveletz(ii,:) = squeeze(irf.nanmean(Ewavelet.p{1,3}([idx(ii)-nc/2+1:idx(ii)+nc/2-1],:),1));
end
specE=struct('t',Ewavelettimes);
specE.f=Ewavelet.f/1000;
specE.p=Ewaveletx+Ewavelety+Ewaveletz;
specE.f_label='';
specE.p_label={'log_{10} E^2','mV^2 m^{-2} Hz^{-1}'};

maxpower = max(max(Ewaveletx+Ewavelety+Ewaveletz));

%% Electrondists
ePDist = mms.get_data('PDe_fpi_brst_l2',Tints2,ic);
ePDistrmb = mms.remove_edist_background(ePDist,'tint',Tints2,'nSecondary', 0.5);

vg = -30e3:1e3:30e3;

nMC = 2e2;

Bvec = Bxyz.resample(ePDist);
Bvec = Bvec/Bvec.abs;

SCpot = SCpot.resample(ePDist);

ePDistrmb = ePDistrmb.convertto('s^3/cm^6');
ef1D = ePDistrmb.reduce('1D',Bvec,'vg',vg,'scpot',SCpot,'nMC',nMC);

ef1Ds = ef1D.specrec('1D_velocity');
ef1Ds.f = ef1Ds.f/1e3;
ef1Ds.p_label={'log_{10}f_{e,||}','(s m^{-4})'};

%% Get 1D reduced distributions
iPDist.data(:,1:10,:,:) = 0;

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

vvecn = f1Dn.depend{1}(1,:);
dv = median(diff(vvecn))*1e3;
idxn = vvecn < -300 & vvecn > -550;
nibulk = squeeze(sum(f1Dn.data(:,idxn),2))*dv;
vnbulk = squeeze(sum(f1Dn.data(:,idxn).*f1Dn.depend{1}(:,idxn)*dv,2))./nibulk;

%% Figure
h=irf_plot(9,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.85;
ywidth = 0.095;
ywidth2 = 0.1;
set(h(1),'position',[0.12 0.975-ywidth xwidth ywidth]);
set(h(2),'position',[0.12 0.975-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.12 0.975-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.12 0.975-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.12 0.975-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.12 0.06+ywidth2*3 xwidth ywidth2]);
set(h(7),'position',[0.12 0.06+ywidth2*2 xwidth ywidth2]);
set(h(8),'position',[0.12 0.06+ywidth2 xwidth ywidth2]);
set(h(9),'position',[0.12 0.06 xwidth ywidth2]);

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

h(1)=irf_panel('Bntt');
irf_plot(h(1),Bntt);
ylabel(h(1),{'B','(nT)'},'Interpreter','tex');
irf_zoom(h(1),'y',[-10 65]);
irf_legend(h(1),{'B_{n}','B_{t1}','B_{t2}'},[0.90 0.94],'fontsize',14)
irf_legend(h(1),'(a)',[0.99 0.94],'color','k','fontsize',14)
title(h(1),'MMS1')

h(2)=irf_panel('ne');
irf_plot(h(2),ne);
ylabel(h(2),{'n_e','(cm^{-3})'},'Interpreter','tex');
irf_legend(h(2),'(b)',[0.99 0.94],'color','k','fontsize',14)

specrecn = f1Dn.specrec('1D_velocity');
specrecn.p_label={'log_{10}f_{i,n}','(s m^{-4})'};
specrecn.p(specrecn.p < 1e-3) = NaN;
h(3)=irf_panel('f1Dn');
irf_spectrogram(h(3),specrecn);
hold(h(3),'on')
irf_plot(h(3),Vintt.x)
hold(h(3),'off')
ylabel(h(3),{'v_n','(km s^{-1})'},'interpreter','tex')
colormap(h(3),cmap)
irf_zoom(h(3),'y',vlimn)
caxis(h(3),[-3 2])
irf_legend(h(3),'(c)',[0.99 0.94],'color','k','fontsize',14)

specrect1 = f1Dt1.specrec('1D_velocity');
specrect1.p_label={'log_{10}f_{i,t1}','(s m^{-4})'};
specrect1.p(specrect1.p < 1e-3) = NaN;
h(4)=irf_panel('f1Dt1');
irf_spectrogram(h(4),specrect1);
hold(h(4),'on')
irf_plot(h(4),Vintt.y)
hold(h(4),'off')
ylabel(h(4),{'v_{t1}','(km s^{-1})'},'interpreter','tex')
colormap(h(4),cmap)
irf_zoom(h(4),'y',vlimt1)
caxis(h(4),[-3 2])
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',14)

specrect2 = f1Dt2.specrec('1D_velocity');
specrect2.p_label={'log_{10}f_{i,t2}','(s m^{-4})'};
specrect2.p(specrect2.p < 1e-3) = NaN;
h(5)=irf_panel('f1Dt2');
irf_spectrogram(h(5),specrect2);
hold(h(5),'on')
irf_plot(h(5),Vintt.z)
hold(h(5),'off')
ylabel(h(5),{'v_{t2}','(km s^{-1})'},'interpreter','tex')
colormap(h(5),cmap)
irf_zoom(h(5),'y',vlimt2)
caxis(h(5),[-3 2])
irf_legend(h(5),'(e)',[0.99 0.94],'color','k','fontsize',14)

irf_plot_axis_align(h(1:5));
irf_zoom(h(1:5),'x',Tint);
irf_timeaxis(h(1:5),'nodate')

h(6)=irf_panel('Bntts');
irf_plot(h(6),Bntt);
ylabel(h(6),{'B','(nT)'},'Interpreter','tex');
irf_zoom(h(6),'y',[-10 65]);
irf_legend(h(6),{'B_{n}','B_{t1}','B_{t2}'},[0.90 0.94],'fontsize',14)
irf_legend(h(6),'(f)',[0.99 0.94],'color','k','fontsize',14)

h(7)=irf_panel('f1Dns');
irf_spectrogram(h(7),specrecn);
ylabel(h(7),{'v_n','(km s^{-1})'},'interpreter','tex')
colormap(h(7),cmap)
irf_zoom(h(7),'y',vlimn)
caxis(h(7),[-3 2])
irf_legend(h(7),'(g)',[0.99 0.94],'color','k','fontsize',14)

ef1Ds.p(ef1Ds.p < 1e-5) = NaN;
h(8) = irf_panel('electron1D');
irf_spectrogram(h(8),ef1Ds,'log');
ylabel(h(8),{'v_{e||}','(10^3 km s^{-1})'},'interpreter','tex')
irf_zoom(h(8),'y',[-20,20])
irf_legend(h(8),'(h)',[0.99 0.94],'color','k','fontsize',14)
colormap(h(8),cmap)
caxis(h(8),[-5 1])

h(9)=irf_panel('Espec');
irf_spectrogram(h(9),specE,'log');
hold(h(9),'on');
irf_plot(h(9),irf.ts_scalar(FpeTS.time,FpeTS.data/1000),'color','w','LineWidth',2.0)
hold(h(9),'off');
irf_legend(h(9),'(i)',[0.99 0.94],'color','w','fontsize',14)
irf_legend(h(9),'f_{pe}',[0.56 0.9],'color','w','fontsize',14)
caxis(h(9),[log10(maxpower)-10 log10(maxpower)]);
ylabel(h(9),{'f','(kHz)'},'fontsize',14,'Interpreter','tex');
colormap(h(9),'jet')

irf_plot_axis_align(h(6:9));
irf_zoom(h(6:9),'x',Tints);

set(h([3:5 7:8]),'xgrid','off','ygrid','off')
set(h([3:5 7:8]),'Color',0.7*[1 1 1]);

irf_plot_zoomin_lines_between_panels(h(5),h(6));
c_eval('irf_pl_mark(h(?),irf_time(Tints,''epochtt>epoch'')'',[255 255 100]/255)',1:5);

set(h(1:8),'fontsize',14);
