%% 
Tint = irf.tint('2023-04-23T23:30:00.00Z/2023-04-24T05:00:00.00Z');
ic = 1;
iPDist = mms.get_data('PDi_fpi_fast_l2',Tint,ic);
ePDist = mms.get_data('PDe_fpi_fast_l2',Tint,ic);
Bxyz=mms.get_data('B_gse_srvy_l2',Tint,ic);

Tints = irf.tint('2023-04-24T03:49:34.00Z/2023-04-24T03:51:02.00Z');
Bxyzb=mms.get_data('B_gse_brst_l2',Tints,ic);
iPDistb = mms.get_data('PDi_fpi_brst_l2',Tints,ic);
ePDistb = mms.get_data('PDe_fpi_brst_l2',Tints,ic);

%%

ePDistomni = ePDist.omni.deflux;
iPDistomni = iPDist.omni.deflux;

ePDistomnib = ePDistb.omni.deflux;
iPDistomnib = iPDistb.omni.deflux;

%% 

%np = mms.db_get_ts('mms1_hpca_brst_l2_moments','mms1_hpca_hplus_number_density',Tint);
%Vp = mms.db_get_ts('mms1_hpca_brst_l2_moments','mms1_hpca_hplus_ion_bulk_velocity',Tint);
%nhep = mms.db_get_ts('mms1_hpca_brst_l2_moments','mms1_hpca_heplus_number_density',Tint);
%Vhep = mms.db_get_ts('mms1_hpca_brst_l2_moments','mms1_hpca_heplus_ion_bulk_velocity',Tint);
%nhepp = mms.db_get_ts('mms1_hpca_brst_l2_moments','mms1_hpca_heplusplus_number_density',Tint);
%Vhepp = mms.db_get_ts('mms1_hpca_brst_l2_moments','mms1_hpca_heplusplus_ion_bulk_velocity',Tint);


np = mms.db_get_ts('mms1_hpca_srvy_l2_moments','mms1_hpca_hplus_number_density',Tint);
Vp = mms.db_get_ts('mms1_hpca_srvy_l2_moments','mms1_hpca_hplus_ion_bulk_velocity',Tint);
nhep = mms.db_get_ts('mms1_hpca_srvy_l2_moments','mms1_hpca_heplus_number_density',Tint);
Vhep = mms.db_get_ts('mms1_hpca_srvy_l2_moments','mms1_hpca_heplus_ion_bulk_velocity',Tint);
nhepp = mms.db_get_ts('mms1_hpca_srvy_l2_moments','mms1_hpca_heplusplus_number_density',Tint);
Vhepp = mms.db_get_ts('mms1_hpca_srvy_l2_moments','mms1_hpca_heplusplus_ion_bulk_velocity',Tint);

%tmpDataObj = dataobj('/Volumes/mms/mms1/hpca/srvy/l2/moments/2023/04/mms1_hpca_srvy_l2_moments_20230424000000_v4.3.5.cdf');
%np = get_ts(tmpDataObj,'mms1_hpca_hplus_number_density');
%Vp = get_ts(tmpDataObj,'mms1_hpca_hplus_ion_bulk_velocity');
%nhep = get_ts(tmpDataObj,'mms1_hpca_heplus_number_density');
%Vhep = get_ts(tmpDataObj,'mms1_hpca_heplus_ion_bulk_velocity');
%nhepp = get_ts(tmpDataObj,'mms1_hpca_heplusplus_number_density');
%Vhepp = get_ts(tmpDataObj,'mms1_hpca_heplusplus_ion_bulk_velocity');
%Thp = get_ts(tmpDataObj,'mms1_hpca_hplus_scalar_temperature');
%Thep = get_ts(tmpDataObj,'mms1_hpca_heplus_scalar_temperature');
%Thepp = get_ts(tmpDataObj,'mms1_hpca_heplusplus_scalar_temperature');

ni = mms.get_data('Ni_fpi_brst_l2',Tints,ic);
Vi = mms.get_data('Vi_gse_fpi_brst_l2',Tints,ic);

%% 
h=irf_plot(8,'newfigure');
%h=irf_figure(540+ic,8);
xSize=750; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.90;
ywidth = 0.11;
set(h(1),'position',[0.085 0.975-ywidth xwidth+0.02 ywidth]);
set(h(2),'position',[0.085 0.975-2*ywidth xwidth+0.02 ywidth]);
set(h(3),'position',[0.085 0.975-3*ywidth xwidth+0.02 ywidth]);
set(h(4),'position',[0.085 0.935-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.085 0.935-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.085 0.935-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.085 0.935-7*ywidth xwidth ywidth]);
set(h(8),'position',[0.085 0.935-8*ywidth xwidth ywidth]);

Bxyzmag = irf.ts_scalar(Bxyz.time,[Bxyz.data Bxyz.abs.data]);

h(1)=irf_panel('idist');
irf_spectrogram(h(1),iPDistomni.specrec,'log');
irf_legend(h(1),'(a)',[0.99 0.98],'color','k','fontsize',14);
set(h(1),'yscale','log');
set(h(1),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(1),'E_{i} (eV)','fontsize',14,'Interpreter','tex');
title(h(1),'MMS1')

h(2)=irf_panel('edist');
irf_spectrogram(h(2),ePDistomni.specrec,'log');
irf_legend(h(2),'(b)',[0.99 0.98],'color','k','fontsize',14);
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(2),'E_{e} (eV)','fontsize',14,'Interpreter','tex');

h(3)=irf_panel('Bxyz');
irf_plot(h(3),Bxyzmag);
ylabel(h(3),{'B (nT)'},'Interpreter','tex');
irf_legend(h(3),{'B_{x}','B_{y}','B_{z}','|B|'},[0.70 0.7],'fontsize',14)
irf_legend(h(3),'(c)',[0.99 0.7],'color','k','fontsize',14)
irf_zoom(h(3),'y',[-70 80])

irf_plot_axis_align(h(1:3));
irf_zoom(h(1:3),'x',Tint);

Bxyzbmag = irf.ts_scalar(Bxyzb.time,[Bxyzb.data Bxyzb.abs.data]);

h(4)=irf_panel('Bxyzb');
irf_plot(h(4),Bxyzbmag);
ylabel(h(4),{'B (nT)'},'Interpreter','tex');
irf_legend(h(4),{'B_{x}','B_{y}','B_{z}','|B|'},[0.9 0.12],'fontsize',14)
irf_legend(h(4),'(d)',[0.99 0.94],'color','k','fontsize',14)
irf_zoom(h(4),'y',[-60 65])

h(5)=irf_panel('idistb');
irf_spectrogram(h(5),iPDistomnib.specrec,'log');
irf_legend(h(5),'(e)',[0.99 0.6],'color','w','fontsize',14);
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(5),'E_{i} (eV)','fontsize',14,'Interpreter','tex');

h(6)=irf_panel('edistb');
irf_spectrogram(h(6),ePDistomnib.specrec,'log');
irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',14);
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3 1e4]);
ylabel(h(6),'E_{e} (eV)','fontsize',14,'Interpreter','tex');

h(7)=irf_panel('ni');
irf_plot(h(7),ni);
hold(h(7),'on')
irf_plot(h(7),np);
irf_plot(h(7),nhepp);
irf_plot(h(7),nhep);
hold(h(7),'off')
set(h(7),'yscale','log')
ylabel(h(7),{'n_i (cm^{-3})'},'Interpreter','tex');
irf_legend(h(7),{'n_{i} ','n_{H+}','n_{He2+}','n_{He+}'},[0.90 0.77],'fontsize',14)
irf_legend(h(7),'(g)',[0.99 0.90],'color','k','fontsize',14)
irf_zoom(h(7),'y',[1e-2 40])
set(h(7),'ytick',[1e-2 1e-1 1e0 1e1]);

h(8)=irf_panel('Vi');
irf_plot(h(8),Vi.x);
hold(h(8),'on')
irf_plot(h(8),Vp.x);
irf_plot(h(8),Vhepp.x);
irf_plot(h(8),Vhep.x);
hold(h(8),'off')
ylabel(h(8),{'V_{i,x} (km s^{-1})'},'Interpreter','tex');
irf_legend(h(8),{'V_{i} ','V_{H+}','V_{He2+}','V_{He+}'},[0.90 0.94],'fontsize',14)
irf_legend(h(8),'(h)',[0.99 0.94],'color','k','fontsize',14)
irf_zoom(h(8),'y',[-700 -250])

irf_plot_axis_align(h(4:8));
irf_zoom(h(4:8),'x',Tints);

colormap(h(1),'jet');
colormap(h(2),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');

irf_plot_zoomin_lines_between_panels(h(3),h(4));
c_eval('irf_pl_mark(h(?),irf_time(Tints,''epochtt>epoch'')'',[230 200 0]/255)',1:3);

set(h(1:8),'fontsize',14);
