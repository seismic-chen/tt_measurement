clear; clc; close all;
% Signal processing package
addpath ./processRFmatlab-master/
processRFmatlab_startup
% Taup package
javaaddpath('./FMI/lib/FMI.jar');
addpath('./FMI/matTaup');
% Plotting package
addpath(fullfile(userpath,'m_map'))
%% load data
load fetchdata.mat;
load rawdata.mat;
%% remove the anomalous traces
xshift=1;
rmsd0=rms(d0);
figure
for n=1:size(d0,2)
    if rmsd0(n)<2000
        plot(d0(:,n)/max(d0(:,n))+(n-1)*xshift,1:size(d0,1),'k');hold on;
    else
        plot(d0(:,n)/max(d0(:,n))+(n-1)*xshift,1:size(d0,1),'r');hold on;
    end
end
index_del=find(rmsd0>2000);
d0(:,index_del)=[];
stla0(index_del)=[];
stlo0(index_del)=[];
startTime(index_del)=[];
stationname(index_del)=[];
% select stations
lonlim = [-116.75 -101.25];
latlim = [32.25 49.75];
tP=[];
tS=[];
keep=stlo0<-100 & stlo0>-117 & stla0>33 & stla0<49;
stla0=stla0(keep);
stlo0=stlo0(keep);
stationname=stationname(keep);
startTime=startTime(keep);
for n=1:length(stla0)
    stla=stla0(n);
    stlo=stlo0(n);
    tt=taupTime('ak135',event_info.PreferredDepth,'P,Pdiff','sta',[stla, stlo],'evt',[evla, evlo]);
    try
        tP(n)=tt(1).time;
    catch me
        tP(n)=NaN;
    end
end
%% filtering traces
dt=1; % note that dt is 1 sec in this example
traces=d0(:,keep);
[nt,nx]=size(traces);
flow = 0.03;
fhigh = 0.125;
traces_filt=zeros(nt,nx);
for n=1:nx
    seis=traces(:,n);
    seis = removeSeisDC(seis);
    seis = removeSeisTrend(seis);
    seis_taper = taperSeis(seis,0.2);
    seis_filter = bandpassSeis(seis_taper,dt,flow,fhigh,3);
    traces_filt(:,n)=seis_filter/max(seis_filter)*0.1;
end
%% measuring the P wave travel time residual using MCCC
tshift=600;  % requesting data 600 sec before origin time
trace_cut=[];
t=(0:nt-1)*dt;
t=repmat(t(:),1,nx);
% xcorr measurement
xwin=[-50 50];
maxlag=5;   % maximum time lag in cross-correlation
t0=tP+tshift;  %phase arrival time
N=size(traces_filt,2);
% upsample
ups = 4;
dns = 1;
traces_filt1 = resample(traces_filt,ups,dns);
dt1=dt/4;
[nt,nx]=size(traces_filt);
t1=[];
for i=1:nx
    t1(:,i)=t(1,n):dt1:t(end,n)+1-dt1;
end
[t_est0, tdelay0, sigma0]=mccc(traces_filt1,t1,t0,xwin,maxlag);
%% measuring the P wave travel time residual using adaptive stacking
ncorr=5;  % number of cross-correlation measurments
is_plot=0;
t_est1=astack(traces_filt1,t1,t0,xwin,maxlag,ncorr,is_plot);
%% compare two travel time picking method
t_est0=t_est0(:);
t0=t0(:);
ttres0=(t_est0-mean(t_est0))-(t0-mean(t0));
ttres1=(t_est1-mean(t_est1))-(t0-mean(t0));
figure;
plot(ttres0,'k-','linewidth',1); hold on;
plot(ttres1,'-')
%% m_map plot
figure;
subplot(131)
set(gcf,'Position',[100 100 1600 800])
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
h=m_scatter(stlo0,stla0,200,t0-mean(t0),'filled');
m_gshhs('l','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k') 
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',10);
m_ruler([.05 .3],.05, 2,'fontsize',8)
title('Predicted travel time','fontsize',18)
hh=colorbar('h');
xlabel(hh,'Time (sec)');
colormap(flipud(jet));
caxis([-40 40])
subplot(132)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
m_scatter(stlo0,stla0,200,-ttres0,'filled')
m_gshhs('l','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k') 
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',10);
m_ruler([.05 .3],.05, 2,'fontsize',8)
title('Relative travel time residual (MCCC)','fontsize',18)
hh=colorbar('h');
xlabel(hh,'Time (sec)');
caxis([-2 2])
subplot(133)
m_proj('lambert','long', lonlim, 'lat', latlim); hold on;
m_scatter(stlo0,stla0,200,-ttres1,'filled')
m_gshhs('l','line','color','k','linewidth',1)
m_gshhs('lb2','line','color','k') 
m_grid('linewidth',2,'tickdir','out',...
    'xaxisloc','bottom','yaxisloc','right','fontsize',10);
m_ruler([.05 .3],.05, 2,'fontsize',8)
caxis([-2 2])
title('Relative travel time residual (adaptive stacking)','fontsize',18)
hh=colorbar('h');
xlabel(hh,'Time (sec)');