function tdelay=astack(din,tin,t0,xwin,maxlag,ncorr,is_plot)
% Feb. 3, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% adpative-stacking based cross-correlation travel time measurements.
% Reference: Rawlinson and Kennett (2004), Rapid estimation of relative and 
% absolute delay times across a network by adaptive stacking, GJI.
% Input:
% din: seismic traces stored in column
% tin: time axis, the same size as traces
% t0: predicted P wave arrival time
% xwin: cross correlation window length, relative to P arrival time
% maxlag: maximum time lag allowed in cross correlation measurement (in sec)
% ncorr: number of cross-correltion measurments
% Output:
% tdelay: the measured time difference of each trace relative to the regional
% average
[nt,nx] = size(din);
t0=t0(:);
tdiff=zeros(nx,1);
trace_ref=[];
tlag=[];
tc=t0;
dt = tin(2,1)-tin(1,1);
nlag=round(maxlag/dt);
for icorr = 1:ncorr
    % cut the waveform around P
    din_cut=[];
    tin_cut=[];
    for i = 1:nx
        start =  tc(i)+xwin(1);
        finish = tc(i)+xwin(2);
        d=din(:,i);
        t=tin(:,i);
        [d_cut, t_cut] = chopSeis(d, t, start, finish);
        din_cut(:,i) = d_cut/max(d_cut);
    end
    if icorr==1
        % save a copy of original data
        din0_cut=din_cut;
    end
    % reference trace
    din_ref(:,icorr) = mean(din_cut,2);
    for i=1:nx
        [C,lags] = xcorr(din_cut(:,i),din_ref(:,icorr),nlag);
        [~,index] = max(C);
        % find the optimal delay time
        tdiff(i) = lags(index)*dt;
        % use the grid search method from Rawlinson and Kennett (2003)?
    end
    tlag=[tlag tdiff];
    tc=t0+sum(tlag,2);
    if is_plot
        figure;
        t=(0:nt-1)*dt;
        imagesc(1:size(din_cut,2),t,[din0_cut din_cut])
        colormap(seismic(3))
        caxis([-1,1])
    end
end
tdelay=t0+sum(tlag,2);