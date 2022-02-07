function [t_est,tdelay, sigma] = mccc(din,tin,t0,xwin,maxlag)
% Feb. 7, 2022, Yunfeng Chen, Global Seismology Group, Zhejiang University
% Travel-time measurement using multi-channel cross-correlation method
% Reference: VanDecar and Cross, Determination of telesiesmic relative
% phase arrival times using multi-channel cross-correlation and leaasst
% squares, 1990, BSSA.
% Input:
% din: seismic traces stored in column
% tin: time axis, the same size as traces
% t0: P wave arrival time
% xwin: cross correlation window length, relative to P arrival time
% maxlag: maximum lag allowed in cross correlation measurement (in sec)
% Output:
% tdelay: timeshif that leads to the optimal cross correlation coefficients
% the strcuture of tdelay is a upper triangle matrix, with the diagnal
% elements represent the auto correlation results
% t_est: least square optimized arrival time
% sigma: uncertainty in travel time measurements

% cut the waveform arround P
[nt,nx] = size(din);
yshift=0;
tp=zeros(nx,nx);
for i = 1:nx
    start = t0(i)+xwin(1);
    finish = t0(i)+xwin(2);
    d=din(:,i);
    t=tin(:,i);
    [d_cut, t_cut] = chopSeis(d, t, start, finish);
    % debug
%     plot(t,d+yshift,'b'); hold on; plot(t_cut,d_cut+yshift,'r','linewidth',1.2)
%     plot(d_cut+yshift,'r'); hold on;
%     yshift=yshift+0.1;
    din_cut(:,i) = d_cut/max(d_cut);
    % time difference in t0
    for j = 1:nx
       tp(i,j)=t0(j)-t0(i); 
    end
end

dt = t(2)-t(1);
nlag = round(maxlag/dt);
[C,lags] = xcorr(din_cut,nlag);
[~,index] = max(C);
I = reshape(index,nx,nx)';
% find the optimal delay time
tdelay = (lags(I(:)))*dt;
tdelay = triu(reshape(tdelay,nx,nx)-tp);
% least squares optimization
A = zeros(nx*(nx-1)/2+1,nx);
countj = 1;
counti = 1;
for i = 1:size(A,1)
    for j = 1:size(A,2)
        if j == counti+countj
           A(i,j) = -1; 
        end
        if j == counti
            A(i,j) = 1;
        end
    end
    if countj == size(A,2)-counti
        counti = counti + 1;
        countj = 0;
    end
    countj = countj+1;
end
A(end,:) = ones(1,nx);
% solve the least square problem to obtain the estimated travel time
k=1;
d=zeros(nx*(nx-1)/2+1,1);
for i = 1:size(tdelay,1)
    for j = 1:size(tdelay,2)
        if j > i
            d(k) = tdelay(i,j);
            k = k+1;
        end
    end
end
t_est = (A'*A)\A'*d;
res=d-A*t_est;
% compute the standard deviation of the time measurements
res1=zeros(nx,nx);
k=0;
for i=1:nx
    for j=1:nx
        if j>i
            k=k+1;
            res1(i,j)=res(k);
        end
    end
end
res1=res1+res1';
sigma=sqrt((sum(res1,2).^2)/(nx-2));
% design weighting function based on residual
% W=[abs(res); 1];
% for i=1:size(A,1)
%    AW(i,:)=A(i,:)*W(i);
%    dW(i,1)=d(i)*W(i);
% end
% t_est_w=(AW'*AW)\AW'*dW;

% % plot the traces
% tnew = (0:nt-1)*dt;
% tnew=tnew(:);
% for i = 1:nx
%     % apply timeshift to the traces
%     tdiff = t_est(1) - t_est(i);
%     data = din(:,n);
%     datashift = fftShift(data,tnew,tdiff);
%     din_shift(:,i) = datashift;
% end
