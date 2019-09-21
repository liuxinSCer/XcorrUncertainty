function [xcorrarray,xdigit,hFIR2,subindspec]=getBootXcorrFiltered(normbootarray,freqrangedata,vfreq)
DF=freqrangedata(2)-freqrangedata(1);
if vfreq==0
    filterFIR=1;% 1==FIR filter; 1: FIR band-pass filter w/ flat response
else
    filterFIR=0;% 0==Narrowband Gaussian filter
    ifreqdistri=round((vfreq-freqrangedata(1))/DF+1);

end
% Gaussian filter freq width (1-sigma)
sigfreq=0.03; % Hz

 freqwin=[freqrangedata(1) freqrangedata(end)];
 
% freqmeasure=0.07:DF*2:0.59;
% freqmeasure=0.07:DF:0.59;
% inddistrifreq=round((freqmeasure-freqrangedata(1) )/DF+1);

% DF=freqrangedata(2)-freqrangedata(1);
nfreqstore=length(freqrangedata);
subindspec=round(freqwin(1)/DF)+1:round(freqwin(2)/DF)+1;


Fsamp=4;Fsamp=2; Fny=Fsamp/2; bdwd=0.02;%bdwd=0.04;bdwd=0.1;
dt=1/Fsamp;

Fbandwidth=0.02; % in Hz; minimum bandwidth between pass to cut-off frequency.
BW=Fbandwidth/Fsamp; % normalized bandwidth by sampling frequency Fsamp.
Mfilter=4/BW; % theoretical number of points for filter kernal to satisfy the bandwidth requirement.
% filterpts=50;
filterpts=Mfilter; %filterpts=100; % number of order for the FIR filter

% freqindwinspecdata=round([freqv(1) freqv(end)]/DFxcorr+1);
% freqindspecdata=freqindwinspecdata(1):freqindwinspecdata(2);

% NFsamp=1;
npts=round(1*Fsamp/DF+1);
halfnpts=(npts+1)/2;
midp=halfnpts; % mid point of correlation time window; just an alias names


xdigit=(-halfnpts+1:halfnpts-1)*dt;
freqwinplot=freqwin;
% freqwinplot(1)=0.03;
fwin=(freqwinplot+[bdwd+0.0 -bdwd])/Fny;



%% FILTER PARAMETERS
% ifreqm=4; % 0.20 Hz
b = fir1(filterpts,fwin);
[hFIR,w]=freqz(b,1,halfnpts);
wfreq=w*Fny/pi;
if filterFIR

    %hFIR2=hFIR.*conj(hFIR);
    hFIR2=abs(hFIR);

else % Gaussian filter
    
%     centerfreq=distrifreq(ifreqm);
    centerfreq=freqrangedata(ifreqdistri)
    
    %         halffreqwidth=min([centerfreq/5  0.10]); % half of the freqwidth
    %         halffreqwidth=min([ 0.05]); % half of the freqwidth for filter
    
    % NOW TRY NARROW-BAND GAUSSIAN FILTER
    gaussnarrow=exp(-(wfreq-centerfreq).^2/(2*sigfreq*sigfreq));
    
    hFIR2=gaussnarrow;
    
    
end



nboots=size(normbootarray,1);
% nboots=min(2000,nboots);

wfxcorrnew=zeros(npts,1);
xcorrarray=zeros(npts,nboots);
for iboot=1:nboots
    xspectra=zeros(npts,1);
    xspectra(subindspec)=normbootarray(iboot,1:end);
    xspectra(1:midp)=xspectra(1:midp).*hFIR2;
    %xspectraratio(1:midp)=xspectraratio(1:midp).*hFIR2;
    % freqindwin=round(freqwinplot/df)+1;
    % xspectra(1:freqindwin(2))=detrend(xspectra(1:freqindwin(2)),'constant');
    
    % THE LINE BELOW ENSURES HERMITIAN SYMMETRY; DISABLEING IT WILL MAKE
    % ANALYTICAL SIGNAL:)
%     xspectra(midp+1:end)=conj(xspectra(midp:-1:2));
    
    wfxcorrnew(:)=ifft(xspectra);
    wfxcorrnew=[wfxcorrnew(midp+1:end); wfxcorrnew(1:midp) ];
    xcorrarray(:,iboot)=wfxcorrnew;
    
end

xcorrarray=xcorrarray';