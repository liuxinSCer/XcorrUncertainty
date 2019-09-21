function [distrixcorr,errvecboth,freqrange]=createXspecMat(wftime,dt,winlen,freqwin,preprocstr,ndays)
%%%%%%
% COMPUTE ALL CROSS-SPECTRA OBSERVATIONS FOR SLIDING WINDOWS
% wftime: npts*2 matrix containing two noise time series
% winlen: window length in sec for sliding windows
% ndays: No. of days to use:: can be any float number
% by Xin Liu
% liuxin@stanford.edu
% All rights reserved, Dec 2017
%%%%%%
ploteachwf = true;

npts=length(wftime(:,1));
fontsize = 10;
% % % SPECIFY NUMBER OF DAYS BY LIMITING THE NUMBER OF DAYS
rightcorner = ndays *(24*60*60/dt);% number of days in use
%  rightcorner=npts/5;  % 5 DAYS0
%    rightcorner=npts; % 25 DAYS
% % % % % % ndays=rightcorner/(24*60*60/dt);

% if applyprewhiten
%     preprocstr='prewhiten';
%
% elseif applyonebitnorm
%
%     preprocstr='onebit'
% else
%     preprocstr='';
% end

if strcmp(preprocstr,'onebit')
    
    applyonebitnorm=true;
elseif  strcmp(preprocstr,'prewhiten')
    applyprewhiten = true;
elseif  strcmp(preprocstr,'both') % both prewhiten & onebit
    applyprewhiten = true;
    applyonebitnorm=true;
else
    applyprewhiten = false;
    applyonebitnorm = false;
    
end
%% NOW COMPUTE XCORR AND STATISTICS
%%% TODO: define the hopping window size and hopping length
Fsamp=1/dt;
% DF=0.01;%0.02;
DF = 1.0/winlen;
DFxcorr=DF/2; % since the xcorr length is twice the segment length
SegLenFreq=Fsamp*1/DF; % LENGTH OF EACH TIME.D SEGMENT in number of points
LenFFT=2*SegLenFreq+1; % number of points for FFT
gapT=winlen*0.2; % gap length is 20% of window length
hopsize=Fsamp*(gapT);%Fsamp*20;Fsamp*0; % hopping size of windows; 20 sec;

freqv=(0:Fsamp/DFxcorr)*DFxcorr;
nfreqv=length(freqv);
% midfreqv=(nfreqv-1)/2+1;

% fourthM=0;% 4th moment
xcorrseg=0;% xcorr
xcorrnoise=0'; % the pure uncorr-noise term
psd1=0;% PSD of sta1
psd2=0;% PSD of sta2

psdnoise1=0; psdnoise2=0;% only the psd for uncorrelated noise
ratioalpha=0;
% stack numerator and denominator separately

% now compute the probablity distribution
% niter=ceil(  npts/SegLenFreq);
niter=ceil(  rightcorner/SegLenFreq);
% distrifreq=[0.2 0.3];
freqrange=freqwin(1):DFxcorr:freqwin(2); % FOR NORMAL USE!
% freqrange=freqwin(1):DFxcorr*5:freqwin(2); % FOR NORMAL USE!
distrifreq=freqrange;
% distrifreq=0.1:0.02:0.7;%0.1:0.01:0.7;% for plotting all spectra curves only
freq0=0; % frequency start from 0Hz
ndistri=length(distrifreq);
inddistrifreq=round((distrifreq-freq0)/DFxcorr)+1;
distrixcorr=zeros(niter,ndistri);
distriratio=zeros(niter,ndistri);% FOR THE NORMALIZED XSPECTRA
distrilower=zeros(niter,ndistri);
% wftime(:,1)=wftime(:,1)/sqrt(sum(wftime(:,1).^2));
% wftime(:,2)=wftime(:,2)/sqrt(sum(wftime(:,2).^2));

%% apply onebit normalization

if applyonebitnorm
    
    wftime=onebitnorm(wftime);% apply onebit norm to each column of matrix wftime;
end

% maxlag=100/dt+1;
% plot(xcorr(wftime(:,1),wftime(:,2),maxlag))

% DEFINE TWO POINTERS FOR THE WINDOW
i1=1;
i2=1+SegLenFreq-1;

% NORMALIZE THE TIME DOMAIN NOISE BY TOTAL POWER


% START A WHILE LOOP ON HOPPING WINDOW!
iter=1;% count iterations
alphatukey=0.1;
tukeyw=tukeywin(SegLenFreq,  alphatukey  ); % use tapered tukey window!
tukeyw=1;% use rect window
while i2<rightcorner
    % for iter=1:niter
    
    
    % produce two orthogonal (uncorrelated) white noise RVs:
    %     rv0=(noiseBsqrt+noiseBstdv*randn(nrv,1)).*exp(1i*rand(nrv,1)*2*pi);% col vector!
    %     rv1=(noiseBsqrt+noiseBstdv*randn(nrv,1)).*exp(1i*rand(nrv,1)*2*pi);% col vector!
    %
    %     upper=staBphase*rv0+deltaAtri*(rv1.*Fbool');% right
    %     lower=staAphase*rv0+deltaAtri*(rv1.*~Fbool');% left
    %     upper=staBphase*rv0+deltaB*(rv1.*Fbool');% right
    %     lower=staAphase*rv0+deltaA*(rv1.*~Fbool');% left
    
    % PRODUCE FFTs OF TR1 AND TR2 SEGMENTS
    tr1=zeros(LenFFT,1);tr2=tr1;
    tr1(1:SegLenFreq)=wftime(i1:i2,1).*tukeyw;
    tr2(1:SegLenFreq)=wftime(i1:i2,2).*tukeyw;
    
    lower=fft(tr1);
    upper=fft(tr2);
    
    
    
    xcorrseg=xcorrseg+conj(lower).*upper;
    %xcorrnoise=xcorrnoise+conj(lowernoise).*uppernoise;
    
    numerator=conj(lower).*(upper);% xcorr
    denominator=abs(lower).*abs(upper);
    
    psd1=psd1+conj(lower).*lower;
    psd2=psd2+conj(upper).*upper;
    
    %     psdnoise1=psdnoise1+conj(lowernoise).*lowernoise;
    %     psdnoise2=psdnoise2+conj(uppernoise).*uppernoise;
    %numerator=conj(lower).*upper.*lower.*conj(upper);
    %fourthM=fourthM+numerator;
    ratioalpha=ratioalpha+numerator./(denominator);
    
    distrixcorr(iter,:)=( conj(lower(inddistrifreq)).*upper(inddistrifreq) )';
    distriratio(iter,:)=( numerator(inddistrifreq)./(denominator(inddistrifreq)) )';
    
    % This is not needed here% distrilower(iter,:)= conj(lower(inddistrifreq));
    
    
    iter
    iter=iter+1;
    i1=i1+SegLenFreq+hopsize;
    i2=i2+SegLenFreq+hopsize;
    
end

niter=iter-1;
xcorrseg=xcorrseg/niter;
xcorrnoise=xcorrnoise/niter;
psd1=psd1/niter;
psd2=psd2/niter;

psdnoise1=psdnoise1/niter;
psdnoise2=psdnoise2/niter;
ratioalpha=ratioalpha/niter;

unitpsd=ones(LenFFT,1);

if applyprewhiten
    sigmaTotal=sqrt((unitpsd.*unitpsd)/niter);% sigmaRe^2+sigmaIm^2=Var(Xspectrum)
    sigmaRe=sqrt(0.5)*sigmaTotal;% sigmaRe ~ sigmaIm
    
else
    
    sigmaTotal=sqrt((psd1.*psd2)/niter);% sigmaRe^2+sigmaIm^2=Var(Xspectrum)
    sigmaRe=sqrt(0.5)*sigmaTotal;% sigmaRe ~ sigmaIm
    
end
%% %%% PLOT THE PSD FUNCITONS
% % figure(3)
% % semilogy(freqv,psd1-psdnoise1,'b')
% % hold on
% % semilogy(freqv,psd2-psdnoise2,'r')
% % semilogy(freqv,psdnoise1,'b--')
% % semilogy(freqv,psdnoise2,'r--')
% %
% % hold off
% %
% % %%% PLOT THE PSD FUNCITONS
% % figure(1)
% % % plot(freqv,real(ratioalpha),'b')
% % % hold on
% % % plot(freqv,imag(ratioalpha),'r')
% % semilogy(freqv,abs((real(ratioalpha))),'b')
% %
% % freqmeasure=[0.1 0.450];
% %
% % indfreqmeasure=round((freqmeasure-freqwin(1))/df)+1;
% %
% % [pks,locs]=findpeaks(abs((real(ratioalpha(indfreqmeasure(1):indfreqmeasure(2))))));
% %
% % ratioxcorrauto=xcorr./sqrt(psd1.*psd2);
% % % [pks,locs]=findpeaks(abs((real(ratioxcorrauto(indfreqmeasure(1):indfreqmeasure(2))))));
% %
% %
% % hold on
% % semilogy(freqv(locs),pks,'ro')
% %
% % % semilogy(freqv,imag(ratioalpha),'r')
% %
% %
% % % semilogy(freqv,psdnoise1,'b--')
% % % semilogy(freqv,psdnoise2,'r--')
% % semilogy(freqv,sqrt(c./freqv).*exp(-2*pi*(1/Qin)*x*freqv.*invc),'k--')
% %
% %
% % npks=length(pks);% number of peaks
% % Gmat=[-2*pi*x*freqv(locs).*invc(locs)  ones(npks,1)];
% % ymat=log(pks)-0.5*log(c(locs)./freqv(locs));
% % Qinvvec=Gmat\ymat;
% %
% % 1./Qinvvec
% % semilogy(freqv(locs),exp(ymat),'gs')
% %
% % hold off
% %
% % xlim([0.1 0.6])
% % ylim([ min(pks) 10])

%% THE DISTRIBUTION

% remove the unused zero samples
distrixcorr=distrixcorr(1:niter,:);% the unwhitened xspectra
distriratio=distriratio(1:niter,:);% the prewhitened xspectra


% suppose we do whitening:
if applyprewhiten
    distrixcorr=distriratio;
    xcorrseg=ratioalpha;
end

distrixcorr=conj(distrixcorr);
%
if ploteachwf
    % figure(4)
    ifreqdistri=10;
    % mx=mean(real(distrixcorr(:,ifreqdistri)));
    % stdx=std(real(distrixcorr(:,ifreqdistri)))/sqrt(niter);
    % hist(real(distrixcorr(:,ifreqdistri)),30)
    % % text(4e-4,500,'Frequency=0.14 Hz')
    % legend(sprintf('%.2f Hz',distrifreq(ifreqdistri)))
    % legend boxoff
    % title('real part')
    % % axis off
    % set(gca,'Fontsize',fontsize+12)
    % set(gcf,'PaperPositionMode','auto');
    %
    % xlabel('cross-spectrum sample value')
    % title({sprintf('Real part: %d samples.Freq=%.2fHz',niter,distrifreq(ifreqdistri));sprintf('Mean=%.2e; Stde=%.2e',mx,stdx)})
    % set(gca,'Fontsize',fontsize)
    % set(gcf,'PaperPositionMode','auto');
    % saveas(gcf,sprintf('1Drealdistri_%.2f Hz_%d_iters%s.pdf',distrifreq(ifreqdistri),niter,preprocstr),'pdf')
    
    figure(3)
    mx=mean(imag(distrixcorr(:,ifreqdistri)));
    stdx=std(imag(distrixcorr(:,ifreqdistri)));
    hist(imag(distrixcorr(:,ifreqdistri)),30);
    title(sprintf('Imaginary part: %d samples.Freq=%.2fHz Mean=%.2e; Stde=%.2e',niter,distrifreq(ifreqdistri),mx,stdx))
    xlabel('cross-spectrum sample value')
    % fourthM=fourthM/niter;
    
    
    % PLOT THE 2D DISTRIBUTION
    % ifreqdistri=5;
    figure(3)
    hist3([real(distrixcorr(:,ifreqdistri)) imag(distrixcorr(:,ifreqdistri))],[50 50 ]);
    xlabel('Real'); ylabel('Imaginary');
    set(gcf,'renderer','opengl');
    set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
    view(2);view(3);
    % legend(sprintf('%.2f Hz',distrifreq(ifreqdistri)))
    % legend boxoff
    title('Joint')
    set(gca,'Fontsize',fontsize)
    set(gcf,'PaperPositionMode','auto');
    
    title({sprintf('Joint: Freq=%.2f Hz',distrifreq(ifreqdistri)); sprintf(' %d observations/%.2f days.',niter,ndays)} )
    % covariance matrix
    
    %    saveas(gcf,sprintf('2Ddistri_%.2f Hz_%d_iters%s.pdf',distrifreq(ifreqdistri),niter,preprocstr),'pdf')
    
    %% PLOT THE RESULTS
    errvec=std(real(distrixcorr))/sqrt(niter);% real; A ROW VECTOR FOR EACH FREQUENCY
    errvecimag=std(imag(distrixcorr))/sqrt(niter);% imag; A ROW VECTOR FOR EACH FREQUENCY
    errvecboth=errvec+1j*errvecimag;
    
    meanvec=mean(distrixcorr);% A ROW VECTOR FOR EACH FREQUENCY
    % xcorrsmooth=medfilt1(xcorr,5);
    % xcorrsmooth=smooth(xcorrseg,3);
    
    h1=figure(2);
    clf(h1)
    % subplot(211)
    normxspec=false;
    if normxspec
        normfactor=smooth( real(errvecboth),5 ); % normalize the spectra by smoothed stderror
        freqerrornorm = errvecboth.' ./ normfactor;
        
        % plot(freqv,fourthM,'r'); hold on;
        plot(freqv,real(xcorrseg),'b'); hold on
        plot(freqv,imag(xcorrseg),'r-')
        errorbar(distrifreq,real(meanvec),errvec,'.')
        errorbar(distrifreq,imag(meanvec),errvecimag,'.')
        
    else
        % plot(freqv,fourthM,'r'); hold on;
        plot(freqv,real(xcorrseg),'b'); hold on
        plot(freqv,imag(xcorrseg),'r-')
        errorbar(distrifreq,real(meanvec),errvec,'.')
        errorbar(distrifreq,imag(meanvec),errvecimag,'.')
        
        
    end
    %  plot(freqv,psd1.*psd1,'y--')
    %  plot(freqv,psd2.*psd2,'g--')
    %plot(freqv,conj(xcorr).*xcorr-((fourthM)-psd1.*psd2),'k')
    %  plot(freqv,((fourthM)-psd1.*psd2)-conj(xcorr).*xcorr,'k')
    % title({'cross-spectra';sprintf('%d observations/%.2f days',niter,ndays)})
    title(sprintf('cross-spectra %d observations/%.2f days',niter,ndays))
    ylabel('Cross-spectra')
    xlabel('Frequency (Hz)')
    % plot(freqv,sqrt((psd1-psdnoise1).^2/niter),'r--')
    
    %%%%% THE LINE BELOW IS USED TO APPROXIMATE THE THEORETICAL ERROR
    %%%%% Var(Xspectrum)==sigRe^2+sigIm^2
    % plot(freqv,sqrt((psd1.*psd2)/niter),'r--')
    
    
    plot([min(freqv) max(freqv)],[0 0],'k--')
    
    
    % plot(freqv,real(xcorrsmooth),'k')
    hold off
    % legend('Real part of averaged cross-spectra','Imaginary part of averaged cross-spectra','Error bar for real part','Error on amplitude overall')
    % legend('Real part of averaged cross-spectra','Imaginary part of averaged cross-spectra','Error bar for real part','Error bar for imag part')
    legend('Real part of averaged cross-spectra','Imaginary part of averaged cross-spectra')
    
    xlim([0.1 0.6])
    set(gca,'Fontsize',fontsize-2)
    set(gcf,'PaperPositionMode','auto');
    %  saveas(gcf,sprintf('NewCrossspectra_%d_iterations%s.eps',niter,preprocstr),'psc2')
    %  saveas(gcf,sprintf('NewCrossspectra_%d_iterations%s.pdf',niter,preprocstr),'pdf')
end