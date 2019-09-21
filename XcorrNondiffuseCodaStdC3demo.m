%%%%% COMPUTE CROSS-CORRELATION OF
%%%%% NON-DIFFUSE COMPONENTS added to FULLY DIFFUSE NOISE FIELDS
%%%%% NEW ADDONS: 
%%%%% Xin Liu, Stanford Univ, 2017
%%%%% ALRIGHTS RESERVED Hahah
clear

symmCompOnly=false;
ploteachwf=1;% SWITCH FOR MAKING PLOTS

AKINORM=0; % NORMALIZE THE X-SPECTRA BY POWER SPECTRA BEFORE GETTING PHASE V
fontsize=18;

initParams
srcNdispersion

nsta=6;%4;
%  CA STATIONS AND 1D DISP: continous nondiffuse src
% wfdatatimeSyndiffuseAniso_HHZ_15
fname=sprintf('C:/Users/liusi/OneDrive/wf_sim/wfdatatimeSyndiffuseAniso_HHZ_%2.1f_days_%d_new.mat',ndays,nsta);
fname=sprintf('C:/Users/liusi/OneDrive/wf_sim/wfdatatimeSyndiffuseAnisoCA_HHZ_%2.1f_days_%d_new.mat',ndays,nsta);
% fname=sprintf('~/scwork/wfdatatimeSyndiffuseAniso_HHZ_%2.1f_days_%d_new.mat',ndays,nsta);

%  CA STATIONS AND 1D DISP::for C3:: continous nondiffuse src
% fname=sprintf('~/scwork/wfdatatimeCANondiff_%2.1f_days_%d.mat',ndays,nsta);
% fname=sprintf('C:/Users/liusi/scwork/wfdatatimeCANondiff_%2.1f_days_%d.mat',ndays,nsta);

% C:\Users\liusi\scwork

%  CA STATIONS AND 1D DISP: 5 discrete nondiffuse src
%  fname=sprintf('~/scwork/wfdatatimeCA_5srcNondiff_%2.1f_days_%d.mat',ndays,nsta);


% % SYNTHETIC DATA & SYNTHETIC DISPERSION
% fname=sprintf('~/scwork/wfdatatimeNondiff_%2.1f_days_%d.mat',ndays,nsta);
% % % % % % % fname=sprintf('C:/Users/liusi/scwork/wfdatatime_%2.1f_days_%d.mat',ndays,nsta);
load(fname)

%%  THE FULLY DIFFUSE NOISE::CONVERT TO PETERSON NOISE SPECTRA ENVELOPE

% MODIFY THE SPECTRAL AMPLITUDE (WHITEN, THEN APPLY Peterson Low Noise Model)

wftimestaPeterson=zeros(nptsS,nsta);

for ista=1:nsta
%     wftimestaPeterson(:,ista)=conv(wftimesta(:,ista),XfreqTimeCentreScaling,'same');
    wftimestaPeterson(:,ista)=wftimesta(:,ista);
end


xlimrange=[-100 100];
xlimrange=[-60 60];
% xlimrange=[-100 100]*10;
%% MAKE PLOTS::

% COAST PARALLEL
sta2=1;sta1=6;
sta2=3;sta1=5;
sta2=3;sta1=4;

% COAST PERPENDICULAR
sta2=1;sta1=2;
sta2=1;sta1=3;
% sta2=3;sta1=2;

% BELOW DEFINE THE VIRTUAL SRC STA:
stav=6; % for CA synthetic data 6 sta;
% stav=4; % synthetic data: 4 stations

% GET STATION STRINGS and DISTANCE
ssta2=stalist{sta2};
ssta1=stalist{sta1};
sstav=stalist{stav};

dist=1e-3*norm(coorddict(ssta1)-coorddict(ssta2));% CONVERT TO km


% MIDDLE POINT OF RAW XCORR:
midpX=nptsS;

%% XCORR SLIDE WINDOW:
ndayselect=5;%5;
preprocstr = '';
% preprocstr = 'prewhiten';
% preprocstr = 'onebit';

freqwinstore=[0.10,0.6];
winlen = 100;%400;%100;%200; % 200 sec window
wftime1=wftimestaPeterson(:,[stav, sta1]);
wftime2=wftimestaPeterson(:,[stav, sta2]);
[distrixcorr1,errvecboth1,freqrange]=createXspecMat(wftime1,dt,winlen,freqwinstore,preprocstr,ndayselect);
normfactor1=smooth( real(errvecboth1),5 ); % normalize the spectra by smoothed stderror
freqerrornorm1 = errvecboth1.' ./ normfactor1;

[distrixcorr2,errvecboth2,freqrange]=createXspecMat(wftime2,dt,winlen,freqwinstore,preprocstr,ndayselect);
normfactor2=smooth( real(errvecboth2),5 ); % normalize the spectra by smoothed stderror
freqerrornorm2 = errvecboth2.' ./ normfactor2;


%% BOOTSTRAP 
getStackedWins=true;

keypair=sprintf('%s-%s',sta1,sta2);
Nbootstrap=2000;
dayrange=[1,ndayselect];
bootfreqarray1=bootstrapblock(distrixcorr1,dayrange,freqrange,Nbootstrap,getStackedWins,keypair);
if symmCompOnly
    normbootarray1=bsxfun(@rdivide,real(bootfreqarray1), normfactor1' );
else
    normbootarray1=bsxfun(@rdivide,(bootfreqarray1), normfactor1' );    
end
bootfreqarray2=bootstrapblock(distrixcorr2,dayrange,freqrange,Nbootstrap,getStackedWins,keypair);
if symmCompOnly
    normbootarray2=bsxfun(@rdivide,real(bootfreqarray2), normfactor2' );
else
    normbootarray2=bsxfun(@rdivide,(bootfreqarray2), normfactor2' );    
end
% make vfreq=0 in order to use FIR band-pass filter for the entire range:
vfreq=0;%0.25;%0; % center frequency for narrow band filter; 

[xcorrarray1,xdigit,fltrcoef,subindspec]=getBootXcorrFiltered(normbootarray1,freqrange,vfreq);
dt=xdigit(2)-xdigit(1);
npts = size(xcorrarray1,2);
midp=(npts+1)/2;
[newDenMat,PSDcenters]=assemDensityImageXcorr(real(xcorrarray1));

[xcorrarray2,xdigit,fltrcoef,subindspec]=getBootXcorrFiltered(normbootarray2,freqrange,vfreq);

calC3=true;
percExclT=0.4;% exclude how much percentage of Xcorr time series data
if vfreq==0 && calC3
    [CCpos,CCneg,xdigitC3]=CalC3Boots(real(xcorrarray2),real(xcorrarray1),midp,dt,percExclT);
end



load('xcorr_stacked_C1')

coefFixAmp=max(meanxcorrnorm)/max(CCpos);
figure(2)
plot(xdigit,meanxcorrnorm,'k')
hold on

plot(xdigitC3,CCpos*coefFixAmp,'-.')
plot(xdigitC3,CCneg*coefFixAmp,'-.')


hold off
xlabel('time (s)')
xlim([-100,100])
% ylim([-1.5, 1])
legend('C1','Xcorr positive coda','Xcorr negative coda')
title('C3 from 120 sec block')
title('C3 from 240 sec block')
title('C3 from 720 sec block')
title('C3 from 1800 sec block')
if symmCompOnly
    stderrTD=estiAmpError(fltrcoef,real(freqerrornorm),subindspec,npts);
else
    stderrTD=estiAmpError(fltrcoef,(freqerrornorm),subindspec,npts);    
end

stderrTD=[stderrTD(midp+1:end); stderrTD(1:midp) ];

figure(5)
clf
if vfreq~=0
   subplot(211) 
end

imagesc(xdigit,PSDcenters,newDenMat)
hold on
plot(xdigit,mean(real(xcorrarray1)) , 'r')

hold off
colormap(flipud(gray))
caxis([0,320])
colorbar
if vfreq~=0
    legend(sprintf('f_c=%.2fHz',vfreq))
end
set(gca,'YDir','normal')

xlim([-100 100])
xlabel('time (s)')
title( sprintf('stacked cross-correlation %d days',ndayselect)  )
%         set(gca,'Fontsize',fontsize)
%         set(gcf,'PaperPositionMode','auto');


    if vfreq~=0
        [newDenMatXcorrEnv,PSDcenters]=assemDensityImageXcorr(abs(xcorrarray));
        subplot(212)
        imagesc(xdigit,PSDcenters,newDenMatXcorrEnv)
        caxis([0  max(max(newDenMatXcorrEnv))*0.25 ])
        xlabel('Time (sec)')
        %legend('real part Xspectra','imag part Xspectra')
        title(sprintf('bootstrap stacked envelope %d days',ndayselect))
        
        
        colormap(flipud(gray))
        %colorarr=[];
        set(gca,'YDir','normal')
        colorbar
        hold on
        plot(xdigit,mean(abs(xcorrarray)),'r-')
        if vfreq
%             legend(sprintf('%s-%s f_c=%.2fHz',sta1,sta2,vfreq))
            legend(sprintf('f_c=%.2fHz',vfreq))
        else
            legend(sprintf('%s-%s',sta1,sta2))
        end
        %     plot(freqrange,log10(errpercentile(1,:) ))
        %     plot(freqrange,log10(errpercentile(2,:) ))
        hold off
        %         pbaspect([2.5,1,1])
        xlim([-100 100])
%                 xlim([-50 50])
%         set(gca,'Fontsize',fontsize)
%         set(gcf,'PaperPositionMode','auto');
    end



figure(6)
subplot(211)
plot(xdigit,std(real(xcorrarray)))
hold on
plot(xdigit,std(imag(xcorrarray)),'r')
plot(xdigit,stderrTD,'k')
hold off
xlabel('time (s)')
ylabel('std of wf')
legend('real part error','imag part error','stderr prediction')
title(sprintf('std for each time point: winLen=%d s',winlen))
subplot(212)
% plot(xdigit,(abs(xcorrarray))./std(real(xcorrarray)))
plot(xdigit,abs(mean(xcorrarray))./std(real(xcorrarray)))
xlabel('time (s)')
ylabel('SNR of wf envelope')
title('SNR for each time point')

%% XCORR THE WHOLE TRACE:

codawin=[100,lentracetotal/2];
codaindp=round(codawin/dt)+midpX;
codaindn=round(-codawin(end:-1:1)/dt)+midpX;

xcorrdiff=xcorr(wftimestaPeterson(:,sta1),wftimestaPeterson(:,sta2));
xcorrnon=xcorr(wftimestanon(:,sta1)+wftimestaPeterson(:,sta1),wftimestanon(:,sta2)+wftimestaPeterson(:,sta2));
% nptswf=length(wftimesta(:,1));
% timedigitsXcorr=(-nptswf+1:nptswf-1)*dt;
% if ploteachwf & 1
%     
%     figure(8)
%     subplot(211)
%     %     plot(timedigitsXcorr,xcorr(wftimesta(:,sta1),wftimesta(:,sta2)))
%     plot(timedigitsXcorr,xcorrdiff)
%     xlim(xlimrange)
%         xlabel('time (sec)')
%    title(sprintf('Xcorr of fully diffuse noise for %s-%s',ssta2,ssta1))
%    set(gca,'FontSize',fontsize)
%    set(gcf,'PaperPositionMode','auto');
%     subplot(212)
%     % ONLY NON=DIFFUSE NOISE
%     %     plot(timedigitsXcorr,xcorr(wftimestanon(:,sta1),wftimestanon(:,sta2)))
%     % XCORR of BOTH DIFF AND NON-DIFF
%     plot(timedigitsXcorr,xcorrnon)
%     
%     xlabel('time (sec)')
%    title(sprintf('Xcorr with non-diffuse noise for %s-%s',ssta2,ssta1))
%    set(gca,'FontSize',fontsize)
%    set(gcf,'PaperPositionMode','auto');
%        xlim(xlimrange)
% 
% %    saveas(gcf,sprintf('XcorrCompareNondiffsim_%s_%s_trueconti.pdf',ssta2,ssta1),'pdf')
% 
% end
% 
% 
% %% NOW TRY TO EXTRACT PHASE / GROUP VELOCITIES!
% 
% 
% freqwinsub=[0.1 0.6];
% 
% lenxcorrT=250; % in sec
% window=[-lenxcorrT/2 lenxcorrT/2];
% nptsW=lenxcorrT/dt+1;% WINDOWED XCORR LENGTH
% midpW=(nptsW+1)/2;% middle point of windowed xcorr
% dfxcorr=1/lenxcorrT;% freq spacing for windowed xcorr
% freqrange=freqwinsub(1):dfxcorr:freqwinsub(2);
% freqvPhasev=freqwinsub(1):dfxcorr:freqwinsub(2);
% 
% 
% if AKINORM
%     
%     % ESTIMATE POWER SPECTRA
%     autodiff=xcorr(wftimestanon(:,sta1)+wftimestaPeterson(:,sta1),...
%         wftimestanon(:,sta1)+wftimestaPeterson(:,sta1));
%     autocorrW=timewin(autodiff,window,dt,midpX);
%     psdP=abs( fft(autocorrW, round(nptsW)) );
%     psdP=smooth(psdP,20);
%     % lenPeterson=length(XfreqTimeCentreScaling);
%     % autoP=zeros(nptsW,1);
%     % autoP(1:(lenPeterson+1)/2)=XfreqTimeCentreScaling((lenPeterson+1)/2:end);
%     % autoP(midpW+1:end)=autoP(midpW:-1:2);
%     % psdP=abs( fft(autoP, round(nptsW)) );
% end
% taperwinW=tukeywin(nptsW,0.1);
% xcorrdiffW=timewin(xcorrdiff,window,dt,midpX);
% xcorrdiffW=shift0(xcorrdiffW,midpW).*taperwinW;
% spectra=fft(xcorrdiffW);
% if AKINORM
%     spectra=spectra./psdP;
% end
% 
% [slowgrdiff,slowphTDdiff,freqmodifiedSYM]= ... 
%     groupvelsym(real(spectra),dfxcorr,freqwinsub,dist);
% 
% xcorrnonW=timewin(xcorrnon,window,dt,midpX);
% xcorrnonW=shift0(xcorrnonW,midpW).*taperwinW;
% 
% spectranon=fft(xcorrnonW);
% 
% if AKINORM
%     spectranon=spectranon./psdP;
% end
% 
% [slowgrnon,slowphTDnon,~]= ... 
%     groupvelsym(real(spectranon),dfxcorr,freqwinsub,dist);
% alambda=1.5;% x times wavelength
% lowfreqbound=min( freqvPhasev(freqvPhasev'>alambda./slowphTDdiff/dist) );
% freqxlim=[max([0.12,lowfreqbound]) freqwinsub(2)];
% % COMPARE DIFFUSE AND NONDIFFUSE PHASE VELOCITY
% linewidth=1.5;
% if ploteachwf
%    figure(5)
%    clf
%    plot(freqvPhasev,1./slowphTDdiff,'LineWidth',linewidth)
%     
%    hold on
%    
%    plot(freqvPhasev,1./slowphTDnon,'-.','LineWidth',linewidth)
%    
%    hold off
%    xlim(freqxlim)
%    xlabel('frequency (Hz)')
%    ylabel('phase velocity (km/sec)')
%    title(sprintf('Station pair %s-%s',ssta2,ssta1))
%    legend('fully diffuse','added non-diffuse noise')
%     set(gca,'FontSize',fontsize)
%     set(gcf,'PaperPositionMode','auto');
%     saveas(gcf,sprintf('PhaseVelCompareNondiffsim_%s_%s_trueconti.pdf',ssta2,ssta1),'pdf')
% 
% 
%     % ESTIMATE dv/v velocity bias
%     dv_v=slowphTDdiff./slowphTDnon-1;
%     figure(6);plot(freqvPhasev,dv_v*100)
%       xlim(freqxlim)
%    xlabel('frequency (Hz)')
%    ylabel('phase velocity relative error (%)')
%    title(sprintf('Relative phase velocity bias for %s-%s',ssta2,ssta1))
%    set(gca,'FontSize',fontsize)
%    set(gcf,'PaperPositionMode','auto');
%    save(sprintf('dataphasevel_%s-%s_trueconti.mat',ssta2,ssta1), ...
%        'freqvPhasev','dv_v','slowphTDnon','slowphTDdiff')
% end
% 
% % plot(freqvDisp,c) % true dispersion curve from model
% %% TRY C3
% if ploteachwf
%     stav=2%1;
%     coda1=xcorr(wftimestanon(:,sta1)+wftimestaPeterson(:,sta1),wftimestanon(:,stav)+wftimestaPeterson(:,stav),'unbiased');
%     coda2=xcorr(wftimestanon(:,sta2)+wftimestaPeterson(:,sta2),wftimestanon(:,stav)+wftimestaPeterson(:,stav),'unbiased');
%     
% %     coda1=xcorr(wftimestaPeterson(:,stav),wftimestaPeterson(:,sta1),'unbiased');   
% %     coda2=xcorr(wftimestaPeterson(:,stav),wftimestaPeterson(:,sta2),'unbiased');
% 
%     C3p=xcorr(coda1(codaindp(1):codaindp(2)), coda2(codaindp(1):codaindp(2)) );
%     C3n=xcorr(coda1(codaindn(1):codaindn(2)), coda2(codaindn(1):codaindn(2)) );
%     
%     nptsC3=(length(C3p)+1)/2;
%     timedigitsC3=(-nptsC3+1:nptsC3-1)*dt;
%     
%     
%     
%     figure(8)
%     subplot(211)
%     %     plot(timedigitsXcorr,xcorr(wftimesta(:,sta1),wftimesta(:,sta2)))
%     plot(timedigitsXcorr,xcorrdiff)
%     xlim(xlimrange)
%     subplot(212)
%     % ONLY NON=DIFFUSE NOISE
%     %     plot(timedigitsXcorr,xcorr(wftimestanon(:,sta1),wftimestanon(:,sta2)))
%     % XCORR of BOTH DIFF AND NON-DIFF
%     plot(timedigitsC3,C3p+C3n)
%     
%     xlim(xlimrange)
% end
% 
% 
% %% PLOT NOISE RECORD
% if ploteachwf
%     
%     figure(9)
%     subplot(211)
%     %     plot(timedigitsXcorr,xcorr(wftimesta(:,sta1),wftimesta(:,sta2)))
%     plot(timedigitsS,wftimestaPeterson(:,sta1))
%     title(sprintf('Fully diffuse noise for %s',ssta2))
%    set(gca,'FontSize',fontsize)
%    set(gcf,'PaperPositionMode','auto');
%          xlim([0 1000]+5000)
%     subplot(212)
%     plot(timedigitsS,wftimestanon(:,sta1))
%         title(sprintf('Non-diffuse noise for %s',ssta2))
%    set(gca,'FontSize',fontsize)
%    set(gcf,'PaperPositionMode','auto');
%          xlim([0 1000]+5000)
%          xlabel('time (s)')
%          ylim([-0.05 0.05])
% end