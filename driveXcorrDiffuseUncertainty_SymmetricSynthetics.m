%%%%% COMPUTE CROSS-CORRELATION OF
%%%%% Diffuse noise
%%%%% NEW ADDONS: 
%%%%%       
%%%%% Xin Liu, Stanford Univ, 2019
%%%%% ALRIGHTS RESERVED 
clear

sim4zerofreq = false; false;%false;%true;% false; % if true, store Xspectra for zero frequency
symmCompOnly=false;


ploteachwf=0;% SWITCH FOR MAKING PLOTS
% clear

AKINORM=1;0; % NORMALIZE THE X-SPECTRA BY POWER SPECTRA BEFORE GETTING PHASE V
fontsize=18;

initParams
% srcNdispersion

% nsta=6;
% nsta=4;
%  CA STATIONS AND 1D DISP: continous nondiffuse src
% wfdatatimeSyndiffuseAniso_HHZ_15
% fname=sprintf('~/scwork/wfdatatimeCAdiffuse_%2.1f_days_%d_new.mat',ndays,nsta);
fname = 'twostaWF.mat';
dataset=load(fname);

wftime= dataset.wf5days;
%%  THE FULLY DIFFUSE NOISE::CONVERT TO PETERSON NOISE SPECTRA ENVELOPE

% MODIFY THE SPECTRAL AMPLITUDE (WHITEN, THEN APPLY Peterson Low Noise Model)


xlimrange=[-100 100];
% xlimrange=[-100 100]*10;
%% MAKE PLOTS::

% COAST PERPENDICULAR
sta2=1;sta1=3; % CHF-SBB2

% GET STATION STRINGS and DISTANCE
ssta2=stalist{sta2};
ssta1=stalist{sta1};

dist=1e-3*norm(coorddict(ssta1)-coorddict(ssta2));% CONVERT TO km


% MIDDLE POINT OF RAW XCORR:
% midpX=nptsS;

%% XCORR SLIDE WINDOW:
ndayselect= 5; %0.25 %0.5 %0.25 %1%15;%5;
preprocstr = '';
% preprocstr = 'prewhiten';
% preprocstr = 'onebit';
if sim4zerofreq
    freqwinstore=[0,0.6];
else
    freqwinstore =  [0.10,0.6];
    freqwinstore =  [0.01,0.6];
end
winlen = 100;%400;%100;%200; % 200 sec window
DF=0.5/winlen;
% wftime=wftimestaPeterson(:,[sta1, sta2]);
[distrixcorr,errvecboth,freqrange]=createXspecMat(wftime,dt,winlen,freqwinstore,preprocstr,ndayselect);
FreqList=round(freqrange/DF);
normfactor=smooth( real(errvecboth),5 ); % normalize the spectra by smoothed stderror
freqerrornorm = errvecboth.' ./ normfactor;
%% BOOTSTRAP 
keypair=sprintf('%s-%s',sta1,sta2);
Nbootstrap=2000;
% dayrange=[1,ndayselect];
getStackedWins=false; % 'false' means get bootstrap
bootfreqarray=bootstrapblocknday(distrixcorr,ndayselect,freqrange,Nbootstrap,getStackedWins,keypair);
if symmCompOnly
    normbootarray=bsxfun(@rdivide,real(bootfreqarray), normfactor' );
else
    normbootarray=bsxfun(@rdivide,(bootfreqarray), normfactor' );    
end
% make vfreq=0 in order to use FIR band-pass filter for the entire range:
vfreq=0.25;%0.25;%0; % center frequency for narrow band filter; 
vfreq=0;

[xcorrarray,xdigit,fltrcoef,subindspec]=getBootXcorrFiltered(normbootarray,freqrange,vfreq);
npts = size(xcorrarray,2);
midp=(npts+1)/2;
[newDenMat,PSDcenters]=assemDensityImageXcorr(real(xcorrarray));

% makeMovieRandomFluct(xcorrarray,xdigit)

% % TRY C3 COMPUTATIONS:
% calC3=false;%false;
% percExclT=0.4;% exclude how much percentage of Xcorr time series data
% if vfreq==0 && calC3
%     CalC3Boots(real(xcorrarray),midp,dt,percExclT)
% end

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
plot(xdigit,mean(real(xcorrarray)) , 'r')

hold off
colormap(flipud(gray))
caxis([0,320])
% caxis([0,100])
colorbar
if vfreq~=0
    legend(sprintf('f_c=%.2fHz',vfreq))
end
set(gca,'YDir','normal')

xlim([-100 100])
title( sprintf('stacked cross-correlation %d days',ndayselect)  )
%         set(gca,'Fontsize',fontsize)
%         set(gcf,'PaperPositionMode','auto');
if vfreq==0
   xlabel('time (s)')
    meanxcorrnorm=mean(real(xcorrarray));
   save('xcorr_stacked_C1.mat','meanxcorrnorm','xdigit') 
end

    if vfreq~=0
        [newDenMatXcorrEnv,PSDcenters]=assemDensityImageXcorr(abs(xcorrarray));
        subplot(212)
        imagesc(xdigit,PSDcenters,newDenMatXcorrEnv)
        caxis([0  max(max(newDenMatXcorrEnv))*0.3 ])
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


%% STORE KEY SPECTRAL DATA FOR FURTHER OPTIMIZATION

meanspec= mean(bootfreqarray);
stderr = errvecboth;
distance=dist*1e3;% convert to meter

enddaysamp=0;
fnamesave=sprintf('syn_%s_%s_Xspecdata_stacked_%ddays_winlen_%ds.mat',ssta2,ssta1,ndayselect,winlen);

% USE FRACTION NUMBER of ndayselect:
if ~sim4zerofreq
    fnamesave=sprintf('syn_%s_%s_Xspecdata_stacked_%.2fdays_winlen_%ds.mat',ssta2,ssta1,ndayselect,winlen);
else
    fnamesave=sprintf('syn_%s_%s_Xspecdata_stacked_%.2fdays_winlen_%ds_zerofreq.mat',ssta2,ssta1,ndayselect,winlen);
end

save(fnamesave,'DF','dt','FreqList','stderrTD','distance','stderr','meanspec','bootfreqarray','enddaysamp')

aa=5;
% end
