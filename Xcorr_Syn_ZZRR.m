%%%%% COMPUTE CROSS-CORRELATION OF
%%%%% NON-DIFFUSE COMPONENTS added to FULLY DIFFUSE NOISE FIELDS
%%%%% Xin Liu, Stanford Univ, 2017
%%%%% ALRIGHTS RESERVED Hahah
clear
 
israyleigh= true; %false;% true; if==true, Rayleigh nondiffuse
                        % if==false, Love nondiffuse

winlenhalf=100;% in sec; cut the long Xcorr into short lengths of winlenhalf*2

ploteachwf=1;% SWITCH FOR MAKING PLOTS

AKINORM=1; % NORMALIZE THE X-SPECTRA BY POWER SPECTRA BEFORE GETTING PHASE V

initParams
srcNdispersion
fontsize=12;%18;

%  CA STATIONS AND 1D DISP: continous nondiffuse src
% fname=sprintf('~/scwork/wfdatatimeCANondiff_%2.1f_days_%d.mat',ndays,nsta);
% fname=sprintf('C:/Users/liusi/scwork/wfdatatimeCANondiff_%2.1f_days_%d.mat',ndays,nsta);
% 
%  CA STATIONS AND 1D DISP: 5 discrete nondiffuse src
%  fname=sprintf('~/scwork/wfdatatimeCA_5srcNondiff_%2.1f_days_%d.mat',ndays,nsta);


% % SYNTHETIC DATA & SYNTHETIC DISPERSION
nsrcnon=1;
degsrc=220; %200;%210;%220;%250;%270;
% fname=sprintf('~/scwork/wfdatatimeNondiff_%2.1f_days_%d.mat',ndays,nsta);
if ~israyleigh
    fname=sprintf('~/scwork/dirwfdata/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHTprojZ_%2.1f_days_%d.mat',degsrc(1),nsrcnon,ndays,nsta);
else
    fname=sprintf('~/scwork/dirwfdata/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHZ_%2.1f_days_%d.mat',degsrc,nsrcnon,ndays,nsta);
end

% fname=sprintf('~/scwork/dirwfdata/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHZ_%2.1f_days_%d.mat',degsrc,nsrcnon,ndays,nsta);

% % % % % % % fname=sprintf('C:/Users/liusi/scwork/wfdatatime_%2.1f_days_%d.mat',ndays,nsta);
load(fname)

if israyleigh
    fname=sprintf('~/scwork/dirwfdata/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHE_%2.1f_days_%d.mat',degsrc,nsrcnon,ndays,nsta);
else
    fname=sprintf('~/scwork/dirwfdata/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHTprojR_%2.1f_days_%d.mat',degsrc,nsrcnon,ndays,nsta);
end
% % % % % % % fname=sprintf('C:/Users/liusi/scwork/wfdatatime_%2.1f_days_%d.mat',ndays,nsta);
datasetE=load(fname);

%%  THE FULLY DIFFUSE NOISE::CONVERT TO PETERSON NOISE SPECTRA ENVELOPE

% MODIFY THE SPECTRAL AMPLITUDE (WHITEN, THEN APPLY Peterson Low Noise Model)

wftimestaPeterson=zeros(nptsS,nsta);

for ista=1:nsta
    wftimestaPeterson(:,ista)=conv(wftimesta(:,ista),XfreqTimeCentreScaling,'same');
end

% % % % % % % wftimestaPeterson=wftimesta;
wftimestaPetersonE=zeros(nptsS,nsta);

for ista=1:nsta
    wftimestaPetersonE(:,ista)=conv(datasetE.wftimesta(:,ista),datasetE.XfreqTimeCentreScaling,'same');
end





xlimrange=[-100 100];
% xlimrange=[-100 100]*10;
%% MAKE PLOTS::

% COAST PARALLEL
sta2=1;sta1=6;
 sta2=3;sta1=5;
% sta2=3;sta1=4;

% COAST PERPENDICULAR
sta2=1;sta1=2;
% sta2=1;sta1=3;
% sta2=3;sta1=2;

% GET STATION STRINGS and DISTANCE
ssta2=stalist{sta2};
ssta1=stalist{sta1};

dist=1e-3*norm(coorddict(ssta1)-coorddict(ssta2));% CONVERT TO km


% MIDDLE POINT OF RAW XCORR:
midpX=nptsS;

codawin=[100,lentracetotal/2];
codaindp=round(codawin/dt)+midpX;
codaindn=round(-codawin(end:-1:1)/dt)+midpX;

xcorrdiff=xcorr(wftimestaPeterson(:,sta1),wftimestaPeterson(:,sta2));
xcorrnon=xcorr(wftimestanon(:,sta1)+wftimestaPeterson(:,sta1),wftimestanon(:,sta2)+wftimestaPeterson(:,sta2));

xcorrdiffE=xcorr(wftimestaPetersonE(:,sta1),wftimestaPetersonE(:,sta2));
xcorrnonE=xcorr(datasetE.wftimestanon(:,sta1)+wftimestaPetersonE(:,sta1), ...
                datasetE.wftimestanon(:,sta2)+wftimestaPetersonE(:,sta2));

if AKINORM
    xcorrdiff=conv(xcorrdiff,XfreqTimewhitencenter,'same');
    xcorrnon=conv(xcorrnon,XfreqTimewhitencenter,'same');
    xcorrdiffE=conv(xcorrdiffE,XfreqTimewhitencenter,'same');
    xcorrnonE=conv(xcorrnonE,XfreqTimewhitencenter,'same');
else
    maxWhitenWavelet=1.5/max(xcorrdiff);
    xcorrdiff=xcorrdiff*maxWhitenWavelet;
    xcorrnon=xcorrnon*maxWhitenWavelet;
    xcorrdiffE=xcorrdiffE*maxWhitenWavelet;
    xcorrnonE=xcorrnonE*maxWhitenWavelet;
end
            
clear datasetE
nptswf=length(wftimesta(:,1));
timedigitsXcorr=(-nptswf+1:nptswf-1)*dt;

%% STORE THE DATA:
timedigitsXcorrshort=shortwindata(timedigitsXcorr,winlenhalf,dt)';
xcorrdiffshort=shortwindata(xcorrdiff,winlenhalf,dt);
xcorrnonshort=shortwindata(xcorrnon,winlenhalf,dt);
xcorrdiffEshort=shortwindata(xcorrdiffE,winlenhalf,dt);
xcorrnonEshort=shortwindata(xcorrnonE,winlenhalf,dt);

% TO STORE: nptsshort*2 matrix; col index 1==Z comp & 2==E comp
xcorr3Cdiff=[xcorrdiffshort,xcorrdiffEshort];
xcorr3C=[xcorrnonshort,xcorrnonEshort];
if ~israyleigh
    lovewstr='Love';
else
    lovewstr='';
end
if AKINORM
    fnamesavedata=sprintf('wfXcorr_synthetics_nondiff_AKINORM_2sta_%dsrc_%ddeg%s.mat',nsrcnon,degsrc(1),lovewstr);
else
    fnamesavedata=sprintf('wfXcorr_synthetics_nondiff_2sta_%dsrc_%ddeg%s.mat',nsrcnon,degsrc(1),lovewstr);
end
distance=dist*1e3;% convert dist to meters
% USEFUL DATA; BUT DON'T SAVE YET... % save(fnamesavedata,'dt','timedigitsXcorrshort','xcorr3Cdiff','xcorr3C','distance')
%% MAKE PLOTS:
if ploteachwf & 1
    figure(8)
    subplot(211)
    %     plot(timedigitsXcorr,xcorr(wftimesta(:,sta1),wftimesta(:,sta2)))
    plot(timedigitsXcorr,xcorrdiff)
    xlim(xlimrange)
    xlabel('time (sec)')
    title(sprintf('Xcorr of fully diffuse noise for %s-%s',ssta2,ssta1))
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');
    subplot(212)
    % ONLY NON=DIFFUSE NOISE
% % % % %         plot(timedigitsXcorr,xcorr(wftimestanon(:,sta1),wftimestanon(:,sta2)))
        plot(timedigitsXcorr, xcorrnon-xcorrdiff)
    % XCORR of BOTH DIFF AND NON-DIFF
%     plot(timedigitsXcorr,xcorrnon)
    
    xlabel('time (sec)')
    title(sprintf('Xcorr with non-diffuse noise for %s-%s',ssta2,ssta1))
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');
    xlim(xlimrange)
    
    %    saveas(gcf,sprintf('XcorrCompareNondiffsim_%s_%s_C3_conti.pdf',ssta2,ssta1),'pdf')
    
end


if ploteachwf & 1
    figure(9)
    subplot(211)
    %     plot(timedigitsXcorr,xcorr(wftimesta(:,sta1),wftimesta(:,sta2)))
    plot(timedigitsXcorr,xcorrdiffE)
    xlim(xlimrange)
    xlabel('time (sec)')
    title(sprintf('Xcorr.E of fully diffuse noise for %s-%s',ssta2,ssta1))
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');
    subplot(212)
    % ONLY NON=DIFFUSE NOISE
    plot(timedigitsXcorr, xcorrnonE-xcorrdiffE)
%         plot(timedigitsXcorr,xcorr(datasetE.wftimestanon(:,sta1),datasetE.wftimestanon(:,sta2)))
    % XCORR of BOTH DIFF AND NON-DIFF
%     plot(timedigitsXcorr,xcorrnonE)
    
    xlabel('time (sec)')
    title(sprintf('Xcorr with non-diffuse noise for %s-%s',ssta2,ssta1))
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');
    xlim(xlimrange)
    
    %    saveas(gcf,sprintf('XcorrCompareNondiffsim_%s_%s_C3_conti.pdf',ssta2,ssta1),'pdf')
    
end
clear datasetE

