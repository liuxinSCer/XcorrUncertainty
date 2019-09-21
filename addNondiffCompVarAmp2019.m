%%%%% ADD NON-DIFFUSE COMPONENTS TO FULLY DIFFUSE NOISE FIELDS
%%%%% Xin Liu, Stanford Univ, 2017
%%%%% ALRIGHTS RESERVED
clear

subdir = 'dirwfdatanew';
% subdir = 'dirwfdata';
ifDebug = false;%  true;
ifLoveNondiff= false ; %false;
ifcompE= false;false;true; false; false;true; % false; %
ifcompN= false;false; % North == Transverse component for station pair Syn1-Syn2
if ifcompN
    ifcompE=false;
end

ploteachwf=0;% SWITCH FOR MAKING PLOTS

initParams
srcNdispersion

% CA STATIONS AND 1D DISP
% fname=sprintf('~/scwork/dirwfdata/wfdatatimeCAdiffuse_%2.1f_days_%d_new.mat',ndays,nsta);

% SYNTHETIC STATION LOCATION & DISPERSION
% fname=sprintf('~/scwork/wfdatatime_%2.1f_days_%d.mat',ndays,nsta);
% if ifcompE || ifcompN
%     fname=sprintf('~/scwork/%s/wfdatatimeSyndiffuse_HHE_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
% else
%     fname=sprintf('~/scwork/%s/wfdatatimeSyndiffuse_HHZ_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
% %     fname='/Users/xinliu/scwork/wfdatatimeCAdiffuse_15.0_days_6_new.mat';
% end

% NEW DATASET AND DIRECTORY:
if ifcompE || ifcompN
    fname=sprintf('~/scwork/%s/wfdatatimeSyndiffuseisoCA_HHE_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
else
    fname=sprintf('~/scwork/%s/wfdatatimeSyndiffuseisoCA_HHZ_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
%     fname='/Users/xinliu/scwork/wfdatatimeCAdiffuse_15.0_days_6_new.mat';
end
% % % % % fname=sprintf('C:/Users/liusi/scwork/wfdatatime_%2.1f_days_%d.mat',ndays,nsta);
if  ~ifDebug
    load(fname)
end

if ploteachwf
    nptswf=length(wftimesta(:,1));
    timedigitsXcorr=(-nptswf+1:nptswf-1)*dt;
    figure(8);
    plot(timedigitsXcorr,xcorr(wftimesta(:,2),wftimesta(:,3)))
end


%% NON-DIFFUSE NOISE PART

% COORDINATES OF Non-DIFFUSE SOURCES
% theta=0:dtheta:2*pi;
dthetanon=pi/180;% SPACING IS 1-degree?
% dthetanon=16*pi/180;% SPACING IS 16-degree?

thetanon=pi:dthetanon:1.5*pi; % azimuth angle, clockwise;has to be row vector!
thetanon=(30+180)/180*pi; % azimuth angle, clockwise;has to be row vector!
degsrc=210:10:240;
degsrc=210:20:230;
degsrc=170:20:210;
degsrc= 220; % 225; %270;%220;%250; # GOOD FOR PAPER!
% degsrc=270;%250; # GOOD
degsrc=250;%250; # GOOD
% NEW TEST MARCH 2019:
degsrc=200:15:230;
ampsrc=[0.2,1,0.2];
% NEW TEST 2:
degsrc=[220,130];
ampsrc=[1, 0.7]*1.5;
% NEW TEST 3:
degsrc=[130,220];
ampsrc=[1.0, 1.0]*1.5;
% NEW TEST 4:
degsrc=[30,130,220];
ampsrc=[.5,.5, 1.0]*1.5;
%
thetanon=degsrc/180*pi;
% thetanon=(270)/180*pi; % azimuth angle, clockwise;has to be row vector!
phinon=pi/2-thetanon;% angle between radius and x-axis; counter-clockwise
nsrcnon=length(thetanon);
srccoordsnon=zeros(nsrcnon,2);% source loc; Y, X coordinates in meters
[srccoordsnon(:,2),srccoordsnon(:,1)]=pol2cart(phinon',R);
% % % % % % % fsrcname='srcNondiff_5srcs.mat';
% fsrcname='srcNondiff_cont.mat';
fsrcname=sprintf('srcNondiff_%dsrc.mat',nsrcnon);

save(fsrcname,'R','phinon','srccoordsnon','nsrcnon')


% polargeo(phinon',R*ones(nsrcnon,1)*1e-3,'ro')

srccoordsnon=bsxfun(@plus, srccoordsnon, origin);% IF origin IS NOT ZERO

% LOAD Non-diffuse SOURCE MAT
Xfreqdata=load('../decomposenoise/Xfreqcomps_CHF.mat');

subfreq=[0.05, 0.59];
subfreqold=[0.1, 0.59];% high freq[0.1, 0.59];% high freq
% svecmat=Xfreqdata.svecmatnew;
nthNondiffcomp=4;
% nthNondiffcomp=1:10;
svecmat=Xfreqdata.svecmatnew(:,nthNondiffcomp);
[nfreqnon,nvec]=size(svecmat);
% FREQUENCY SPACING: DETERMINS WINDOW LENGTH ()
% dfXfreq=0.01;% in Hz :: 100 sec window
dfXfreq=diff(subfreqold)/(nfreqnon-1); %

% ZERO PADDING:
nzero = round( (subfreqold(1) - subfreq(1))/ dfXfreq );
svecmat=[zeros(nzero,1);svecmat];
nfreqnon = nfreqnon + nzero;
%

timeX=round(1/dfXfreq);% in sec
nptsX=timeX/dt+1;
midpX=(nptsX+1)/2;
% DEFINE FREQ VALUES FOR Xfreq components
freqvXfreq=subfreq(1):dfXfreq:subfreq(2);
freqindXfreq=round(freqvXfreq/dfXfreq)+1;
nfreqXfreq=length(freqvXfreq);
%
xgaussparam=fitGaussian(svecmat, dfXfreq, freqvXfreq);
newGaussian=@(x0) x0(3)*exp(-(freqvXfreq-x0(1)).^2/(2*x0(2)*x0(2)));
svecmat = newGaussian(xgaussparam);
svecmat=svecmat(:);
%
% DEFINE FILTER TRANSITION BAND WIDTH AND FIR FILTER POINTS NO.
FbandwidthX=0.02;bdwdX=FbandwidthX; % in Hz; minimum bandwidth between pass to cut-off frequency.
BW=FbandwidthX/Fsamp; % normalized bandwidth by sampling frequency Fsamp.
MfilterX=4/BW; % theoretical number of points for filter kernal to satisfy the bandwidth requirement.
% filterpts=50;
filterptsXfreq=MfilterX; %filterpts=100; % number of order for the FIR filter
fwinX=(subfreq+[bdwdX+0.0 -bdwdX])/Fny;
b = fir1(filterptsXfreq,fwinX);
[hFIR,w]=freqz(b,1,midpX);
wfreqX=w*Fny/pi;
hFIR2Xfreq=abs(hFIR);% MAKE IT A ZERO-PHASE ACAUSAL FILTER
% hFIR2Xfreq=(hFIR);% MAKE IT A CAUSAL FILTER

%%% PLOT FIR FILTER RESPOSNE BELOW
figure(2);semilogy(wfreqX,hFIR2Xfreq)

XfreqTimeArr=zeros(nptsX,nvec);

timedigitsXfreq=(-midpX+1:midpX-1)*dt;

for ivec=1:nvec
    specX=complex( zeros(nptsX,1) );
    
    specX(freqindXfreq)=svecmat(:,ivec);
    specX(1:midpX)=specX(1:midpX).*hFIR2Xfreq;
    %specG(midpG+1:end)=conj(specG(midpG:-1:2));
    
    XfreqTime=ifft(specX,'symmetric');
    
    XfreqTimeCentre=[XfreqTime(midpX+1:end);XfreqTime(1:midpX)];
    XfreqTimeCentre=XfreqTimeCentre.*tukeywin(nptsX,0.10);% MULTIPLY TAPER WINDOW
    
    XfreqTimeArr(:,ivec)=XfreqTimeCentre;
    
    figure(3);plot(timedigitsXfreq,XfreqTimeCentre)
    
end

if ploteachwf
    
   figure(12);imagesc(1:nvec-1,timedigitsXfreq,XfreqTimeArr(:,1:end-1)) 
   figure(12);imagesc(1:nvec,timedigitsXfreq,XfreqTimeArr(:,1:end)) 
   ylabel('time (s)')
   xlabel('Cross-frequency Vector No.')
   colorbar
end

%% MAKE PLOTS FOR THE CROSS-FREQ COMPONENT:

corrmatSyn=assemCorrmatFromSvec(svecmat,nfreqnon);

fontsize=15;

figure(1)


imagesc(freqvXfreq,freqvXfreq,corrmatSyn)
axis square
colormap jet
colorbar
xlabel('frequency (Hz)')
ylabel('frequency (Hz)')
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');


figure(5)
plot(freqvXfreq,svecmat(:,1))
set(gca,'xticklabel',[])
xlabel('frequency (Hz)')
legend(' (P_{nondiff}/P_{diff})^{0.5}')
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');



%% COMPUTE POWER RATIO OF NON-DIFFUSE/DIFFUSE NOISE FOR EACH COMPONENT

Pratiofreq=sum(svecmat.^2)/nfreqnon;% 2-NORM SQUARE FOR EACH COLUMN (Xfreq component)

% CALCULATE FULLY DIFFUSE NOISE POWER
tapercornernptshalf=round(0.5*tapercospercent*nptsS);

Pdiff=nrv*4*0.2^2*sum(hFIRSrc.^2)*0.95/nptsS;
Pdiff=Pdiff/nsrcnon;% CALCULATE THE AVERAGE REFERENCE POWER FOR EACH NON-DIFFUSE SRC

% Pdiff=norm(wftimesta(tapercornernptshalf:end-tapercornernptshalf,:)).^2;
% % % % % Pdiff=zeros(nsta,1);
% % % % % for ista=1:nsta
% % % % %     Pdiff(ista)=norm(wftimesta(tapercornernptshalf:end-tapercornernptshalf,ista),2)^2;% L-2 NORM OF EACH COLUMN (one-station SIGNAL)
% % % % % end

%% COMPUTE POWER SPECTRA FOR THE PETERSON (1993) MODEL

% noiseBsqrt=1;noiseBstdv=0.2;%noiseBstdv=0.05;
veloc=zeros(nfreqXfreq,1); % for velocity
NOISE_MODEL='NLNM'; % PETTERSON NEW LOW NOISE MODEL

for l=1:nfreqXfreq
    [~ ,veloc(l), ~]=PetersonNoiseModel(1/freqvXfreq(l),NOISE_MODEL);
    %accel(l)=1; veloc(l)=1; displ(l)=1;
    
    veloc(l)=(10e18)*10.^(veloc(l)/10);% convert to (m/sec)^2
end


noisePetersonStdv=sqrt(0.5*veloc);
%% COMPUTE rescaling freq coefficients

C=1./sqrt(diag(svecmat*svecmat'+eye(nfreqnon)));
    specX=complex( zeros(nptsX,1) );
    
    specX(freqindXfreq)=C.*noisePetersonStdv;c
    specX(1:midpX)=specX(1:midpX).*hFIR2Xfreq;
    %specG(midpG+1:end)=conj(specG(midpG:-1:2));
    
    XfreqTime=ifft(specX,'symmetric');
    
    XfreqTimeCentreScaling=[XfreqTime(midpX+1:end);XfreqTime(1:midpX)];
    XfreqTimeCentreScaling=XfreqTimeCentreScaling.*tukeywin(nptsX,0.10);% MULTIPLY TAPER WINDOW
    
% define the anti-Peterson (whitening) wavelet
    specXwhiten=complex( zeros(nptsX,1) );
    
    specXwhiten(freqindXfreq)=1./noisePetersonStdv.^2;
    specXwhiten(1:midpX)=specXwhiten(1:midpX).*hFIR2Xfreq;
    %specG(midpG+1:end)=conj(specG(midpG:-1:2));
    
    XfreqTimewhiten=ifft(specXwhiten,'symmetric');
    
    XfreqTimewhitencenter=[XfreqTimewhiten(midpX+1:end);XfreqTimewhiten(1:midpX)];
    XfreqTimewhitencenter=XfreqTimewhitencenter.*tukeywin(nptsX,0.10);% MULTIPLY TAPER WINDOW


    figure(3);plot(timedigitsXfreq,XfreqTimeCentreScaling)

    figure(4);plot(timedigitsXfreq,XfreqTimewhitencenter)

%% GENERATE NON-DIFFUSE NOISE SOURCES & MAKE THE NOISE WF
rng('default');
rng(1); % seed==1 :: SET UP A REPEATABLE RANDOM NO. GENERATOR

% INIT ARRAY TO STORE NON-DIFFUSE NOISE TIME SERIES
wftimestanon=zeros(nptsS,nsta);

vecEast=[0,1.0];% East direction

            
% parfor isrc=1:nrv
for isrc=1:nsrcnon    % RND SRC
    
%     specSrc=( zeros(nptsS,1) );
    %specSrc(freqindS)=noiseBstdv*randn(nfreqS,1)+1i*noiseBstdv*randn(nfreqS,1);
    % THE LINE BELOW FILTERS THE NOISE SPECTRA
%     specSrc(1:midpS)=specSrc(1:midpS).*hFIRSrc;
    %specSrc(midpS+1:end)=conj(specSrc(midpS:-1:2));
%     Srcoft=ifft(specSrc,'symmetric');
    
%     Srcoft=Srcoft.*taperwinS;% apply taper window
    %Srcoft= zeros(nptsS,1);
    % DEFINE AVERAGE TIME SPACING BETWEEN TWO repeating XFREQ COMPONENTS
    avgTimeSpacingXcomp=100;% in sec
    Nrepeatnon=lentracetotal/avgTimeSpacingXcomp;
    %%%%% TODO: CREATE WAVEFORM WITH RANDOMLY DISTRIBUTED X-FREQ NOISE COMPONENTS
    Srcoft =nondiffSrcGenerator(XfreqTimeArr,Pratiofreq,dt,nptsS,Nrepeatnon,Pdiff);
    Srcoft=Srcoft.*taperwinS;% apply taper window
    % MAKE NONDIFFUSE NOISE SMALLER:
    coefNondiff=1;%0.5
    coefNondiff=ampsrc(isrc);
    Srcoft=Srcoft*coefNondiff;
    % MODIFY THE SPECTRAL AMPLITUDE (WHITEN, THEN APPLY Peterson Low Noise Model)
    Srcoft=conv(Srcoft,XfreqTimeCentreScaling,'same');
    figure(5);plot(Srcoft)
    
    isrccoord=srccoordsnon(isrc,:);
    
    tic
    % LOOP OVER ALL STATIONS AND ADD GENERATE NON-DIFFUSE NOISE:
    for ista=1:nsta
    %parfor ista=1:nsta
        sta=stalist{ista};
        % get G()
        icoord=coorddict(sta);% LOOK UP sta's coordinates
        
        
        specG=complex( zeros(nptsG,1) );
        dist=norm(icoord-isrccoord);% in meters
        specG(freqindDisp)=exp(-1i*2*pi*dist*freqvDisp.*invc);
        specG(1:midpG)=specG(1:midpG).*hFIR2;
        %specG(midpG+1:end)=conj(specG(midpG:-1:2));
        
        Goft=ifft(specG,'symmetric');
        %         figure(2);plot(timedigitsG,Goft)
        if ifcompE
            % PROJECT TO THE EAST COMPONENT:
            srcsta=icoord-isrccoord;% vector from isrc to ista
            
            projcos=dot( vecEast,srcsta)/norm(srcsta);
            noisewfnon=projcos*conv(Srcoft,Goft,'same');
        elseif ifcompN
            % PROJECT TO THE EAST COMPONENT:
            srcsta=icoord-isrccoord;% vector from isrc to ista
            srcsta90degshift=[srcsta(2), -srcsta(1)];
            projcos=dot( vecEast,srcsta90degshift)/norm(srcsta);
            noisewfnon=projcos*conv(Srcoft,Goft,'same');            
        else
            if ifLoveNondiff
                 noisewfnon=0.0;               
            else
                noisewfnon=conv(Srcoft,Goft,'same');
            end
        end
        %         noisewf=conv(Srcoft.*taperwinS,Goft);
        wftimestanon(:,ista)=wftimestanon(:,ista)+noisewfnon;
    end
    toc

    
end

% fnamesave=sprintf('~/scwork/wfdatatimeCA_5srcNondiff_%2.1f_days_%d.mat',ndays,nsta);
if ~ifcompE
    if ifcompN
        fnamesave=sprintf('~/scwork/%s/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHTprojR_%2.1f_days_%d.mat',subdir,degsrc(1),nsrcnon,ndays,nsta);
    else
        if ifLoveNondiff
            fnamesave=sprintf('~/scwork/%s/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHTprojZ_%2.1f_days_%d.mat',subdir,degsrc(1),nsrcnon,ndays,nsta);
        else
            fnamesave=sprintf('~/scwork/%s/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHZ_%2.1f_days_%d.mat',subdir,degsrc(1),nsrcnon,ndays,nsta);
        end
        
    end
else
    fnamesave=sprintf('~/scwork/%s/wfdatatimeSyn_%ddeg_%dsrcNondiff_HHE_%2.1f_days_%d.mat',subdir,degsrc(1),nsrcnon,ndays,nsta);
end
save(fnamesave,'wftimesta','wftimestanon',...
    'timedigitsS','dfS','dt','nptsS','freqwin','freqlower','freqhigher', ...
    'coorddict','srccoords','XfreqTimeCentreScaling','XfreqTimewhitencenter','timedigitsXfreq','timedigitsXfreq','stalist')

