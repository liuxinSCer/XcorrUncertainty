%%% produce the frequency domain equivalent of
%%% the time domain long record (many days)
%%% XIN LIU
%%% STANFORD, 2017

%%

clear
rng('shuffle')

subdir = 'dirwfdatanewnew';

ifNoniso=false; %true; % if true, use anisotropic source distributions

ifcompE=true; %false; % false;

plotmap=1;
initParams
% ddist=20; % delta distance; km
% dist=300:ddist:580;% in  meters
% niter=2000;% # of iterations


%% COORDINATES OF DIFFUSE SOURCES
% theta=0:dtheta:2*pi;
theta=(0:nrv-1)/(nrv-1)*(2*pi); % azimuth angle, clockwise;has to be row vector!
phi=pi/2-theta;% angle between radius and x-axis; counter-clockwise
srccoords=zeros(nrv,2);% source loc; Y, X coordinates in meters
[srccoords(:,2),srccoords(:,1)]=pol2cart(phi',R);
srccoords=bsxfun(@plus, srccoords, origin);
nonisoArr=2-cos(phi);











srcNdispersion

%% FORMAL SIMULATION BEGINS
% % % % % % % wftimeA=zeros(nfreq,1);
% % % % % % % wftimeB=zeros(nfreq,1);
wftimesta=zeros(nptsS,nsta);
timedigitsS=((1:nptsS)-1)*dt;

% wftimesta=zeros(nptsS+nptsG-1,nsta);
%  timedigitsS=((1:(nptsS+nptsG-1) )-1)*dt;

% % DEFINE NOISE SPECTRA STANDARD DEV for diffuse noise
% noiseBsqrt=1;noiseBstdv=0.2;%noiseBstdv=0.05;

% onevector=ones(nrv,1);
rng('default');
rng(1); % seed==1 :: SET UP A REPEATABLE RANDOM NO. GENERATOR

vecEast=[0,1.0];% East direction

% parfor isrc=1:nrv
for isrc=1:nrv    % RND SRC
    if ifNoniso
        noiseB=noiseBstdv*nonisoArr(isrc);
    else
        noiseB=noiseBstdv;
    end
    specSrc=( zeros(nptsS,1) );
    specSrc(freqindS)=noiseB*randn(nfreqS,1)+1i*noiseB*randn(nfreqS,1);
%     specSrc(freqindS)=noiseBstdv*randn(nfreqS,1)+1i*noiseBstdv*randn(nfreqS,1);
    specSrc(1:midpS)=specSrc(1:midpS).*hFIRSrc;
    %specSrc(midpS+1:end)=conj(specSrc(midpS:-1:2));
    Srcoft=ifft(specSrc,'symmetric');
    Srcoft=Srcoft.*taperwinS;% apply taper window
    isrccoord=srccoords(isrc,:);
    tic
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
            noisewf=projcos*conv(Srcoft,Goft,'same');
        else

            noisewf=conv(Srcoft,Goft,'same');
        end
        %         noisewf=conv(Srcoft.*taperwinS,Goft);
        wftimesta(:,ista)=wftimesta(:,ista)+noisewf;
    end
    toc
    %     figure(6);plot(timedigitsS,wftimesta(:,1))
end












%% SAVE ALL PARAMETERS

if ifNoniso
    % BELOW IS FOR ANISOTROPIC NOISE SOURCE DISTRIBUTION
    if ~ifcompE
        fnamesave=sprintf('~/scwork/%s/wfdatatimeSyndiffuseAnisoCA_HHZ_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
    else
        fnamesave=sprintf('~/scwork/%s/wfdatatimeSyndiffuseAnisoCA_HHE_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
    end
    
else
    % % BELOW IS FOR ISOTROPIC NOISE SOURCE DISTRIBUTION
    if ~ifcompE
        fnamesave=sprintf('~/scwork/%s/wfdatatimeSyndiffuseisoCA_HHZ_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
    else
        fnamesave=sprintf('~/scwork/%s/wfdatatimeSyndiffuseisoCA_HHE_%2.1f_days_%d_new.mat',subdir,ndays,nsta);
    end
end
% fnamesave=sprintf('~/scwork/wfdatatimeSyndiffuse_%2.1f_days_%d_new.mat',ndays,nsta);

save(fnamesave,'wftimesta',...
    'timedigitsS','dfS','dt','nptsS','freqwin','freqlower','freqhigher', ...
    'coorddict','srccoords','stalist')

nptswf=length(wftimesta(:,1));
timedigitsXcorr=(-nptswf+1:nptswf-1)*dt;
figure(11);
plot(timedigitsXcorr,xcorr(wftimesta(:,1),wftimesta(:,2)))
return
