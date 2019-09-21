function bootfreqarray=bootstrapblock(freqarrayselect,ndayactual,freqrangedata,Nbootstrap,getStackedWins,keypair)
% Nbootstrap=4000;% NUMBER OF BOOTSTRAP 
% getStackedWins: if true, return stackedday matrix
ploteachfig=1;0;
DF=freqrangedata(2)-freqrangedata(1);
nfreqstore=length(freqrangedata);
AKInorm=1;
fontsize=16;
% dayrangeinterval=dayrange(1):dayrange(2);
% ndaysselect=length(dayrangeinterval);
% ndayactual=ndaysselect;% this is the actual number of days!
% BELOW: 24*15 means block size is 1/15 of an hour: that is 4 minutes!
ndaysselect=24*15*ndayactual;%6000;% temporarily make this bin size much smaller!
% ndaysselect=24*30*ndayactual;%6000;% temporarily make this bin size much smaller!
% ndaysselect=24*5*ndayactual;%6000;% temporarily make this bin size much smaller!
% ndaysselect=24*2*ndayactual;%6000;% temporarily make this bin size much smaller!
kbinsperday=ndaysselect/ndayactual;% how many windows per day
bb={};% stored daily xspecs;
%%% NOW STACK NEIGHBORING COLUMNS INTO DAYS
ixcorr=1;nxcorr=1;% ONLY ONE CROSS-CORRELATION
nrec=size(freqarrayselect,1);
binsize= floor( nrec/ndaysselect);
%binsize=1;% FOR VERY SMALL bin SIZE
stackedday=complex(zeros(ndaysselect,nfreqstore),0);
for iday=1:ndaysselect
    rowinds=1+binsize*(iday-1);
    rowinde=min([binsize*iday, nrec]);
    row= mean( freqarrayselect(rowinds:rowinde,:),1 );
    stackedday(iday,:)=row;
    
end
bb{ixcorr}=stackedday;





rng('shuffle') % random number generator
meanboots=complex(zeros(nfreqstore,nxcorr));
stdbootsreal=zeros(nfreqstore,nxcorr);
stdbootsimag=zeros(nfreqstore,nxcorr);
% stores all BOOTSTRAP RESULTS
% bootfreqarray=complex(zeros(nfreqstore,Nbootstrap));

freqwin=[freqrangedata(1) freqrangedata(end)];
% freqmeasure=0.07:DF*2:0.59;
freqmeasure=0.10:DF:0.59;
freqmeasure=freqrangedata(1):DF:freqrangedata(2);
inddistrifreq=round((freqmeasure-freqrangedata(1) )/DF+1);
nstatday=5;% LOOP OVER HOW MANY DAYS! can be a fraction!
nseg=ceil( ndayactual/nstatday);
% nseg=floor( ndayactual/nstatday);

snrpFarray=[];
% for ixcorr=1:nxcorr
 %for iseg=1:nseg % FOR CUMULATIVE DAYS
  for iseg=nseg:nseg  % FOR THE MAX NUMBER OF DAYS ONLY..
    ixcorr=1; % for first crosscorrelation pair; test purposes only
    iday=iseg*nstatday;
    if iday>ndayactual
       iday=ndayactual; 
    end
    ibinsselect=iday*kbinsperday;
    % init the matrix that stores all bootstrap samples
    bootsarr=complex(zeros(nfreqstore,Nbootstrap),0);
    iday
%     parfor iboot=1:Nbootstrap
    for iboot=1:Nbootstrap
        display iboot
        % GENERATE A RANDOM SAMPLE OF GIVEN NUMBER OF DAYS
        %%IMPORTANT%
        %dayrange=6*ones(ndays,1); % I'VE TESTED UP TO DAY 75 FOR JF ARRAY
        rndint=randi(ibinsselect,[ibinsselect 1]);
        % IF I WANT JUST SIMPLE STACK
        % rndint=1:ndaysselect;
        % BELOW GET THE XSPECTRA OF ALL THREE CORR PAIRS
        onedatspec=stackRandSampleithcorr(rndint,bb,ixcorr); % bb{ixcorr} is ndays*nfreqs matrix
        % NOW ADD CURRENT STACKED BOOTSTRAP SAMPLE SPECTRA TO THE MATRIX
        bootsarr(:,iboot)=onedatspec;
        
        
        % ANALYZE THE ATTENUATION ON TRIPLET
        %         Qvalue=getattenuationLowFreq(stalist,tripletdatspec,autocorrfreq, ...
        %             windowlength,correctedphaseSlow(:,1:nxcorr),freqrange,subfreqwin,Fsamp,distarr,freqrangedata([1 end]),measureINtime,figdir)
        %
        %
        %         Qinvstore(iboot,:)=1./Qvalue;
    end
    meanboots(:,ixcorr)=mean(bootsarr,2);% MEAN OF Nbootstrap samples
    stdbootsreal(:,ixcorr)=std(real(bootsarr),0,2);
    stdbootsimag(:,ixcorr)=std(imag(bootsarr),0,2);
%     snrpF=freqcausalSNR(conj(meanboots(:,ixcorr)),stdbootsreal(:,ixcorr),stdbootsimag(:,ixcorr), ...
%         freqrangedata,inddistrifreq,AKInorm,freqwin);
%     snrpFarray=[snrpFarray, snrpF];
    if(iseg==nseg)
        bootfreqarray=bootsarr';
    end
  end



if getStackedWins
    
   bootfreqarray=stackedday; 
end
%% PLOT ALL SNR RESULTS
niter=size(freqarrayselect,1);
colordef='jet';
iterarray=(1:nseg)*nstatday;
% figure(7);imagesc(freqmeasure,iterarray,snrpTDarray')
% % colormap parula
% colormap jet
% colorbar
% set(gca,'YDir','normal')
% % set(gca,'YTickLabel',iterarray)
% xlabel('frequency (Hz)')
% ylabel('number of days')
% title('SNR from time domain narrowband cross-correlation')
% set(gca,'Fontsize',fontsize)
% set(gcf,'PaperPositionMode','auto');
% saveas(gcf,sprintf('SNR_timedomain_%s-%s__%d_iterations%s.eps',sta1,sta2,niter,colordef),'psc2')
% saveas(gcf,sprintf('SNR_timedomain_%s-%s__%d_iterations%s.pdf',sta1,sta2,niter,colordef),'pdf')


% snrpFarray(:,1)=0;
if ploteachfig
    figure(8);imagesc(freqmeasure,iterarray,(snrpFarray'))
    % colormap parula
    colormap jet
    colorbar
    set(gca,'YDir','normal')
    % ylim([0 10000])
    xlabel('frequency (Hz)')
    ylabel('number of days')
    title({'SNR from frequency domain';'mean/std cross-spectrum'})
    set(gca,'Fontsize',fontsize-1)
    set(gcf,'PaperPositionMode','auto');
%     saveas(gcf,sprintf('SNR_bootnew_freqdomain_%s__%d_iterations%s.eps',keypair,niter,colordef),'psc2')
    saveas(gcf,sprintf('SNR_bootnew_freqdomain_%s__%d_iterations%scausal.pdf',keypair,niter,colordef),'pdf')
end
 %save(sprintf('dataSNR_bootstrapnew_%s__%d_iterationsNewcausal.mat',keypair,niter),'snrpFarray','bootfreqarray','keypair','freqmeasure','iterarray','Nbootstrap','nseg')

aa=1;
