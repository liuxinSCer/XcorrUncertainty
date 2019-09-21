function onecorrdat=stackRandSampleithcorr(rndint,bb,ixcorr)
%% REUTRNS ONLY ONE COLUMN VECTOR AS THE STACKED SPECTRA FOR ixcorr component
ndays=length(rndint);
nfreqs=size(bb{1},2);

onecorrdat=complex(zeros(nfreqs,1),0);
% for ixcorr=1:nxcorr
    onecorrdat=mean(bb{ixcorr}(rndint,:))';
    
    
    
% end