function  [newDenMat,PSDcenters]=assemDensityImageXcorr(xcorrarrayselected)
% call the hist function and compute hist bins for each
% column in indexfreqarrayselected
maxval=max(max(xcorrarrayselected));
minval=min(min(xcorrarrayselected));
nbins=1000;nbins=500;nbins=200;
% nbins=1000; % (for coda only ...)
xbinlog=minval:(maxval-minval)/nbins:maxval;
%hist(indexfreqarrayselected);
 [newDenMat,PSDcenters]=hist(xcorrarrayselected,xbinlog);

% newDenMatlog=[]