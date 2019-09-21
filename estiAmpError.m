function stderrTD=estiAmpError(fltrcoef,freqerrornorm,subindspec,npts)

midp=length(fltrcoef);
npts=midp*2-1;

spechalf=zeros(midp,1);

% FILTER THE ERROR TIME SERIES
spechalf(subindspec)=freqerrornorm;
spechalf=spechalf.*fltrcoef;

vecnum=(0:npts-1)';
% USE COMPLEX CONJUGATE BELOW: original:
varvec= (spechalf'*spechalf)*(1+cos(2*pi*vecnum/npts)*8/(pi*pi))/(npts*npts);

% EMPERICAL UPPER LIMIT:
% varvec= (spechalf'*spechalf)*2*(1+cos(2*pi*vecnum/npts)*8/(pi*pi))/(npts*npts);

stderrTD=sqrt(varvec/2); 
% stderrTD=sqrt(varvec/1);
aa=5;