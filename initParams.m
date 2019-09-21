% % % % % INITIALIZE ALL FOURIER TRANFORM PARAMETERS
% % % % % AND STATIONS LOCATIONS
% % % % % XIN LIU, STANFORD UNIV, APRIL 2017


% IF TRUE, USE CA STATIONS; OTHERWISE: SYNTHETIC LOCATIONS
useCAsta=false; %true;%true; %false;%1;

R=600*1e3; % radius of sources circle; in meters
nrv=600; % No. of diffuse noise sources
%% DEFINE SAMPLING FREQUENCY
Fsamp=2; dt=1/Fsamp;
Fny=Fsamp/2;

if ~useCAsta
    %% DEFINE STATION LIST AND COORDINATES
    stanames='syn1 syn2 syn3 syn4';
    stalist=strsplit(stanames);
    coorddict=containers.Map();
    nsta=length(stalist);
    % DEFINE CARTETIAN COORDINATE SYSTEM
    origin=[0,20]*1e3;% Y, X coordinates in meters
    % 4 LINES BELOW ARE FOR PREVIOUS SYNTHETIC GEOMETRY...
%     coorddict('syn1')=[0,-20]*1e3+origin;% Y, X coordinates in meters
%     coorddict('syn2')=[0,20]*1e3+origin;% Y, X coordinates in meters
%     coorddict('syn3')=[-20,0]*1e3+origin;% Y, X coordinates in meters
%     coorddict('syn4')=[20,0]*1e3+origin;% Y, X coordinates in meters

%     % 4 LINES BELOW ARE FOR NEW SYNTHETIC GEOMETRY: 60 km distance ...
    coorddict('syn1')=[0,-30]*1e3+origin;% Y, X coordinates in meters
    coorddict('syn2')=[0,30]*1e3+origin;% Y, X coordinates in meters
    coorddict('syn3')=[-30,0]*1e3+origin;% Y, X coordinates in meters
    coorddict('syn4')=[30,0]*1e3+origin;% Y, X coordinates in meters
else
    %% DEFINE CI station coordinates
    % (CONVERT UNIT TO METERS)
    stanames='CHF LMR2 SBB2 IPT ADO BFS'; % 6 stations
% % % % % % %     stanames='CHF LMR2 SBB2 BFS IPT ADO EDW2 TA2 LUG MWC CJM LUC2 RRX BBR SVD CLT FHO';
%     stanames='CHF SBB2 LMR2 ADO BFS IPT EDW2 TA2 LUG MWC CJM LUC2 RRX BBR SVD CLT FHO';
    stalist=strsplit(stanames);
    initMapProj % SET UP MAP PROJECTION AND STATION COORDINATES
end
%% DEFINE LENGTH OF NOISE REC
ndays=15;% number of days
% ndays=0.10;% number of days
lentracetotal= 60*60*24*ndays; % in sec
% df=0.005;
df=1.0/lentracetotal;
freqwin=[0.05 1.0];% DEFINE THE FREQUENCY WINDOW
% freqv=0.05:df:0.8; % in hertz
freqv=freqwin(1):df:freqwin(2); % in hertz
freqv=freqv(:);
nfreq=length(freqv);

% omega=2*pi*freqv;