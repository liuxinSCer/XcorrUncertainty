clear

diffuseOnly = false;true; false;

nsrcnon=3;

fsrcname='srcNondiff_continuous.mat';
fsrcname='srcNondiff_1src.mat';
% % % % % fsrcname='srcNondiff_5srcs.mat';

fsrcname=sprintf('srcNondiff_%dsrc.mat',nsrcnon);


load(fsrcname,'R','phinon','srccoordsnon','nsrcnon')

simtypetag='SYN';%'CA';% 'CA' or 'SIM'
nrv=1000;
%% COORDINATES OF DIFFUSE SOURCES
% theta=0:dtheta:2*pi;
% theta=(0:nrv-1)/(nrv-1)*(2*pi); % azimuth angle, clockwise;has to be row vector!
theta=rand(nrv,1)*2*pi;
phi=pi/2-theta;% angle between radius and x-axis; counter-clockwise
ringwidth=50;% km
radial=(R*1e-3-ringwidth/2)+ringwidth*rand(nrv,1);
% srccoords=zeros(nrv,2);% source loc; Y, X coordinates in meters
% [srccoords(:,2),srccoords(:,1)]=pol2cart(phi',R);
% srccoords=bsxfun(@plus, srccoords, origin);
figure(10)
polarplot(phi,radial,'b.')
if ~diffuseOnly
    hold on
    polarplot(phinon,R*1e-3+ringwidth,'ro','MarkerFaceColor','r')
    hold off
end
pax=gca;
angles = 0:45:360;
pax.RTickLabel=[];
pax.ThetaTick = angles;
labels = {'E','NE','N','NW','W','SW','S','SE'};
pax.ThetaTickLabel = labels;
if ~diffuseOnly
legend('fully diffuse noise source','non-diffuse (physical) source')
else
legend('diffuse noise source')
end
fontsize=18;
legend boxoff
    set(gca,'FontSize',fontsize)
    set(gcf,'PaperPositionMode','auto');
        saveas(gcf,sprintf('Source_polarplot_%s_%dsrcnondiff.pdf',simtypetag,nsrcnon),'pdf')

