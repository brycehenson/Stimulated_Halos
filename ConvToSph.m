function [halo_radial,halo_azm,halo_elev]=ConvToSph(single_halo)
%======================================90char=============================================
%+++++++++++++++++++++++++++++++++++++++ABOUT+++++++++++++++++++++++++++++++++++++++++++++
%this function converts centered TXY data into radial cord

% ref system
% bec's are at the poles
% azimuthal angle is arround the equator (zero to 2 pi)
%
%          y pi
%            |
%            |
% pi/2 --------------x 3 pi/4
%            |
%            |
%            0,2pi

% inclination is angle from the poles (which is zero inclination)(zero to pi)
% this was chosen to reduce the wrapping problems arround zero/2 pi 
% 

halo_radial=sqrt(single_halo(:,1).^2+single_halo(:,2).^2+single_halo(:,3).^2);
halo_azm=atan2(single_halo(:,2),single_halo(:,3))+pi;
halo_elev=acos(single_halo(:,1)./halo_radial);
end
