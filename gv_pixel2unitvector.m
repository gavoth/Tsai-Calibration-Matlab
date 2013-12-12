function u = gv_pixel2unitvector(cal, pixin)
% function to find unit vectors pointing in the direction corresponding to
%  pixel coordinates. 
%
% inputs:
%    cal  -- calibrated camera parameters (structure with .T, .Tinv, .R etc
%   pixin   -- array(N by 2) containing pixel coordinates

% turn a pixel into a unit vector.  The position of a particle in 3D space
% is along the ray given by the point cal.Tinv + lambda*(unit vector) where
% lambda is any real number.

%create arrays 
L=length(pixin(:,1));
pos = zeros(L,3);
u=zeros(L,3);

%First convert pixel number into position on the image plane

pos(:,1)=(pixin(:,1)-cal.Npixw/2-cal.Noffw)*cal.wpix;
pos(:,2)=(-pixin(:,2)+cal.Npixh/2-cal.Noffh)*cal.hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative

%now include distortion.  k1 seems to be defined using the radial distance
%after distortion rather than using the undistorted radius (because calibProj_Tsai iterates) so we
%don't iterate going this way. (We have distored (pixel) coordinates)  I still have some questions about this
%since http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/DIAS1/
%puts the (1+k1*r_d^2) term as multiplying the distorted position rather
%than dividing as we need to do here to match calibProj_Tsai.
%This is all conventions in the end though--just need to be consistent.

radius2 = pos(:,1)^2 + pos(:,2)^2;
pos(:,1:2) = pos(:,1:2)/(1+cal.k1*radius2);

%now scale horizontal coordinates by the effective focal length and assign
%the vertical coordinate.
pos(:,1:2)=pos(:,1:2) .* (cal.T(3) / cal.f_eff);
pos(:,3)=cal.T(3);

%To get the unit vector, we also need to map the pinhole which is a common point on the ray for any
%pixel.  This is the same calculation as for pixin=[0 0], but the initial z
%coordinate is zero rather than cal.T(3)
 
zpoint = [0, 0, 0]; 
zpoint = (zpoint - (cal.T)') * (cal.Rinv)';

for i=1:L
    pos(i,:)=(pos(i,:) - (cal.T)') * (cal.Rinv)'; % I found this by reversing the steps in calibProj_Tsai
    u(i,:)=(pos(i,:)-zpoint)/norm(pos(i,:)-zpoint); %here we subtract zpoint and normalize to a unit vector
end





