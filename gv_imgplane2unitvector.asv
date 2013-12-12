function u = gv_imgplane2unitvector(calarray, Rinv, dat2din)

% function to find unit vectors pointing in the direction corresponding to
% 2d image plane coordinates. 
%Similar to gv_pixel2unitvector but accepts a calibration matrix rather
%than a structure, and accepts undistorted image plane coordinates rather than 
% pixels.  These differences are needed to use dyanmic calibration.
%
% inputs:
%    calarray  -- 8 by ncams array with calibration parameters
%   Rinv  -- Inverse rotation matrix--could be calculated from calarray,
%       but this is faster
%   dat2din   -- array(N by 2) containing undistorted image plane
%       coordinates

% turns a 2d position into a unit vector.  The position of a particle in 3D space
% is along the ray given by the point cal.Tinv + lambda*(unit vector) where
% lambda is any real number.

%create arrays 
L=size(dat2din,1);
pos = zeros(L,3);
u=zeros(L,3);

pos(:,1:2)=dat2din;
%include distortion.  k1 seems to be defined using the radial distance
%after distortion
%rather than using the undistorted radius (because calibProj_Tsai iterates) so we
%don't iterate going this way.  I still have some questions about this
%since http://homepages.inf.ed.ac.uk/rbf/CVonline/LOCAL_COPIES/DIAS1/
%puts the (1+k1r_d^2) term as multiplying the distorted position rather than dividing as calibProjTsai
%does.  It is all conventions in the end though--just need to be
%consistent.

radius2 = pos(:,1)^2 + pos(:,2)^2;
pos(:,1:2) = pos(:,1:2)/(1+calarray(8)*radius2);

%Now scale by the effective focal length and the z coordinates

pos(:,1:2)=pos(:,1:2) .* (calarray(6) / calarray(7));
pos(:,3)=calarray(6);

%we also need to map the pinhole which is a common point on the ray for any
%pixel.  This is the same calculation as for pixin=[0 0], but the initial z
%coordinate is zero rather than cal.T(3)
 
zpoint = [0, 0, 0]; 
zpoint = (zpoint - (calarray(4:6))') * (Rinv)';

for i=1:L
    pos(i,:)=(pos(i,:) - (calarray(4:6))') * (Rinv)'; % I found this by reversing the steps in calibProj_Tsai
    u(i,:)=(pos(i,:)-zpoint)/norm(pos(i,:)-zpoint); %here we subtract zpoint and normalize to a unit vector
end





