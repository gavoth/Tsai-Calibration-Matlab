function R = gv_angles2rotmat(angles)
%This function calculates the angles that correspond to a rotation matrix.
%It is needed so that dynamic calibration can be performed on the minimum
%number of parameters rather than on the full rotation matrix.  

%inputs:
%       angles: three angles that define the rotation matrix--theta, phi,
%       and mu where mu is first rotation about z axis, phi is second
%       rotation around y axis and theta is third rotation around z axis.
%       rotout:  best match rotation matrix
%outputs: 
%       R :  a rotation matrix 


%some checks 
%determinant = det(rot) %should be -1 (for some reason this calibration uses reflection and rotation matrix
%check1=(rot(1,3)^2 + rot(2,3)^2+rot(3,3)^2) % should be 1
%check2 = (rot(3,1)^2 + rot(3,2)^2+rot(3,3)^2) %should be 1

theta = angles(1);
phi = angles(2);
mu = angles(3);

rotztheta=[cos(theta), -sin(theta), 0; +sin(theta), cos(theta), 0; 0 0 1];
rotzmu=[cos(mu), -sin(mu), 0; +sin(mu), cos(mu), 0; 0 0 1];
rotyphi=[cos(phi), 0, sin(phi); 0 1 0; -sin(phi), 0, cos(phi)];
reflect = [1 0 0;0 1 0;0 0 -1];

R = reflect*rotztheta*rotyphi*rotzmu;

