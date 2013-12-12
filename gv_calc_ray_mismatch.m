function [ray3mismatch,h,points3D] = gv_calc_ray_mismatch(calarray, cam2d)

%calarray will be a 2D matrix nparameters by ncams.  Parameters need to be
%ordered--so I'll choose
%1 theta
%2 phi
%3 mu
%4 Tx
%5 Ty
%6 Tz
%7 f_eff
%8 k1
ncams=3;

for icam = 1:ncams
    R(:,:,icam) = gv_angles2rotmat(calarray(1:3,icam));
    Rinv(:,:,icam) = inv(R(:,:,icam));
    Tinv(:,icam) = Rinv(:,:,icam) * (-1* calarray(4:6,icam));
end

npoints=size(cam2d,1);
ray3mismatch=zeros(npoints,1); %average deviation over all cameras
h = zeros(npoints, ncams); %deviation from each camera
point3D= zeros(npoints,3);

for np=1:npoints 
M = zeros(3,3);
pM = zeros(3,ncams);
u = zeros(3, ncams);
for icam = 1:ncams
        % then find the unit vector designated by the camera position and
        % the 2D pixel coordinates. 
    u(:,icam) = gv_imgplane2unitvector(calarray(:,icam), Rinv(:,:,icam), cam2d(np,:,icam));
    uM = eye(3) - u(:,icam) * (u(:,icam))';
    pM(:,icam) = uM * Tinv(:,icam);
    M = M + uM;
end
% find the point minimizing the distance from all rays
p = M \ sum(pM,2);  % sums pm x together for all three cameras.  Makes a column vector, then does inv(M)*sum(pM,2)
%find the distances from each ray.
for icam = 1:ncams
    temp = p - ((p') * u(:,icam)) * u(:,icam) - pM(:,icam);
    h(np,icam) = sqrt(temp' *temp);
end
points3D(np,:)=p;
ray3mismatch(np) = sqrt(mean(h(np,:).*h(np,:)));
end





