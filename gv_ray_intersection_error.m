
ncams=3;
calibimgsize=[256 256];
filepath='C:\Documents and Settings\gvoth\My Documents\data\cornell\calib\jan16\calibpoints\';

% known parameter of the Phantom cameras
camParaknown.Npixh = 256;
camParaknown.Npixw = 256;
camParaknown.hpix = 0.022;	% pixel size (mm)
camParaknown.wpix = 0.022;	% pixel size (mm)

for icam = 1:ncams
	x = myload(strcat(filepath,'calibpointspos.cam', num2str(icam-1), '.dat'), '#');
	[camParaCalib(icam) err Xout] = calib_Tsai(x(:,1:2), x(:,3:5), camParaknown, calibimgsize);
end

realpoint=[0 0 0; 0 0 0; 0.1 0 -0.1];

% calculate the distance between camera axes
M = zeros(3,3);
pM = zeros(3,ncams);
u = zeros(3, ncams);
for icam = 1:ncams
    % a unit vector pointing along camera axes
    cam2d=calibProj_Tsai(camParaCalib(icam), (realpoint(:,icam))');
    u(:,icam) = gv_pixel2unitvector(camParaCalib(icam), cam2d);
    uM = eye(3) - u(:,icam) * (u(:,icam))';
    pM(:,icam) = uM * camParaCalib(icam).Tinv;
    M = M + uM;
end
% point of axes crossing
p = M \ sum(pM,2);  % sums pm x together for all three cameras.  Makes a column vector
h = zeros(1, ncams);
for icam = 1:ncams
    temp = p - ((p') * u(:,icam)) * u(:,icam) - pM(:,icam);
    h(icam) = sqrt(temp' *temp);
end
hbar = sqrt(mean(h.*h));
p
h
hbar