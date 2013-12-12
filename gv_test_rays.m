%calibration unit vector test

%first find the calibration
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

cal = camParaCalib(1);

zin = [0 0 0];
pixin = [0 0 cal.f_eff];

pos = pixin .* cal.T(3)/cal.f_eff;
pos = pos- (cal.T)';
pixout = pos * (cal.Rinv)';
pixout

pos = zin .* cal.T(3)/cal.f_eff;
pos = pos- (cal.T)';
zout = pos * (cal.Rinv)';
zout

uvect=(pixout - zout)/norm(pixout-zout);
u = cal.Rinv * [0 0 1]';
u = u/norm(u);

uvect'
u




