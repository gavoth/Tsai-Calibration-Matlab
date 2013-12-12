
clear all;
ncams=3;
calibimgsize=[256 256];


filepath='C:\Documents and Settings\gvoth\My Documents\data\cornell\calib\jan16\calibpoints\';
camParaknown.Npixh = 256;
camParaknown.Npixw = 256;
camParaknown.hpix = 0.022;	% pixel size (mm)
camParaknown.wpix = 0.022;	% pixel size (mm)
for icam = 1:ncams
	x = myload(strcat(filepath,'calibpointspos.cam', num2str(icam-1), '.dat'), '#');
	[camParaCalib(icam) err Xout] = calib_Tsai(x(:,1:2), x(:,3:5), camParaknown, calibimgsize);
end
 
%for icam=1:ncams
%    camParaCalib(icam).k1=0;
%    camParaCalib(icam).k1star=0;
%end
%camParaCalib(2) = gv_makecalib(0, -pi/2, -pi/2, 0, 0, 0);
%camParaCalib(1) = gv_makecalib(+pi/4, -pi/4, -pi/2, 0, 0, 0);
%camParaCalib(3) = gv_makecalib(-pi/4, -pi/2, -pi/2, 0, 0, 0);


npoints=3^3;
point3d=zeros(npoints, 3);
cnt=0;
for xi=-1:1
    for yi=-1:1
        for zi=-1:1
            cnt=cnt+1;
            point3d(cnt,:)=[xi,yi,zi];
        end
    end
end

cam2d=zeros(npoints,2,ncams);


%get 2D positions in pixel plane coordinates (so that Noff and Npix
%parameters do not have to be passed through the optimization)
for icam = 1:ncams
  %To get exact results we first need to project the input rotation matrix
  %onto a true rotation matrix--the input from calib_Tsai is typically not quite normalized
      camParaCalib(icam).R = gv_angles2rotmat(gv_rotmat2angles(camParaCalib(icam).R));
      camParaCalib(icam).Rinv = inv(camParaCalib(icam).R);
      camParaCalib(icam).T= -1* camParaCalib(icam).R * camParaCalib(icam).Tinv; % well keep Tinv and make T match Tinv
      cam2d(:,:,icam)=calibProj_Tsai(camParaCalib(icam), point3d);
      cam2d(:,1,icam)=(cam2d(:,1,icam)-camParaCalib(icam).Npixw/2-camParaCalib(icam).Noffw)*camParaCalib(icam).wpix;
      cam2d(:,2,icam)=(-cam2d(:,2,icam)+camParaCalib(icam).Npixh/2-camParaCalib(icam).Noffh)*camParaCalib(icam).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
end

%Now stuff calibration data into an array needed by fminsearch
%calfixed is the fixed camera.  calopt is the remaining cameras whos
%parameters will be optimized
%If a camera other than 1 is fixed, then gv_dynamic_fitfunc must be changed
%also.
% calfixed=zeros(8,1);
% calfixed(1:3)=gv_rotmat2angles(camParaCalib(1).R);
% calfixed(4:6)=camParaCalib(1).T;
% calfixed(7)=camParaCalib(1).f_eff;
% calfixed(8)=camParaCalib(1).k1;
% 
% calopt=zeros(8,2);
% calopt(1:3,1)=gv_rotmat2angles(camParaCalib(2).R);
% calopt(4:6,1)=camParaCalib(2).T;
% calopt(7,1)=camParaCalib(2).f_eff;
% calopt(8,1)=camParaCalib(2).k1;
% 
% calopt(1:3,2)=gv_rotmat2angles(camParaCalib(3).R);
% calopt(4:6,2)=camParaCalib(3).T;
% calopt(7,2)=camParaCalib(3).f_eff;
% calopt(8,2)=camParaCalib(3).k1;

calinitial=zeros(8,3);
for icam = 1:ncams
  calinitial(1:3,icam)=gv_rotmat2angles(camParaCalib(icam).R);
  calinitial(4:6,icam)=camParaCalib(icam).T;
  calinitial(7,icam)=camParaCalib(icam).f_eff;
  calinitial(8,icam)=camParaCalib(icam).k1;
end

%The parameters not in calopt are fixed.  Changing this section requires
%matching changes in gv_dynamic_fitfunc
calfixed=calinitial;
calopt=zeros(6,1);
for icam = 2:2
  calopt(1:3,icam-1)=gv_rotmat2angles(camParaCalib(icam).R);
  calopt(4:6,icam-1)=camParaCalib(icam).T;
%  calopt(7, icam-1)=camParaCalib(2).f_eff;
 % calopt(7,icam-1)=camParaCalib(icam).k1;
end


calorig = calopt;
gv_dynamic_fitfunc(calopt,calfixed,cam2d)
calopt(1,:)=calopt(1,:)-0.05;
calopt(2,:)=calopt(2,:)+0.05;
calopt(3,:)=calopt(3,:)+0.05;
calopt(4,:)=calopt(4,:)-1;
calopt(5,:)=calopt(5,:)+1;
calopt(6,:)=calopt(6,:)+5;
%calfixed(4,3)=calfixed(4,3)+0.2
%calfixed(5,3)=calfixed(5,3)-0.2
gv_dynamic_fitfunc(calopt,calfixed,cam2d)


fmin_options.Display='iter';
fmin_options.MaxFunEvals=4000;
fmin_options.TolFun=1e-15;
fmin_options.TolX=1e-6;
calout = fminsearch(@(x) gv_dynamic_fitfunc(x,calfixed,cam2d), calopt, fmin_options);

calorig
calopt
calout


