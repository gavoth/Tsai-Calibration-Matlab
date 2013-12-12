% This code optimizes the parameters for a camera calibration based
% on the Tsai model.  It starts with an approximate calibration (currently
% from a calibpointspos.cam file created by PTVsetupPrep.m), and then reads 
% in a dist3D.dat file containing 2D pixel coordinates 
% of matched particles.  It then does a nonlinear optimization to find the 
% calibration parameters that minimize the mismatch between rays of matched particles.  
% Currently we leave the first camera
% where it was calibrated and adjust the position and rotation matrix for the
% other two cameras to minimize the mean ray intersection distance. 
%
% inputs:
%   calibpointspos.cam       --  file containing positions of mask points
%                           that are used to find the approximate
%                           calibration.
%   dist3D.dat     --  file containing 2D pixel coordinates of matched
%                       particles. 
% hardcoded parameters
%   dist3D_path and dist3D_stem -- designates the dist3D file to use.
%   npointsfit  --sets the number of 2D positions used for the
%           optimization.  More points gives a better fit, but it is slower.
%           Something like 200-800 works pretty well.
%   MaxFunEvals  -- Sets the maximum number of evaluations during the
%           minimization.   Larger numbers make it slower but more likely to
%           find the best minimum.  Something like 300 works pretty well for me, but will
%           depend on your data.
%
% outputs:
%   .cfg file (PTVSetup.cfg) determined by fname_cfg at the end of the code
%
%  called functions:
%   calib_Tsai  (with subfunctions 'myload' and 'importfile') 
%   gv_rotmat2angles  -- converts a rotation matrix to 3 angles.
%   gv_angles2rotmat  -- converts 3 angles to a rotation matrix.
%   gv_dynamic_fitfunc --repacks parameters and calls gv_calc_ray_mismatch
%   gv_calc_ray_mismatch  -- calculates the distances between lines of sight
%   gv_imgplane2unitvector -- called by gv_calc_ray_mismatch
%   gv_write_calib_cfg  -- writes the PTVSetup.cfg file
%
% Work left to do:
%   -- make the option of loading the initial calibration from a PTVSetup.cfg file rather
%   than using a calibpointspos.cam file to generate it
%   -- store the calibration error in the same formate that PTVSetup.cfg
%   previoulsy used;  currently it stuffs the camera's ray mismatch in mm into both err_x and err_y while leaving err_t unchanged 
%   -- make the option of using a better minimization routine if the user
%   has access to the optimization toolbox.
%   --ultimately it would be nice if the track files contained the 2D pixel
%   coordinates that created that 3D particle.  Then this code could just
%   read in a track file and optimize the calibration.
%   -- It would also be nice to clean up the calarray data passing so that
%   this code and gv_dynamic_fitfunc do not have to be hacked to change the number of cameras or which
%   parameters are optimized.  
%
% calarray:
%   contains the primary calibration data for all 3 cameras.  Format is 8 x ncams matrix with the 8 entries for each camera: 
%  3 angles for rotation matrix (R), 3 positions for position of camera (T), effective_focal length (f_eff), distortion (k1) 
%
%  If you are new to the Tsai calibration, I recommend trying the
%  script gv_tsai_calibration_sandbox.m
%
% Author:  Greg Voth, Wesleyan University, Spring 2008, gvoth@wesleyan.edu

ncams=3;
calibimgsize=[256 256];


%calibpoints_filepath='C:\Documents and Settings\gvoth\My Documents\data\cornell\calib\jan16\calibpoints\';
calibpoints_filepath='C:\Documents and Settings\gvoth\My Documents\data\cornell\calib\jan07\3positions\';

dist3D_path='C:\Documents and Settings\gvoth\My Documents\data\cornell\dist3D\';
dist3D_stem='ptv5_20kHz_55um_0_dist';
%dist3D_stem='dist3D_ptv3_222um_525Hz_0';
dist3D_fname=strcat(dist3D_path,dist3D_stem, '.dat')
%dist3D_fname=strcat(dist3D_path,dist3D_stem, '.dat')

camParaknown.Npixh = 256;
camParaknown.Npixw = 256;
camParaknown.hpix = 0.022;	% pixel size (mm)
camParaknown.wpix = 0.022;	% pixel size (mm)
for icam = 1:ncams
	x = myload(strcat(calibpoints_filepath,'calibpointspos.cam', num2str(icam-1), '.dat'), '#');
	[camParaCalib(icam) err Xout] = calib_Tsai(x(:,1:2), x(:,3:5), camParaknown, calibimgsize);
end

%fname_cfg=strcat(dist3D_path,'PTVSetup_firstguess_',dist3D_stem,'.cfg');
%gv_write_calib_cfg(camParaCalib, ncams, fname_cfg);


%for our data with 200mm lenses and 1cm field of view, the radial
%distortion does not seem to affect the calibration quality at all, so it
%is set to zero to simplify things.  This code should work fine
%with these four lines commented out.  In the current form, it would leave
%the distortion as a fixed parameter, but it could be optimized as well.
for icam=1:ncams
   camParaCalib(icam).k1=0;
   camParaCalib(icam).k1star=0;
end

 
%Read in the dist3D files containing the matched camera positions from all
%3 cameras

%S = uiimport(dist3D_fname)
dist3D = importdata(dist3D_fname);
%The dist3D files contain x3d, y3d, z3d, matchingerror, xcam1, ycam1,
%xcam1proj, ycam1proj, xcam2, ycam2, xcam2proj, ycam2proj, xcam3, ycam3,
%xcam3proj, ycam3proj   

[n2d,ncols]=size(dist3D);

%Plot the histogram of mismatches as determined by the matching code that
%created the dist3D file.
figure(6)
hist(dist3D(:,4),80);
title('Matching error in dist3D data');
xlabel('error (mm)');


%Pick the points to use for the fit.  Because the data is tracks, the
%positions are highly correlated and it is best to choose positions spaced
%throughout the file.  
npointsfit=300;
step = fix(n2d/npointsfit);
chosen=[1:step:n2d];
dist3D_chosen=dist3D(chosen(1:npointsfit),:);  %Not very elegant--it just needs to select npointsfit elements from throughout dist3D 


cam2d=zeros(npointsfit,2,ncams);

%Choose another set of points to check the calibration against at the end.
chosen_check=[10:step:n2d];
dist3D_chosen_check=dist3D(chosen_check(1:npointsfit),:);  
cam2d_check=zeros(npointsfit,2,ncams);

%repack the 2d camera coordinate data into a cam2d array.  
for icam = 1:ncams
    cam2d(:,:,icam)=dist3D_chosen(:,5+(icam-1)*4:6+(icam-1)*4);
    cam2d_check(:,:,icam)=dist3D_chosen_check(:,5+(icam-1)*4:6+(icam-1)*4);
end

%show the particles from camera 1
figure(1)
plot(cam2d(:,1,1),cam2d(:,2,1), '.g');
title('Points used as seen by the first camera')
xlabel('x position (pixels)')
ylabel('y position (pixels)')

%Now camera coordinates need to be scaled and shifted
for icam = 1:ncams
      cam2d(:,1,icam)=(cam2d(:,1,icam)-camParaCalib(icam).Npixw/2-camParaCalib(icam).Noffw)*camParaCalib(icam).wpix;
      cam2d(:,2,icam)=(-cam2d(:,2,icam)+camParaCalib(icam).Npixh/2-camParaCalib(icam).Noffh)*camParaCalib(icam).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
      cam2d_check(:,1,icam)=(cam2d_check(:,1,icam)-camParaCalib(icam).Npixw/2-camParaCalib(icam).Noffw)*camParaCalib(icam).wpix;
      cam2d_check(:,2,icam)=(-cam2d_check(:,2,icam)+camParaCalib(icam).Npixh/2-camParaCalib(icam).Noffh)*camParaCalib(icam).hpix;  %vertical coordinate needed to be switched in sign to make it work--I thought the relection in the rotation matrix took care of this--but it doesn't work without this sign negative
end

%Now stuff the calibration data into an the arrays needed by fminsearch.
%calinitial and calarray contain the calib parameters for all 3 cameras.
%fminsearch wants an array with only the optimization parameters, so we pass
%those separately, and then stuff the optimizing parameters back into the full
%calarray inside gv_dynamic_fitfunc.  If you change which parameters are
%optimized you need to change gv_dynamic_fitfunc so that it stuffs the
%correct array elements.

calinitial=zeros(8,ncams);
for icam = 1:ncams
  calinitial(1:3,icam)=gv_rotmat2angles(camParaCalib(icam).R);
  calinitial(4:6,icam)=camParaCalib(icam).T;
  calinitial(7,icam)=camParaCalib(icam).f_eff;
  calinitial(8,icam)=camParaCalib(icam).k1;
end
calarray=calinitial;

     %plot initial distribution of mismatches--not necessary, but
     %convenient to do here.
     [dist3, dist1,points3D]=gv_calc_ray_mismatch(calarray,cam2d);
        figure(2);
        hist(dist3,50);
        title('initial mismatch distribution using initial calibration');
        xlabel('mismatch (mm)')
        allcams_initial_mismatches=mean(dist3)
        
 dist3D_chosen(1:10,1:4) 
 points3D(1:10,:)
 
%The parameters not in calopt are fixed.  Currently we only fit R and T,
%the rotation matrix and the position of the cameras (the rotation matrix
%is represented by 3 angles).
calopt=zeros(6,2);
for icam=2:3
  calopt(1:3,icam-1)=gv_rotmat2angles(camParaCalib(icam).R);
  calopt(4:6,icam-1)=camParaCalib(icam).T;
  %calopt(7, 1)=camParaCalib(icam).f_eff;
 %calopt(8,1)=camParaCalib(icam).k1;
 end

%gv_dynamic_fitfunc(calopt,calfixed,cam2d)
fmin_options.Display='iter'; %change to 'final' to display only the final mean distance.
fmin_options.MaxFunEvals=600;
calout = fminsearch(@(x) gv_dynamic_fitfunc(x,calarray,cam2d), calopt, fmin_options);

calarray(1:6,2:3)=calout(1:6,1:2);
[dist3, dist1]=gv_calc_ray_mismatch(calarray,cam2d);

%display results of initial optimization
figure(3);
hist(dist3,50);
title('mismatch distribution after first optimization');
xlabel('mismatch (mm)')
allcams=mean(dist3)
cam1mismatch=mean(dist1(1,:))
cam2mismatch=mean(dist1(2,:))
cam3mismatch=mean(dist1(3,:))

%choose good matches
good_matches=(dist3 < 0.04);
cam2d_good=cam2d(good_matches,:,:);


  calopt=calarray(1:6,2:3);
  %By using calopt=calarray(1:8,2:3) here (and using gv_tmp_dynamic_fitfunc) this
  %can optimize all 8 parameters, but it doesn't give any better fit, so I
  %am using the simpler 12 parameter fit as above but without the bad matches.

calout = fminsearch(@(x) gv_dynamic_fitfunc(x,calarray,cam2d_good), calopt, fmin_options);
calarray(1:6,2:3)=calout(1:6,1:2);

%Check the quality of the final calibration
[dist3, dist1]=gv_calc_ray_mismatch(calarray,cam2d_good);

figure(4);
hist(dist3,50);
title('mismatch distribution after optimization on good matches');
xlabel('mismatch (mm)')
allcams=mean(dist3)
cammismatch(1)=mean(dist1(1,:))
cammismatch(2)=mean(dist1(2,:))
cammismatch(3)=mean(dist1(3,:))


[dist3, dist1]=gv_calc_ray_mismatch(calarray,cam2d_check);

figure(5);
hist(dist3,50);
title('final mismatch distribution of data not used for the optimization');
xlabel('mismatch (mm)')
allcams_check_with_mismatches=mean(dist3)


% Refill the camParaCalib structure and write it to a .cfg file.
%  (There turns out to be a very small change in the first camera even though it was not
%    dynamically calibrated since the rotation matrix is projected onto a matrix with 
%     determinant exactly -1, so it changes just a bit.  The new R needs to be kept since the 
%     dynamic calibration was done for this rotation matrix)
for icam = 1:ncams
  camParaCalib(icam).R=gv_angles2rotmat(calarray(1:3,icam));
  camParaCalib(icam).Rinv=inv(camParaCalib(icam).R);
  camParaCalib(icam).T=calarray(4:6,icam);
  camParaCalib(icam).Tinv=camParaCalib(icam).Rinv * (-1* camParaCalib(icam).T);
  camParaCalib(icam).f_eff = calarray(7,icam);
  camParaCalib(icam).k1 = calarray(8,icam);
  camParaCalib(icam).err_x=cammismatch(icam);
  camParaCalib(icam).err_y=cammismatch(icam);
end

fname_cfg=strcat(dist3D_path,'PTVSetup_optimized_',dist3D_stem,'.cfg');
%gv_write_calib_cfg(camParaCalib, ncams, fname_cfg);

