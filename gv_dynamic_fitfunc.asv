function deviation = gv_dynamic_fitfunc(calopt, calarray, cam2d)

% This function is for use with dynamic calibration.  It combines the calibration
% parameters into a single array and then calls the routine to calculate the distances
%  between the rays corresponding to given 2d image plane coordinates.  
%
% inputs:
%   calopt       --  parameters that are being optimized.  Since fminsearch
%                     wants a matrix containing only the optimized
%                     parameters, these are passed separately and stuffed back
%                     into the full calibration matrix in this function.
%   calarray     --  full 8 by 3 calibration matrix.  Columns are for each camera.  1:3 are angles, 4:6 are T, 7 is effective focal length (f_eff) and 8 is distortion (k1) 
%   cam2d        --  N by 2 by 3 array containing 2d image plane coordinates of
%                       N matched particles.   Last index is camera id.  Middle index is coordinate (h or w).
%
% outputs:
%   deviation    -- average of the 3D ray mismatch. 
%
%  called functions:
%    gv_calc_ray_mismatch

calarray(1:6,2:3)=calopt(1:6,1:2);   %here we assume camera 1 is fixed.  If this changes, it needs to be changed here also.

ray3mismatch=gv_calc_ray_mismatch(calarray, cam2d);

deviation = mean(ray3mismatch);




