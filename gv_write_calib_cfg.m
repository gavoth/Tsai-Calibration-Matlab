function gv_write_calib_cfg(camParaCalib, ncams, fname)

% Now write the configuration file
fid = fopen(fname, 'w');
fprintf(fid, '# PTV experiment configuration file\n');
fprintf(fid, '\n %d\t# ncams\n\n', ncams);
for icam = 1:ncams
	fprintf(fid, '######## cam %d ########\n\n', icam-1);
	fprintf(fid, '%d\t\t\t# Npix_x\n', camParaCalib(icam).Npixw);
	fprintf(fid, '%d\t\t\t# Npix_y\n', camParaCalib(icam).Npixh);
	fprintf(fid, '%11.8f\t\t\t# pixsize_x (mm)\n', camParaCalib(icam).wpix);
	fprintf(fid, '%11.8f\t\t# pixsize_y (mm)\n', camParaCalib(icam).hpix);
	fprintf(fid, '%12.8f\t\t# effective focal length (mm)\n', camParaCalib(icam).f_eff);
	% Note the sign change for k1, because its meaning is different in calib_Tsai and the stereomatching code
	fprintf(fid, '%15.8e\t\t# radial distortion k1 (1/pixel)\n', -(camParaCalib(icam).k1));	
	fprintf(fid, '%15.8e\t\t# tangential distortion p1 (1/pixel)\n', camParaCalib(icam).p1);
	fprintf(fid, '%15.8e\t\t# tangential distortion p2 (1/pixel)\n', camParaCalib(icam).p2);
	fprintf(fid, '0.0\t\t\t# x0\n');
	fprintf(fid, '0.0\t\t\t# y0\n');
	fprintf(fid, '# rotation matrix R\n');
	fprintf(fid, '%12.8f\t%12.8f\t%12.8f\n', (camParaCalib(icam).R)');
	fprintf(fid, '# translation vector T\n');
	fprintf(fid, '%15.8f\n', camParaCalib(icam).T);
	fprintf(fid, '# inverse rotation matrix Rinv\n');
	fprintf(fid, '%12.8f\t%12.8f\t%12.8f\n', (camParaCalib(icam).Rinv)');
	fprintf(fid, '# inverse translation vector Tinv\n');
	fprintf(fid, '%15.8f\n', camParaCalib(icam).Tinv);
	fprintf(fid, '# rms distance between particle centers found on image plane and their projections\n');
	fprintf(fid, '# %15.8f\t\t# err_x\n', camParaCalib(icam).err_x);
	fprintf(fid, '# %15.8f\t\t# err_y\n', camParaCalib(icam).err_y);
	fprintf(fid, '# %15.8f\t\t# err_t\n', camParaCalib(icam).err_t);
	
end

fprintf(fid, '\n');
fprintf(fid, '\n');
fprintf(fid, '##### laser beam #####\n');
fprintf(fid, '\n');
fprintf(fid, '0\t# finite volume illumination\n');
fprintf(fid, '-500.0\t# illum_xmin\n');
fprintf(fid, '500.0\t# illum_xmax\n');
fprintf(fid, '12.0\t# illum_ymin\n');
fprintf(fid, '16.5\t# illum_ymax\n');
fprintf(fid, '-100.0\t# illum_zmin\n');
fprintf(fid, '100.0\t# illum_zmax\n');
fprintf(fid, '\n');
fprintf(fid, '#### parameter for 3D matching ####');
fprintf(fid, '\n');
fprintf(fid, '1.0\t# mindist_pix (pixel)');
fprintf(fid, '\n');
fprintf(fid, '0.08\t# maxdist_3D (mm)');

fclose(fid);