function camParaCalib = checkalign(filepath, calibimgsize, imgsize, ncams)
% check alignment of cameras to see if further tuning is needed
%
% Inputs:
%   filepath    --  path to particle center files
%   calibimgsize     --  size of the calibration images, [Nx, Ny]
%   imgsize     --  size of the images that will be used in experiments, [Nx, Ny]
%   ncams       --  number of cameras
%
% Outputs:
%   No output in current version.
%

% example: 
%
% filepath = '/Users/haitao/mnt/raid/hd2/calibration/070111';
% calibimgsize = [800 600];
% imgsize = [256 256];
% ncams = 3;


if filepath(end) ~= '\'
    filepath = strcat(filepath, '\');
end

% known parameter of the Phantom cameras
camParaknown.Npixh = 256;
camParaknown.Npixw = 256;
camParaknown.hpix = 0.022;	% pixel size (mm)
camParaknown.wpix = 0.022;	% pixel size (mm)

for icam = 1:ncams
	x = myload(strcat(filepath,'calibpointspos.cam', num2str(icam-1), '.dat'), '#');
	[camParaCalib(icam) err Xout] = calib_Tsai(x(:,1:2), x(:,3:5), camParaknown, calibimgsize);
	X(icam).Ximg = Xout(:,1);
	X(icam).Yimg = Xout(:,2);
	X(icam).x3D = Xout(:,3);
	X(icam).y3D = Xout(:,4);
	X(icam).z3D = Xout(:,5);
	err_img(icam).err_x = err(:,1);
	err_img(icam).err_y = err(:,2);
end

% calculate the projection of camera axes on other cameras' image plane
Zmax = 1000;
% projection of other camera's axis on one camera's image plane
Ndisp = 100;
xorig_proj = zeros(ncams, ncams);
yorig_proj = zeros(ncams, ncams);
xaxes_proj = zeros(Ndisp, ncams, ncams);
yaxes_proj = zeros(Ndisp, ncams, ncams);
for i = 1:ncams
    for j = 1:ncams
        if( j ~= i)
            % origin of cam j's center on cam i's image plane
            P100 = camParaCalib(i).R * camParaCalib(j).Tinv + camParaCalib(i).T;
            x100 = camParaCalib(i).f_eff/P100(3)*P100(1);
            y100 = camParaCalib(i).f_eff/P100(3)*P100(2);
            xorig_proj(i, j) = x100;
            yorig_proj(i, j) = y100;
            % one point on cam 1's axis on cam 0's image plane
            P101 = camParaCalib(i).R * (camParaCalib(j).Rinv * [0 0 Zmax]' + camParaCalib(j).Tinv) + camParaCalib(i).T;
            x101 = camParaCalib(i).f_eff/P101(3)*P101(1);
            y101 = camParaCalib(i).f_eff/P101(3)*P101(2);
            xaxes_proj(:,i,j) = reshape(linspace(x100, x101, Ndisp), Ndisp, 1);
            yaxes_proj(:,i,j) = reshape(linspace(y100, y101, Ndisp), Ndisp, 1);
        end
    end
end

% project some points on 3 cameras
L = 30;
Np = 10000;
xp = (rand(Np,1)-0.5)*L;
yp = (rand(Np,1)-0.5)*L;
zp = (rand(Np,1)-0.5)*L;
% check laser beam illumination
nlaser = 2;
laserdir = zeros(nlaser, 3);
laserdir(1,:) = [1 0 0];
laserdir(2,:) = [cos(pi/8) -sin(pi/8) 0];
Rlaser = [200 200];
a = [xp yp zp];
ind = logical(zeros(Np,1));
for i=1:nlaser
    laserdir(i,:) = laserdir(i,:)/sqrt(sum(laserdir(i,:).^2));
    d = zeros(Np,1);
    for j=1:Np
        d(j) = sum((a(j) - sum(a(j).*laserdir(i,:))*laserdir(i,:)).^2);
    end
    ind = ind | ( d <= (Rlaser(i))^2 );
end
a = a(ind,:);
xp = a(:,1);
yp = a(:,2);
zp = a(:,3);
Np = length(xp);

xPimg = zeros(Np,ncams);
yPimg = zeros(Np,ncams);
imglimit_512 = [[-256 256]*camParaknown.wpix [-256 256]*camParaknown.hpix];
imglimit_256 = [[-128 128]*camParaknown.wpix [-128 128]*camParaknown.hpix];

for i = 1:Np
    for icam = 1:ncams
        p0 = camParaCalib(icam).R * ([xp(i) yp(i) zp(i)]') + camParaCalib(icam).T;
        xPimg(i, icam) = camParaCalib(icam).f_eff/p0(3)*p0(1);
        yPimg(i, icam) = camParaCalib(icam).f_eff/p0(3)*p0(2);
    end
end
ind512 = logical(zeros(Np, ncams));
ind256 = logical(zeros(Np, ncams));
ind512_all = logical(ones(Np, 1));
ind256_all = logical(ones(Np, 1));
for icam = 1:ncams
    ind512(:,icam) = (xPimg(:,icam) > imglimit_512(1)) & (xPimg(:,icam) < imglimit_512(2)) & (yPimg(:,icam) > imglimit_512(3)) & (yPimg(:,icam) < imglimit_512(4));
    ind256(:,icam) = (xPimg(:,icam) > imglimit_256(1)) & (xPimg(:,icam) < imglimit_256(2)) & (yPimg(:,icam) > imglimit_256(3)) & (yPimg(:,icam) < imglimit_256(4));
    ind512_all = ind512_all & ind512(:,icam);
    ind256_all = ind256_all & ind256(:,icam);
end

% calculate the distance between camera axes
M = zeros(3,3);
pM = zeros(3,ncams);
u = zeros(3, ncams);
for icam = 1:ncams
    % a unit vector pointing along camera axes
    u(:,icam) = camParaCalib(icam).Rinv * [0 0 1]';
    u(:,icam) = u(:,icam) / sqrt(sum(u(:,icam) .* u(:,icam)));
    uM = eye(3) - u(:,icam) * (u(:,icam))';
    pM(:,icam) = uM * camParaCalib(icam).Tinv;
    M = M + uM;
end
% point of axes crossing
p = M \ sum(pM,2)  % sums pm x together for alll three cameras.  Makes a column vector
h = zeros(1, ncams);
for icam = 1:ncams
    temp = p - ((p') * u(:,icam)) * u(:,icam) - pM(:,icam);
    h(icam) = sqrt(temp' *temp);
end
hbar = sqrt(mean(h.*h));
[h hbar]


% efficiency at 512
eff512 = (sum(ind512_all)*ones(1,ncams))./sum(ind512, 1)
eff256 = (sum(ind256_all)*ones(1,ncams))./sum(ind256, 1)

dkgreen = [0 0.5 0];
brown = [0.5 0 0];
camcolor = {'r', 'g', 'b', 'm', 'c', 'k'};
figure;
% plot figures in 2 rows, the number of columns depends on ncams
nfrow = 2;
nfcol = ceil(ncams/nfrow);
for icam = 1:ncams
    subplot(nfrow, nfcol, icam)
    plot(xPimg(ind512_all, icam), yPimg(ind512_all, icam), '.', 'Color', dkgreen);
    hold on
    plot(xPimg(ind256_all, icam), yPimg(ind256_all, icam), '.', 'Color', brown);
    for j = 1:ncams
        if (j ~= icam)
            plot(xaxes_proj(:, icam, j), yaxes_proj(:, icam, j), '-.', 'Color', char(camcolor(j)));
            plot(xorig_proj(icam, j), yorig_proj(icam, j), 'o', 'Color', char(camcolor(j)));
        end
    end
    plotrect(imglimit_512, 'k-');
    plotrect(imglimit_256, 'k-');
    axis([-300 300 -300 300]*0.022);
    title(strcat('projection on camera ', num2str(icam-1)));
end
% plot3(xp(ind512), yp(ind512), zp(ind512), 'c.');
% hold on
% plot3(xp(ind256), yp(ind256), zp(ind256), 'm.');
% hold off
% plot(xp(ind256), yp(ind256), 'm.');
% title('particles seen by all 3 cameras in 3D');

% plot the projection of camera axes on xy, xz and yz plane
Zp_cam = [0 0 Zmax]';
Zp_world = zeros(3, ncams);
for icam = 1:ncams
    Zp_world(:,icam) = camParaCalib(icam).Rinv * Zp_cam + camParaCalib(icam).Tinv;
end
% camera viewing prism edges
prismEdges_256 = zeros(3, 4, ncams);
for icam = 1:ncams
    xcorner = imgsize(1)/2*camParaCalib(icam).wpix;
    ycorner = imgsize(2)/2*camParaCalib(icam).hpix;
    feff = camParaCalib(icam).f_eff;
    Zv1_cam = [xcorner ycorner feff]' * (Zmax/feff);
    Zv2_cam = [-xcorner ycorner feff]' * (Zmax/feff);
    Zv3_cam = [-xcorner -ycorner feff]' * (Zmax/feff);
    Zv4_cam = [xcorner -ycorner feff]' * (Zmax/feff);
    prismEdges_256(:, 1, icam) = camParaCalib(icam).Rinv * Zv1_cam + camParaCalib(icam).Tinv;
    prismEdges_256(:, 2, icam) = camParaCalib(icam).Rinv * Zv2_cam + camParaCalib(icam).Tinv;
    prismEdges_256(:, 3, icam) = camParaCalib(icam).Rinv * Zv3_cam + camParaCalib(icam).Tinv;
    prismEdges_256(:, 4, icam) = camParaCalib(icam).Rinv * Zv4_cam + camParaCalib(icam).Tinv;
end
% calculate the radius of the inscribe sphere
% approximated by the minimum distance to camera view prism surfaces
Rsphere = 1.E6;
for icam = 1:ncams
    Ocam = camParaCalib(icam).Tinv;
    A = prismEdges_256(:, :, icam);
    h1 = dist2plane(Ocam, A(:,1), A(:,2), p);
    h2 = dist2plane(Ocam, A(:,2), A(:,3), p);
    h3 = dist2plane(Ocam, A(:,3), A(:,4), p);
    h4 = dist2plane(Ocam, A(:,4), A(:,1), p);
    Rsphere = min([Rsphere, h1, h2, h3, h4]);
end
Rsphere
% magnification
for icam=1:ncams
    dist = sqrt(sum((camParaCalib(icam).T).^2));
    str=sprintf('cam %d each pixel corresponding to: [%6.3f %6.3f] mm in space', ...
                icam, camParaCalib(icam).wpix*dist/camParaCalib(icam).f_eff, ...
                        camParaCalib(icam).hpix*dist/camParaCalib(icam).f_eff);
    disp(str);
end

% a parameter along the lines
Ndisp = 500;
t = [0:1/Ndisp:1]';
Zp_disp = zeros(3, length(t), ncams);
prismEdgesdisp_256 = zeros(3, length(t), 4, ncams);
for icam = 1:ncams
    for i = 1:length(t)
        Zp_disp(:, i, icam) = (camParaCalib(icam).Tinv + t(i)*(Zp_world(:,icam) - camParaCalib(icam).Tinv))';
        prismEdgesdisp_256(:, i, 1, icam) = (camParaCalib(icam).Tinv + t(i)*(prismEdges_256(:,1,icam) - camParaCalib(icam).Tinv))';
        prismEdgesdisp_256(:, i, 2, icam) = (camParaCalib(icam).Tinv + t(i)*(prismEdges_256(:,2,icam) - camParaCalib(icam).Tinv))';
        prismEdgesdisp_256(:, i, 3, icam) = (camParaCalib(icam).Tinv + t(i)*(prismEdges_256(:,3,icam) - camParaCalib(icam).Tinv))';
        prismEdgesdisp_256(:, i, 4, icam) = (camParaCalib(icam).Tinv + t(i)*(prismEdges_256(:,4,icam) - camParaCalib(icam).Tinv))';
    end
end

% projection of camera axes and edges of viewing zone on to planes in 3D space
figure;
markersize = 6;
subplot(131);
for icam = 1:ncams
    plot(Zp_disp(1,:,icam), Zp_disp(2,:,icam), '-.', 'Color', char(camcolor(icam)));
    hold on;
    plot(prismEdgesdisp_256(1,:,1,icam), prismEdgesdisp_256(2,:,1,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(1,:,2,icam), prismEdgesdisp_256(2,:,2,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(1,:,3,icam), prismEdgesdisp_256(2,:,3,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(1,:,4,icam), prismEdgesdisp_256(2,:,4,icam), '-', 'Color', char(camcolor(icam)));
end
plot(p(1), p(2), 'k+', 'MarkerSize', markersize);
%circle(p(1), p(2), Rsphere, 'k--');
hold off
title('projection on to XY');
subplot(132);
for icam = 1:ncams
    plot(Zp_disp(1,:,icam), Zp_disp(3,:,icam), '-.', 'Color', char(camcolor(icam)));
    hold on;
    plot(prismEdgesdisp_256(1,:,1,icam), prismEdgesdisp_256(3,:,1,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(1,:,2,icam), prismEdgesdisp_256(3,:,2,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(1,:,3,icam), prismEdgesdisp_256(3,:,3,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(1,:,4,icam), prismEdgesdisp_256(3,:,4,icam), '-', 'Color', char(camcolor(icam)));
end
plot(p(1), p(3), 'k+', 'MarkerSize', markersize);
%circle(p(1), p(3), Rsphere, 'k--');
hold off
title('projection on to XZ');
subplot(133);
for icam = 1:ncams
    plot(Zp_disp(2,:,icam), Zp_disp(3,:,icam), '-.', 'Color', char(camcolor(icam)));
    hold on;
    plot(prismEdgesdisp_256(2,:,1,icam), prismEdgesdisp_256(3,:,1,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(2,:,2,icam), prismEdgesdisp_256(3,:,2,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(2,:,3,icam), prismEdgesdisp_256(3,:,3,icam), '-', 'Color', char(camcolor(icam)));
    plot(prismEdgesdisp_256(2,:,4,icam), prismEdgesdisp_256(3,:,4,icam), '-', 'Color', char(camcolor(icam)));
end
plot(p(2), p(3), 'k+', 'MarkerSize', markersize);
%circle(p(2), p(3), Rsphere, 'k--');
hold off
title('projection on to YZ');

% plot image plane error statistics
binsize = 0.01;
bin = [-2:binsize:2];
figure;
for icam = 1:ncams
	errx = err_img(icam).err_x;
	erry = err_img(icam).err_y;
	Px = hist(errx, bin)/length(errx)/binsize;
	Py = hist(erry, bin)/length(erry)/binsize;
	subplot(ncams, 1, icam);
	plot(bin, Px, 'r-', bin, Py, 'b-');
	title(strcat('image plane error for cam', num2str(icam-1)));
	h = legend(sprintf('x: \\mu_x = %8.5f, rms = %8.5f', mean(errx), camParaCalib(icam).err_x), ...
		sprintf('y: \\mu_y = %8.5f, rms = %8.5f', mean(erry), camParaCalib(icam).err_y));
	set(h, 'FontSize', 10);
	xlabel('Pixel');
	ylabel('PDF');
	% write projected particle coordinates to files for later access
	fid = fopen(strcat('pointsimgpos.cam',num2str(icam-1),'.dat'), 'w');
	fprintf(fid, '%14.7f\t%14.7f\t%14.7f\t%14.7f\n', [X(icam).Ximg X(icam).Yimg errx erry]');
	fclose(fid);
end


% plot 3D error
% first find matches, since we know both the projection and 3D information of particles, matching is easier
xstep = 1.27;	% x-distance between dots
ystep = 2;	% y-distance between dots
zstep = 2;	% z-distance
Nx = 5;
Ny = 9;
Nz = 9;
ind = zeros(ncams+1, Nx, Ny, Nz);
Nmax = 0;
for icam = 1:ncams
	x3D = X(icam).x3D;
	y3D = X(icam).y3D;
	z3D = X(icam).z3D;
	N = length(x3D);
	Nmax = max(N, Nmax);
	for n = 1:N
		i = floor(x3D(n)/xstep + floor(Nx+1)/2 + 0.5);
		j = floor(y3D(n)/ystep + floor(Ny+1)/2 + 0.5);
		k = floor(z3D(n)/zstep + floor(Nz+1)/2 + 0.5);
		ind(icam, i, j, k) = n;
		ind(ncams+1, i, j, k) = ind(ncams+1, i, j, k) + 1;  %checking if all cameras see this particle--if so this will end up equal to ncams.
	end
end

Xp3D = zeros(Nmax,3);
Ximg = zeros(Nmax,ncams*2);
N3D = 0;
for i = 1:Nx
for j = 1:Ny
for k = 1:Nz
	if (ind(ncams+1, i, j, k) == ncams)
		N3D = N3D + 1;
		i1 = ind(1, i, j, k);
		Xp3D(N3D,:) = [X(1).x3D(i1) X(1).y3D(i1) X(1).z3D(i1)];
		for icam = 1:ncams
			n = ind(icam, i, j, k);
			Ximg(N3D, (icam-1)*2+1:(icam-1)*2+2) = [X(icam).Ximg(n) X(icam).Yimg(n)];
		end
	end
end
end
end
Xp3D = Xp3D(1:N3D,:);
Ximg = Ximg(1:N3D,:);
[xmatch, dist3D] = matchpoints(ncams, camParaCalib, Ximg);
errx = xmatch(:,1) - Xp3D(:,1);
erry = xmatch(:,2) - Xp3D(:,2);
errz = xmatch(:,3) - Xp3D(:,3);
errr = sqrt(errx.^2+erry.^2+errz.^2);

binsize = 0.05;
pixsize = Rsphere*2/min(imgsize);
bin = [-1:binsize:1]*pixsize;
Px = hist(errx, bin)/length(errx)/(binsize*pixsize);
Py = hist(erry, bin)/length(erry)/(binsize*pixsize);
Pz = hist(errz, bin)/length(errz)/(binsize*pixsize);
Pr = hist(errr, bin+pixsize)/length(dist3D)/(binsize*pixsize);
Ph = hist(dist3D, pixsize+bin)/length(dist3D)/(binsize*pixsize);
figure;
subplot(211);
plot(bin, Px, 'r-', bin, Py, 'b-', bin, Pz, 'g-');
title('3D error');
h = legend(sprintf('x: \\mu_x = %8.5f, rms = %8.5f', mean(errx), sqrt(mean(errx.^2))), ...
	sprintf('y: \\mu_y = %8.5f, rms = %8.5f', mean(erry), sqrt(mean(erry.^2))), ...
	sprintf('z: \\mu_z = %8.5f, rms = %8.5f', mean(errz), sqrt(mean(errz.^2))));
set(h, 'FontSize', 10);
xlabel('mm');
ylabel('PDF');
subplot(212);
plot(bin+pixsize, Pr, 'r-', bin+pixsize, Ph, 'b-');
title('3D distance in matching and 3D error magnitude');
h = legend(sprintf('errr: \\mu = %8.5f, rms = %8.5f', mean(errr), sqrt(mean(errr.^2))), ...
	sprintf('h: \\mu = %8.5f, rms = %8.5f', mean(dist3D), sqrt(mean(dist3D.^2))) );
set(h, 'FontSize', 10);
xlabel('mm');
ylabel('PDF');

% write projected particle coordinates to files for later access
fid = fopen(strcat('points3Dpos.dat'), 'w');
fprintf(fid, '%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\t%14.7e\n', ...
	[Xp3D xmatch dist3D Ximg]');
fclose(fid);

% Now write the configuration file
fid = fopen(strcat(char(filepath), 'PTVSetup_', num2str(imgsize(1)), '.cfg'), 'w');
fprintf(fid, '# PTV experiment configuration file\n');
fprintf(fid, '\n %d\t# ncams\n\n', ncams);
for icam = 1:ncams
	fprintf(fid, '######## cam %d ########\n\n', icam-1);
	fprintf(fid, '%d\t\t\t# Npix_x\n', imgsize(1));
	fprintf(fid, '%d\t\t\t# Npix_y\n', imgsize(2));
	fprintf(fid, '%11.8f\t\t\t# pixsize_x (mm)\n', camParaCalib(icam).wpix);
	fprintf(fid, '%11.8f\t\t# pixsize_y (mm)\n', camParaCalib(icam).hpix);
	fprintf(fid, '%12.8f\t\t# effective focal length (mm)\n', camParaCalib(icam).f_eff);
	% Note the sign change for k1, because its meaning is different in calib_Tsai and the stereomatching code
	fprintf(fid, '%15.8e\t# radial distortion k1 (1/pixel)\n', (camParaCalib(icam).k1));	
	fprintf(fid, '%15.8e\t# tangential distortion p1 (1/pixel)\n', camParaCalib(icam).p1);
	fprintf(fid, '%15.8e\t# tangential distortion p2 (1/pixel)\n', camParaCalib(icam).p2);
	fprintf(fid, '0.0\t\t# x0 (mm)\n');
	fprintf(fid, '0.0\t\t# y0 (mm)\n');
    % Parameters for camera movement correction
	fprintf(fid, '0.0\t\t# x0_offset (pixel)\n');
	fprintf(fid, '0.0\t\t# y0_offset (pixel)\n');
	fprintf(fid, '1.0\t\t# x_rot (cos(theta))\n');
	fprintf(fid, '0.0\t\t# y_rot (sin(theta))\n');
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
	fprintf(fid, '\n');
end
% laser beam parameters
fprintf(fid, '##### laser beam #####\n');
fprintf(fid, '\n');
fprintf(fid, '0\t\t# finite volume illumination\n');
fprintf(fid, '0\t\t# illum_xmin\n');
fprintf(fid, '0\t\t# illum_xmax\n');
fprintf(fid, '0\t\t# illum_ymin\n');
fprintf(fid, '0\t\t# illum_ymax\n');
fprintf(fid, '0\t\t# illum_zmin\n');
fprintf(fid, '0\t\t# illum_zmax\n');
fprintf(fid, '\n');
% 3D matching parameters
fprintf(fid, '#### thresholds for 3D matching ####\n');
fprintf(fid, '\n');
fprintf(fid, '1.5\t\t# mindist_pix (pixel)\n');
fprintf(fid, '0.1\t\t# maxdist_3D (mm)\n');
fclose(fid);
