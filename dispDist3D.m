function corrPar = dispDist3D(fname)
h = load(fname);
imgsize = 256;
[N M] = size(h)
ncams = floor(M/7);
dmax = 0.05;
binsize = 0.005;
nbins = 100;
X3D = h(:,1:3);		% true 3D coordinate in world coordinates
X2D_seen = zeros(N, 2*ncams);
X2D = zeros(N, 2*ncams);
vdist = zeros(N, 3*ncams);
hdist = zeros(N, 1);
corrPar = zeros(4, ncams);
for icam = 1:ncams
	X2D_seen(:, [1:2]+(icam-1)*2) = h(:,[4:5]+(icam-1)*7) - imgsize/2;
	X2D(:, [1:2]+(icam-1)*2) = h(:, [6:7]+(icam-1)*7) - imgsize/2;
	vdist(:, [1:3]+(icam-1)*3) = h(:,[8:10]+(icam-1)*7);
end
hdist = sqrt(sum(vdist.^2, 2)/ncams);
hbin = [0:binsize:max(hdist)];
Ph = hist(hdist, hbin);
ind = logical(hdist <= dmax);


clear h;
nf_PDF = figure;
subplot(224);
loglog(hbin, Ph, 'r-');
linestyle = {'r-', 'g-', 'b-', 'k-'};
for icam = 1:ncams	
	% plot PDF of displacement components
	xbinsize = min(binsize, (max(vdist(:,1+(icam-1)*3))-min(vdist(:,1+(icam-1)*3)))/nbins);
	xbin = [min(vdist(:,1+(icam-1)*3)):xbinsize:max(vdist(:,1+(icam-1)*3))];
	ybinsize = min(binsize, (max(vdist(:,2+(icam-1)*3))-min(vdist(:,2+(icam-1)*3)))/nbins);
	ybin = [min(vdist(:,2+(icam-1)*3)):ybinsize:max(vdist(:,2+(icam-1)*3))];
	zbinsize = min(binsize, (max(vdist(:,3+(icam-1)*3))-min(vdist(:,3+(icam-1)*3)))/nbins);
	zbin = [min(vdist(:,3+(icam-1)*3)):zbinsize:max(vdist(:,3+(icam-1)*3))];

	Px = hist(vdist(:,1+(icam-1)*3), xbin);
	Py = hist(vdist(:,2+(icam-1)*3), ybin);
	Pz = hist(vdist(:,3+(icam-1)*3), zbin);
	
	figure(nf_PDF);
	subplot(221);
	semilogy(xbin, Px, char(linestyle(icam)));
	hold on;
	indx = logical(abs(xbin) <= dmax);
	[Pmax, Ix] = max( Px(indx));
	xbin = xbin(indx);
	xmin = xbin(max(1, Ix-max(1, floor(0.02/xbinsize))));
	xmax = xbin(min(length(xbin), Ix+max(1, floor(0.02/xbinsize))));
	semilogy([xmin xmin], [1 Pmax], strcat(char(linestyle(icam)), '-'));
	semilogy([xmax xmax], [1 Pmax], strcat(char(linestyle(icam)), '-'));
	
	subplot(222);
	semilogy(ybin, Py, char(linestyle(icam)));
	hold on;
	indy = logical(abs(ybin) <= dmax);
	[Pmax, Iy] = max( Py(indy));
	ybin = ybin(indy);
	ymin = ybin(max(1, Iy-max(1, floor(0.02/ybinsize))));
	ymax = ybin(min(length(ybin), Iy+max(1, floor(0.02/ybinsize))));
	semilogy([ymin ymin], [1 Pmax], strcat(char(linestyle(icam)), '-'));
	semilogy([ymax ymax], [1 Pmax], strcat(char(linestyle(icam)), '-'));
	
	
	subplot(223);
	semilogy(zbin, Pz, char(linestyle(icam)));
	hold on;
	indz = logical(abs(zbin) <= dmax);
	[Pmax, Iz] = max( Pz(indz));
	zbin = zbin(indz);
	zmin = zbin(max(1, Iz-max(1, floor(0.02/zbinsize))));
	zmax = zbin(min(length(zbin), Iz+max(1, floor(0.02/zbinsize))));
	semilogy([zmin zmin], [1 Pmax], strcat(char(linestyle(icam)), '-'));
	semilogy([zmax zmax], [1 Pmax], strcat(char(linestyle(icam)), '-'));
	
	ind = logical((vdist(:,1+(icam-1)*3)>=xmin) & (vdist(:,1+(icam-1)*3)<=xmax) ...
			& (vdist(:,2+(icam-1)*3)>=ymin) & (vdist(:,2+(icam-1)*3)<=ymax) ...
			& (vdist(:,3+(icam-1)*3)>=zmin) & (vdist(:,3+(icam-1)*3)<=zmax) );
	x2D = X2D(ind,[1:2]+(icam-1)*2);
	x2D_seen = X2D_seen(ind, [1:2]+(icam-1)*2);
	
	n = length(x2D);
	% fit rotation and translation, model as solid body rotation
	A = zeros(2*n, 4);
	xc = mean(x2D, 1);
	xc_seen = mean(x2D_seen, 1);
	x = x2D(:,1) - xc(1);
	y = x2D(:,2) - xc(2);
	xp = x2D_seen(:,1) - xc_seen(1);
	yp = x2D_seen(:,2) - xc_seen(2);
	r2 = xp.^2 + yp.^2;
	rho = sum(sqrt(x.^2+y.^2).*sqrt(r2))/sum(r2)
	rcos = sum((xp.*x + yp.*y).*r2)/sum(r2.^2);
	rsin = sum((x.*yp - xp.*y).*r2)/sum(r2.^2);
	% rsin = sqrt(rho.^2 - rcos.^2);
	a(1) = xc(1) - ( rcos*xc_seen(1) + rsin*xc_seen(2));
	a(2) = xc(2) - (-rsin*xc_seen(1) + rcos*xc_seen(2));
	a(3) = rcos;
	a(4) = rsin;
	a
	corrPar(:,icam) = a';
	
	% % fit only translation
	% A = zeros(2*n, 2);
	% b = [x2D(:,1) - x2D_seen(:,1); x2D(:,2) - x2D_seen(:,2)];
	% A = [[ones(n,1) zeros(n,1)]; [zeros(n,1) ones(n,1)]];
	% [U s V] = svd(A,0);
	% x = U'*b;
	% minsv = 1E-10;
	% a = zeros(2,1);
	% for i=1:2
	%	a = a + x(i)/s(i,i)*V(:,i);
	% end
	% a = [a; 1; 0;]
	% corrPar(:,icam) = a;
	
% 	figure;
% 	plot(x2D(:,1), x2D(:,2), 'bx');
% 	hold on;
% 	plot(x2D_seen(:,1), x2D_seen(:,2), 'r.');
	x2D_corr = [a(3)*x2D_seen(:,1)+a(4)*x2D_seen(:,2)+a(1), -a(4)*x2D_seen(:,1)+a(3)*x2D_seen(:,2)+a(2)];
	xerr = sqrt(mean((x2D(:,1)-x2D_corr(:,1)).^2));
	yerr = sqrt(mean((x2D(:,2)-x2D_corr(:,2)).^2));
	[xerr yerr]
% 	plot(x2D_corr(:,1), x2D_corr(:,2), 'ro');
% 	for i = 1:n
% 		plot([x2D_seen(i,1) x2D_corr(i,1)], [x2D_seen(i,2) x2D_corr(i,2)], 'r-');
% 	end
% 	title(strcat('camera ', num2str(icam)));
% 	xlabel('x (pixel)');
% 	ylabel('y (pixel)');
end
figure(nf_PDF);
subplot(221); xlabel('{\Delta}x'); ylabel('P({\Delta}x)');
subplot(222); xlabel('{\Delta}y'); ylabel('P({\Delta}y)');
subplot(223); xlabel('{\Delta}z'); ylabel('P({\Delta}z)');
subplot(224); xlabel('h'); ylabel('P(h)');

