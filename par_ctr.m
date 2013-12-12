function [xc, yc, Ap, varargout] = par_ctr(Iimg, threshold, Amin, ctrfind_method, flagshow, area)
% This script finds particle centers and particle size
%
% syntax:
%	[xc, yc, Ap] = par_ctr(Iimg, threshold, Amin, ctrfind_method, flagshow, area);
%

if( (nargin < 5) | (nargin > 6) )
	error('par_ctr.m takes 5 or 6 input parameters. Type help for details');
end

[npixy, npixx] = size(Iimg);
%Ibright = uint8(ones(npixy, npixx)*255);
%Ainv = imsubtract(Ibright, Iimg);

% only process the area indicated
if(nargin == 5)
	% Ainv = uint8(Iimg);
	Ainv = Iimg;
else
	Ainv = zeros(npixy, npixx);
	Ainv(area(3):area(4), area(1):area(2)) = Iimg(area(3):area(4), area(1):area(2));
%	for j=1:npixy
%	for i=1:npixx
%		if( (i>=area(1)) & (i<=area(2)) & (j>=area(3)) & (j<=area(4)) )
%			Ainv(j,i) = uint8(Iimg(j,i));
%			Ainv(j,i) = Iimg(j,i);
%		else
%			Ainv(j,i) = uint8(0);
%		end
%	end
%	end
end

B = double(Ainv)>threshold;

if strcmp(lower(flagshow), 'show') == 1
	figure;
	subplot(211)
	imshow(uint8(Iimg));
	if(nargin == 6)
		hold on;
		plotrect(area, 'c-');
	end
	subplot(212)
	imshow(B);
	title('thresholded image');
end

if strcmp(lower(ctrfind_method), 'com') == 1
	[pseg, Nseg] = seqlabel(Ainv, threshold, Amin);
	npar = length(Nseg);
	par = zeros(npar,3);
	% find particle size and center of mass
	nstart = 0;
	for j = 1:npar
		pamaskPos_x.datr(j,3) = Nseg(j);
		ind = [nstart+1:nstart+Nseg(j)]';
		weight = sum(pseg(ind,3));
		par(j,1) = sum(pseg(ind,2).*pseg(ind,3))/weight;
		par(j,2) = sum(pseg(ind,1).*pseg(ind,3))/weight;
		nstart = nstart+Nseg(j);
	end
	Ap = par(:,3);
end
if strcmp(lower(ctrfind_method), 'gaussianfit') == 1
	myeps = 1.E-10;
	[pIseg, pCseg] = seqlabel(Ainv, threshold, 'center', Amin);
	pIseg = max(double(pIseg), myeps);
	[npar, dummy] = size(pCseg);
	par = zeros(npar,2);
	% particle size
	Ap = pCseg(:,3);
	% find particle center by gaussian fit
	% note in the output from seqlabel_par, pIseg is in (j, i) format,
	% but pCseg is [i, j]
	for j = 1:npar
		Ip = pIseg(:,:,j);
		par(j,2) = pCseg(j, 1) + 0.5*log(Ip(2,3)/Ip(2,1))/(log((Ip(2,2)*Ip(2,2))/(Ip(2,1)*Ip(2,3))));
		par(j,1) = pCseg(j, 2) + 0.5*log(Ip(3,2)/Ip(1,2))/(log((Ip(2,2)*Ip(2,2))/(Ip(1,2)*Ip(3,2))));
	end
end

if strcmp(lower(flagshow), 'show') == 1
		% plot the particle centers
		subplot(211)
		hold on
		for j = 1:npar
			plot(par(j,2), par(j,1), 'r+');
			%text(par(j,2), par(j,1), num2str(j), 'Color', 'b');
		end
		hold off
		subplot(212)
		hold on
		for j = 1:npar
			plot(par(j,2), par(j,1), 'r+');
			%text(par(j,2), par(j,1), num2str(j), 'Color', 'b');
		end
		hold off
end
xc = par(:,2);
yc = par(:,1);

if nargout == 4
	varargout(1) = {B};
end
