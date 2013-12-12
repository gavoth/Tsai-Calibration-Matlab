function [pAB, pBA, tprojAB, tz, tpairAB] = gv_pairs_newpairs(camSpecA, posA, camSpecB, posB, mindist)
% search for pairs of ray of sight that intersect each other

% project rays from camera A onto the image plane of camera B
t1 = cputime;

% first, position of the perspective center of camera A on the image plane
% of camera B
epszv = 1.E-10;
XcvA_B = ((camSpecB.T) * ([camSpecA.xc camSpecA.yc camSpecA.zc 1]'))';
% Note that it is possible for the projection of the center of camera A to
% lie at infinity on the image plane of camera B
if abs(XcvA_B(3)) > epszv 
    XcpA_B = XcvA_B(1:2) / XcvA_B(3) * (camSpecB.d_imgplane);
else
    if XcvA_B(3) > 0
        XcpA_B = XcvA_B(1:2) * (2.E+12) * (camSpecB.d_imgplane);
    else
        XcpA_B = XcvA_B(1:2) * (-2.E+12) * (camSpecB.d_imgplane);
    end
end
% distance between camera center images
dCC = sqrt(sum(XcpA_B.*XcpA_B));
% unit vector pointing from camera A center to camera B center
uCC = -XcpA_B / dCC;
% determine the line (on the image plane) that is normal to the line 
% connecting both camera centers
lnorm = [XcpA_B(2) -XcpA_B(1)] / dCC;

% the transform materix that transfer particle images on camera A to the
% view coordinates of camera B
T_AB = (camSpecB.T) * (camSpecA.Tinv);

% transfer and project particle images of camera A onto the image plane of
% camera B
NA = length(posA);
XppA_B = zeros(NA,2);
XpvA = zeros(NA,2);
% convert units, coordinates given in pos are in pixels
XpvA(:,1) = (posA(:,1) / (camSpecA.Npix) -0.5) * camSpecA.wCCD;     % image position, in cm
XpvA(:,2) = (posA(:,2) / (camSpecA.Npix) -0.5) * camSpecA.hCCD;
for i=1:NA
    XpvA_B = (T_AB * ([XpvA(i,:) camSpecA.d_imgplane 1]'))';
    if abs(XpvA_B(3)) > epszv 
        XppA_B(i,:) = XpvA_B(1:2) / XpvA_B(3) * (camSpecB.d_imgplane);
    else
        if XpvA_B(3) > 0
            XppA_B(i,:) = XpvA_B(1:2) * (1.E+10) * (camSpecB.d_imgplane);
        else
            XppA_B(i,:) = XpvA_B(1:2) * (-1.E+10) * (camSpecB.d_imgplane);
        end
    end
end

t2 = cputime;
tprojAB = t2 - t1;

% project particle images of A onto the line lnorm, from the center of
% camera A. lppA is in units of cm
for i = 1:NA
    vtemp = XppA_B(i,:) - XcpA_B(1:2);
    lppA(i) = sum(vtemp.*lnorm) / sum(vtemp.*uCC) * dCC;
end

% project particle images of B onto the line lnorm, from the center of
% camera A. lppB is in units of cm
NB = length(posB);
for i = 1:NB
    XppB = (posB(i,1:2) / (camSpecB.Npix) - 0.5)* (camSpecB.hCCD);
    vtemp = XppB - XcpA_B(1:2);
    lppB(i) = sum(vtemp.*lnorm) / sum(vtemp.*uCC) * dCC;
end

t3 = cputime;
tz = t3-t2;

% finally, pairing particle images that have nearly the same coordniates on
% lnorm
%
% method 1: a brute force algorithm, which is O(N^2)
% pAB = zeros(NA,1);
% pBA = zeros(NB,1);
% lmin = mindist / (camSpecA.Npix) * camSpecA.hCCD;
% npairB = zeros(NB,1);
% for i = 1:NA
%     lAmin = lppA(i) - lmin;
%     lAmax = lppA(i) + lmin;
%     npairA = 0;
%     for j = 1:NB
%         if lppB(j) >= lAmin & lppB(j) <= lAmax
%             npairA = npairA +1;
%             pAB(i,npairA) = j;
%             npairB(j) = npairB(j) + 1;
%             pBA(j,npairB(j)) = i;
%         end
%     end
% end

% % method 2: use a serach grid with 1/2 pixel, which is O(N)
% % first indexing lppB in a grid with size lgrid
pAB = zeros(NA,1);
pBA = zeros(NB,1);
npairB = zeros(NB,1);
lBmax = max(lppB);
lBmin = min(lppB);
dgrid = 1/2;                                    % grid size in pixel
lgrid = dgrid / (camSpecB.Npix) * camSpecB.hCCD;    % grid size in cm
Ngrid = floor( (lBmax-lBmin)/lgrid ) + 1;
nBind = zeros(Ngrid,1);
for i=1:NB
    jB = floor((lppB(i)-lBmin)/lgrid) + 1;
%     if isnan(jB) | isinf(jB)
%         X = [i lBmin lppB(i) lgrid posB(i, 1:2)]
%     else
        nBind(jB) = nBind(jB) + 1;
        indB(jB, nBind(jB)) = i;
%     end
end
% search through lppA for near neighbors
Nsearch = ceil(mindist/dgrid);
lsearch = mindist / (camSpecB.Npix) * camSpecB.hCCD;
NBind = size(indB,2);
for i=1:NA
    npairA = 0;
    jA = floor((lppA(i)-lBmin)/lgrid) + 1;
    jmin = max(1, jA - Nsearch);
    jmax = min(Ngrid, jA + Nsearch);
    lAmin = lppA(i) - lsearch;
    lAmax = lppA(i) + lsearch;
    for j = jmin:jmax
        for nB = 1:NBind
            iB = indB(j,nB);
            if (iB > 0) & (lppB(iB) >= lAmin) & (lppB(iB) <= lAmax)
                npairA = npairA +1;
                pAB(i,npairA) = iB;
                npairB(iB) = npairB(iB) + 1;
                pBA(iB,npairB(iB)) = i;
            end
        end    
    end     % end for j=jmin:jmax
end     % end for i=1:NA% pAB = zeros(NA,1);


% tpairAB = etime(clock, t3);
tpairAB = cputime - t3;