function [pseg, Nseg] = seqlabel(A, th, flagout, Amin)
% This function segments image A to regions using sequential labelling.
%
% Syntax:
%	[pseg, Nseg] = seqlabel(A, th);
%	[pseg, Nseg] = seqlabel(A, th, Amin);
%	[pseg, Nseg] = seqlabel(A, th, flagout);
%	[pseg, Nseg] = seqlabel(A, th, flagout, Amin);
%
% input:
%	A       --  image
%   	th      --  threshold
%   	flagout --  output format, flagout == 'center' means to output 3x3 regions around the maxima for gaussian fit, 
%		default is 'all', which outputs all pixels belong to the segment.
%   	Amin    --  minimum area of segments to be output, default is 1
%
% output:
%	if flagout == 'all'
%   	pseg    --  pixels of region segments, pseg has 3 columns: [i j I]
%   	Nseg    --  number of pixels in each segment. For example, pixels belong to 
%               segment i are stored in pseg(Nseq(i-1)+1:Nseg(i)) with the 
%               convention that Nseq(0) = 0
%	if flagout == 'center'
%	pseg 	--	3x3 region segments, pIseg is a 3D matrixs: pIseg(i,j,p) is
%               	the intensity of particle p at (i,j), note the particle
%               	center intensity is pIseg(2,2,p).
%   	Nseg   --  	particle center coordinates
%
Nseg = [];

if (nargin < 2) | (nargin > 4)
	error('syntax error for seqlabel.m, type help seqlabel for details.'); 
end
if nargin == 2
	flagout = 'all';
	Amin = 1;
end
if nargin == 3
	if ischar(flagout)
		Amin = 1;
	else
		Amin = flagout;
		flagout = 'all';
	end
end
if nargin == 4
	if ischar(Amin)
		temp = Amin;
		Amin = flagout;
		flagout = temp;
	end
end
flagout = lower(flagout);
if ~(strcmp(flagout, 'all') | strcmp(flagout, 'center'))
	error('wrong input parameter flagout. type help for details');
end
 
[N, M] = size(A);
Aaug = zeros(N+2, M+2);
% first pass, label 4-connected pixels in right+downward directions
label = 0;
equivbuf = [];
seg = [];
for j = 2:N+1
    for i = 2:M+1
        if A(j-1, i-1) > th   % only label pixels with intensity above threshold
            leftlabel = Aaug(j, i-1);
            upperlabel = Aaug(j-1, i);
            if leftlabel > 0        % if left pixel is already labelled, check upper pixel
                if upperlabel == 0      % if upper pixel is not labelled, set current pixel to leftlabel
                    Aaug(j, i) = leftlabel;
                else                    % if upper pixel is also labelled, check consistency
                    if leftlabel == upperlabel      % if the two labels agree, no problem
                        Aaug(j, i) = leftlabel;
                    else                            % if they don't agree, set current pixel to the smaller label
                        if leftlabel < upperlabel
                            minlabel = leftlabel;
                            maxlabel = upperlabel;
                        else
                            minlabel = upperlabel;
                            maxlabel = leftlabel;
                        end
                        Aaug(j,i) = minlabel;
                        % add minlabel into the equivalence table of maxlabel
                        neq = equivbuf(maxlabel, 1);    % first element in equivalence table saves the # of equivalent labels
                        found = 0;
                        if neq > 0
                            for ieq=2:neq+1
                                if equivbuf(maxlabel, ieq) == minlabel
                                    found = 1;
                                end
                            end
                        end
                        if found == 0
                            equivbuf(maxlabel, 1) = neq+1;
                            equivbuf(maxlabel, neq+1+1) = minlabel;
                        end
                    end     % end of leftlabel == upperlabel
                end     % end if upperlabel == 0
            else    % if left pixel is not labelled
                if upperlabel > 0   % if upper pixel is labelled, use upperlabel
                    Aaug(j,i) = upperlabel;
                else    % if upper pixel is not labelled, use a new label and set its equivalence table to empty
                    label = label + 1;
                    Aaug(j,i) = label;
                    equivbuf(label, 1) = 0;
                end
            end     % end if leftlabel > 0
        end     % end if Aaug(j,i) > 0
    end     % end i=1:M
end     % end j=1:N


% Now go over equivalence tables from largest label to smallest, rearrange
% them so that it contains only the smallest eqivalent label
for i = label:-1:1
    if equivbuf(i,1) > 1    % rearrange only those have more than one equivalent label
        neq = equivbuf(i,1);
        eqlabels = equivbuf(i, 2:neq+1);
        minlabel = min(eqlabels);   % find the minimal label
        for ieq = 1:neq     % search each label on the list
            ilabel = eqlabels(ieq);
            if ilabel ~= minlabel
                neql = equivbuf(ilabel,1);
                if neql == 0    % if the equivlance table of ilabel is empty, add minlabel in
                    equivbuf(ilabel,1:2) = [1, minlabel];
                else        % if not empty, then check whether ilabel is in the table
                    found = 0;
                    for ieql = 2:neql+1
                        if equivbuf(ilabel, ieql) == minlabel;
                            found = 1;
                        end
                    end
                    if found == 0
                        equivbuf(ilabel, 1) = neql+1;
                        equivbuf(ilabel, neql+1+1) = minlabel;
                    end
                end     % end if neql == 0
            end     % end if ilabel ~= minlabel
        end     % end for ieq = 1:neq
        % after passing the minlabel to all labels on its list, we can
        % shrink the equivalence table of current label to contain only the
        % smallest label. Because the table now contains only one label, we 
        % don't need to store its lengthe either. 
        equivbuf(i,1) = minlabel;
    else
        if equivbuf(i,1) == 1   % for those who has only one equivlance label, remove the length saved
            equivbuf(i,1) = equivbuf(i,2);
        end
    end     % end if equivbuf(i,2) > 1
end
equivbuf = equivbuf(:,1);

% Second pass of equivalence table to eliminate all indirect-reference
for i = 1:label
    eqlabel = equivbuf(i);
    % if the equivalence label itself is equivalent to another label, 
    % search down the reference list until a non-reference label is found
    if (eqlabel > 0) & (equivbuf(eqlabel) > 0)     
        while eqlabel > 0
            minlabel = eqlabel;
            eqlabel = equivbuf(eqlabel);
        end
        equivbuf(i) = minlabel;
    end
end

% Walk through the image to output segments
Numseg = 0;
ind = zeros(label, 1);  % to save index of label to segment number
for j = 1:N
    for i = 1:M
        if A(j,i) > th
            label = Aaug(j+1, i+1);
            if equivbuf(label) == 0     % if this is a stand-alone label
                if ind(label) == 0  % if this label hasn't been saved
                    Numseg = Numseg +1;
                    ind(label) = Numseg;
                    temp(Numseg, 1) = 1;  % first element of temp saves the # of pixels belong to this segment
                    temp(Numseg, 2) = (j-1)*M + i;    % save the 1D order of the pixel
                else
                    iseg = ind(label);
                    temp(iseg, 1) = temp(iseg,1) +1;
                    temp(iseg, temp(iseg,1)+1) = (j-1)*M + i;
                end
            else
                eqlabel = equivbuf(label);
                iseg = ind(eqlabel);
                temp(iseg, 1) = temp(iseg, 1) + 1;
                temp(iseg, temp(iseg,1) + 1) = (j-1)*M + i;
            end
        end     % end if A(j,i) > 0
    end
end

% organize segments into the list pseg
segcount = 0;
pcount = 0;
pseg = zeros(20000,3);
for iseg = 1:Numseg
    np = temp(iseg,1);
    if np >= Amin
        segcount = segcount +1;
        Nseg(segcount) = np;
    for ip = 2:np+1
        i = mod(temp(iseg,ip), M);
        if i == 0
            i = M;
            j = floor(temp(iseg, ip)/M);
        else
            j = floor(temp(iseg,ip)/M) + 1;
        end
        pcount = pcount+1;
        pseg(pcount, 1) = i;
        pseg(pcount, 2) = j;
        pseg(pcount, 3) = A(j,i);
    end
    end
end

% Finally, if necessary, find max intensity in each region and output 3x3 pixels aroud the maximum
if strcmp(flagout, 'center');
	Np = length(Nseg);
	Aaug(2:N+1, 2:M+1) = A;
	istart = 0;
	for k=1:Np
		pix = pseg(istart+1:istart+Nseg(k),:);
		istart = istart+Nseg(k);
		Imax = pix(1,3);
		for p=1:Nseg(k)
			if pix(p,3) >= Imax
			Imax = pix(p,3);
			imax = pix(p,1);
			jmax = pix(p,2);
			end
		end
		% output 3x3 regions, if the particle is near the edge, then output
		% zeros to fill the region
		pIseg(:,:,k) = Aaug(jmax:jmax+2, imax:imax+2);
		pCseg(k,:) = [imax jmax];
	end
	pseg = pIseg;
	Nseg = [pCseg Nseg'];
end
