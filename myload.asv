function A = myload(fname, commentsymbol)
% This function is an extension of matlab's load() function.
% It can accept any symbols as comment indicator (not necessarily '%'),
% e.g, myload(fname, '#') can read files compatible with gnuplot.
%
% Syntax:
%	X = myload(fname);
%	X = myload(fname, commentsymbol);
%
% If commentsymbol is missing, then the built-in matlab function load() is
% called, and the comment is indicated by '%'.
%
if (nargin > 2) | (nargin < 1)
	error('myload takes 1 or 2 input parameters. Type help myload for details.');
end

if nargin == 1
	commentsymbol = '%';
end
if ~ischar(commentsymbol) | length(commentsymbol) > 1
    error('Input variable commentsymbol must be char, NOT string');
end

if commentsymbol == '%'     % use matlab's own load()
    A = load(fname);
else
    fid = fopen(fname, 'r');
    
    nextline = fgetl(fid);
    nrow = 1;
    while ischar(nextline)
        row = [];
        [str nextline] = strtok(nextline);
        while ~isempty(str) & str(1) ~= commentsymbol
            row = [row str2num(str)];
            [str nextline] = strtok(nextline);
        end
        if ~isempty(row)
            if nrow == 1
                A = row;
                nAcol = size(row,2);
            else
                if nAcol == size(row,2)
                    A = [A; row];
                else
                    error('Data file must contain the same number of columns.');
                end
            end
            nrow = nrow + 1;
        end
        nextline = fgetl(fid);
    end
    fclose(fid);
end
