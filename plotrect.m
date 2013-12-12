function h = plotrect(rectarea, varargin)
% This function plot a rectangle in current figure.
%
% Syntax:
%   plotrect(rectarea);
%   plotrect(rectarea, Linespec, ...);
%   plotrect(rectarea, 'PropertyName', PropertyValue, ...);
%   h = plotrect(...);
%
% Inputs:
%   rectarea    --  [xmin, xmax, ymin, ymax]
%
% Output:
%   h   --  a column vector of 4 handles corresponding to lines consisting
%           of the rectangle, starts from the left-vertical one and goes
%           counterclockwise
%
hold on;
xmin = rectarea(1);
xmax = rectarea(2);
ymin = rectarea(3);
ymax = rectarea(4);
if nargin == 1
plot([xmin xmin xmax xmax xmin], [ymax ymin ymin ymax ymax]);
else
plot([xmin xmin xmax xmax xmin], [ymax ymin ymin ymax ymax], char(varargin{1:nargin-1}));
end
hold off;