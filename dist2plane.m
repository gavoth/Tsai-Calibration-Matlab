function [h, varargout] = dist2plane(O, A, B, P)
% Function DIST2PLANE calculates the distance from point P to a plane
% defined by 3 points O, A, and B. It returns the distance and optionally
% the projection of P on the plane.
%
% syntax:
%   h = dist2plane(O, A, B, P);
%   [h, Q] = dist2plane(O, A, B, P);
%
% Inputs:
%   O, A, B --  3 points on the plane
%   P       --  The point of interest
%
% Outputs:
%   h   --  distance from P to plane OAB
%   Q   --  projection of P on OAB
%
if nargin ~= 4
    error('Function DIST2PLANE takes 4 arguments. Type help for details.');
end
if nargout > 2
    error('Function DIST2PLANE gives at most 2 outputs.');
end

% find the unit normal of plane OAB
l1 = A - O;
l2 = B - O;
l = zeros(3,1);
l(1) = l1(2)*l2(3) - l1(3)*l2(2);
l(2) = l1(3)*l2(1) - l1(1)*l2(3);
l(3) = l1(1)*l2(2) - l1(2)*l2(1);
n = l/sqrt(sum(l.*l));

% calculate distance and projection
h = abs(sum((P-O).*n));
if nargout == 2
    varargout{1} = P - sum((P - O).*n)*n;
end