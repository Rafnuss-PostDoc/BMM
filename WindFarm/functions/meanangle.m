function [out] = meanangle(in,dim,sens)
% MEANANGLE will calculate the mean of a set of angles (in radian) based
% on polar considerations.
%
% Usage: [out] = meanangle(in,dim)
%
% in is a vector or matrix of angles (in radian)
% out is the mean of these angles along the dimension dim
%
% If dim is not specified, the first non-singleton dimension is used.
%
% A sensitivity factor is used to determine oppositeness, and is how close
% the mean of the complex representations of the angles can be to zero
% before being called zero.  For nearly all cases, this parameter is fine
% at its default (1e-12), but it can be readjusted as a third parameter if
% necessary:
%
% [out] = meanangle(in,dim,sensitivity)
%
% Written by J.A. Dunne, 10-20-05
%
if nargin<3
    sens = 1e-12;
end
if nargin<2
    ind = find(size(in)>1, 1 );
    if isempty(ind)
        %This is a scalar
        out = in;
        return
    end
    dim = ind;
end
if in>2*pi
    warning('Convert input from degree to radian')
    in = deg2rad(in);
end
mid = nanmean(exp(1i*in),dim);
out = atan2(imag(mid),real(mid));
out(abs(mid)<sens) = nan;