function [features labels] = pinwheel( radial_std, tangential_std, ...
                                       num_classes, num_per_class, rate)
%
% [features labels] = PINWHEEL( radial_std, tangential_std, num_classes,
%                               num_per_class, rate )
% 
% This function generates a "pinwheel" data set.  It has as many arms as
% classes.  It generates them by taking Gaussian distributions,
% stretching them and then rotating them appropriately.  The centers are
% equidistant around the unit circle.
%
% INPUT:
%   - radial_std:     the standard deviation in the radial direction
%   - tangential_std: the standard deviation in the tangential direction
%   - num_classes:    how many arms and classes to generate
%   - num_per_class:  how many of each class to generate
%   - rate:           how many radians to turn per exp(radius)
%
% OUTPUT:
%   - features: the 2d locations in space
%   - labels:   the actual class labels
%
% Reasonable usage example:
%  >> X = pinwheel(0.3, 0.3, 3, 1000, 0.25);
%  >> plot(X(:,1), X(:,2), '.');
%
% Copyright: Ryan Prescott Adams, 2008
% This is released under the GNU Public License.
% http://www.gnu.org/licenses/gpl-2.0.txt
%

  
  % Find the equidistant angles.
  rads = linspace(0, 2*pi, num_classes+1);
  rads = rads(1:end-1);
  
  features = randn([ num_classes*num_per_class 2 ]) ...
      * diag([tangential_std radial_std]) ...
      + repmat([1 0], [num_classes*num_per_class 1]);
  labels   = vec(repmat([1:num_classes], [num_per_class 1]) ...
                 .* ones([num_per_class num_classes]));
  angles   = rads(labels)' + rate*exp(features(:,1));

  rowSize = size(angles);
  %for i=1:rows(angles) changed that bit
  
% This would probably be faster if vectorized.
  for i=1:rowSize(1)
    features(i,:) = features(i,:) ...
        * [ cos(angles(i)) -sin(angles(i)) ; ...
            sin(angles(i))  cos(angles(i))];
  end
end

% Y = VEC(x)  Given an m x n matrix x, this produces the vector Y of length
%   m*n that contains the columns of the matrix x, stacked below each other.
%
% See also mat.

function x = vec(X)

%
% This file is part of SeDuMi 1.1 by Imre Polik and Oleksandr Romanko
% Copyright (C) 2005 McMaster University, Hamilton, CANADA  (since 1.1)
%
% Copyright (C) 2001 Jos F. Sturm (up to 1.05R5)
%   Dept. Econometrics & O.R., Tilburg University, the Netherlands.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% Affiliation SeDuMi 1.03 and 1.04Beta (2000):
%   Dept. Quantitative Economics, Maastricht University, the Netherlands.
%
% Affiliations up to SeDuMi 1.02 (AUG1998):
%   CRL, McMaster University, Canada.
%   Supported by the Netherlands Organization for Scientific Research (NWO).
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA
% 02110-1301, USA
%

[m n] = size(X);
x = reshape(X,m*n,1);
end
                        
