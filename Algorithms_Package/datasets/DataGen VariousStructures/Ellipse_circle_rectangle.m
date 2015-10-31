% X = ellipse_circle_rectangle(n, d, stype, varargin)
% This function creates a group of n d-dimensional points with certain
% structure in space.
% Arguments:
%   - stype: 'ball' for spherical structures, and 'rectangle' for parallilograms
%
% Optional arguments:
%   - side_ratio: for a ball it is the ratio of the inner radius to the outer (0<side_ratio<=1).
%            It controls the slice of the shell of the ball in which the data will be produced. 
%               side_ratio = 1: points on n-ball's surface
%               side_ratio = 0: concrete n-ball
%            For a rectagle it creates a weird effect (probably not usefull, so set side_ratio=0).
%   - density_factor: it controlls the distribution of points in structure's area.
%           For a ball 
%               density_factor = 0: perfect uniform circle
%               density_factor < 1: denser in periphery
%               density_factor > 1: denser in center area
%           For a rectangle
%               density_factor = 0: perfect uniform rectangle
%               density_factor < 1: denser FAR from the corner of smallest coordinates
%               density_factor > 1: denser CLOSE to the corner of smallest coordinates
%   - dim_var:  the matrix containing the d scaling factors for each dimension.
%               This is usefull to create uniform ovals or rectangles with arbitary sides' ratio.
%               If dim_var has more elements than d, then only the first d are considered.
%               If dim_var has less elements than d, then it is repeated to cover d. 
%   - rotate:   it is a matrix that contains tripplets of the form
%                    [dim1 dim2 theta1, dim3 dim4 theta ...]
%               where the first two are the dims wrt which the rotation of theta will be done.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argyris Kalogeratos, April 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function X = ellipse_circle_rectangle(n, d, stype, varargin)

%-- parse input arguments --%
[found, side_ratio, varargin] = parsepar(varargin, 'side_ratio');
if (~found), side_ratio = 0; end

[found, density_factor, varargin] = parsepar(varargin, 'density_factor');
if (~found), density_factor = 1; end

[found, dim_var, varargin] = parsepar(varargin, 'dim_var');
if (~found), dim_var = ones(1,d); end

[do_rotate, rotate] = parsepar(varargin, 'rotate');
if (isempty(rotate)), do_rotate = 0; end

if (do_rotate)
    if (d == 2), 
        theta = rotate(3); rot_dims = [1 2];
        RotMatrix = multidim_rotmatrix (d, rot_dims, theta); % it creates: [ cos(theta) -sin(theta); sin(theta) cos(theta)];
    else % read tripplets dim1, dim2, theta
        RotMatrix = speye(d,d);
        for i=1:numel(rotate)/3,
            rot_dims = (i-1)*3 + 1;
            theta = rotate(rot_dims+2);
            rot_dims = rotate([rot_dims, rot_dims+1]);
            
            RotMatrix = RotMatrix * multidim_rotmatrix (d, rot_dims, theta);
        end
    end
end

if (~isrow(dim_var)), dim_var = dim_var'; end;

% not all dim variances are defined
if (length(dim_var) > d)
    disp('More dim vars are defined!\n');
    dim_var = dim_var(1:d);
elseif (length(dim_var) < d) % repeat dim var to cover all dimensions
    dim_var = repmat(dim_var, 1, floor(d/numel(dim_var)));
    if (length(dim_var) ~= d)
        dim_var = [dim_var; dim_var(1:d-length(dim_var))];
    end
end

if (strcmp(stype, 'ball'))
    % Marsaglia(1972) method
        X = randn(n,d);                           % use normal distribution
        norms = sqrt(sum(X.^2,2));                % compute the radius from center
        X = bsxfun(@rdivide, X, norms);           % multiply with 1/radius to get points on surface
        r = side_ratio + (1-side_ratio) * rand(n,1).^density_factor;
        X = bsxfun(@times, X, r.^(1/d));          % multiply with u^(1/d)
elseif (strcmp(stype, 'rectangle'))
    X = side_ratio + (1-side_ratio) * rand(n,d).^density_factor;
else
    error('Unknown cluster structure!');
end

% make the center be the multidim zero vector (but not std=1)
X = bsxfun(@minus, X, sum(X,1)/size(X,1));

% apply scaling in each dimension
X = bsxfun(@times, X, dim_var);

%rotate data
if (do_rotate)
     X = X * RotMatrix;
end
