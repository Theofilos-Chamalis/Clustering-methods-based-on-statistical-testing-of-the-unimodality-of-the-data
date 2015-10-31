
function X = ellipse_circle_rectangle(type, n, d, side_ratio, density_factor, dim_var, theta_rotation)

% type = 'rectangle';    % 'circular'    'rectangle'
% n = 1000;
% side_ratio = 0;
% density_factor = 1;
% dim_var = [1 1 1];
% theta_rotation = 0.5; % 0 ... 2*pi
% d=3;

% density_factor == 0: perfect circle
% density_factor < 1 : denser in periphery 
% density_factor < 1 : denser in center area 

if (d==2)
     RotMatrix = [ cos(theta_rotation) -sin(theta_rotation); sin(theta_rotation) cos(theta_rotation)];
end

if (~isrow(dim_var)), dim_var = dim_var'; end;

if (length(dim_var) ~= d)
    %disp('error');
    if (~isrow(dim_var)), dim_var = dim_var'; end
    dim_var = repmat(dim_var, 1, floor(d/numel(dim_var)));
    if (length(dim_var) ~= d)
        dim_var = [dim_var; dim_var(1:d-length(dim_var))];
    end
end

d_inc = false;

if (strcmp(type, 'circular'))
    if (mod(d,2) == 1), d_inc = true; d = d + 1; dim_var = [dim_var dim_var(end)]; end
    X = zeros(n,d);
    for i=1:ceil(d/2),
        theta = 2*pi*rand(n,1);
        r = sqrt(side_ratio^2 + (1-side_ratio) * rand(n,1).^density_factor);
        d1=(i-1)*2+1;  d2=d1+1;
        X(:,d1:d2) = [dim_var(d1) * r.*cos(theta), dim_var(d2) * r.*sin(theta)];
    end
elseif (strcmp(type, 'rectangle'))
    X = side_ratio + (1-side_ratio) * rand(n,d).^density_factor;
    X = bsxfun(@times, X, dim_var);
end

% make the center be the multidim zero vector
X = bsxfun(@minus, X, .5 * ones(n,1));

%rotate data
if (d==2)
     X = X * RotMatrix;
else % multidim rotation
%     Xaver = sum(X,1) / n;
end

if (d_inc == true) % remove one dimension
    X = X(:, 1:end-1);
end

%dbg:
%plot(X(:,1), X(:,2), '.');
%axis equal; axis tight;