function X = spiral(n, d, rounds, distortion_mode, sigma) 

% rounds = 10;
% n = 1000;
% d = 2;
% sigma = 0.2;
% distortion_mode = 0;

% just plot a spiral line
%t = linspace(0, rounds*pi, line_steps);
%x = t .* cos(t);
%y = t .* sin(t);

s = rand(n,1) * rounds*pi;
X = [s .* cos(s)    s .*sin(s)];

% standardize data
meanX = sum(X,1) / size(X,1);
stdX = std(X);
X = (X - repmat(meanX, n,1)) ./ repmat(stdX,n,1);

if (distortion_mode == 1) % uniform distortion
    X = X + sigma * randn(n,d);
else
    X = X + sigma * rand(n,d);
end
plot(X(:,1),X(:,2), '.')