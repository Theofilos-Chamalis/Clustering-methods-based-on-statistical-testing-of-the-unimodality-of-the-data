function [X,C,Types] = Varclusterstructures(n,k,d,varargin)
%Varclusterstructures - generator of a c-separated mixture of various 
%                       cluster structures (shapes)
%
%[X,C] = Varclusterstructures(n,k,d,varargin)
%  n - sample size
%  k - number of components
%  d - dimension
%
%  Optional:
%  balance - the balance of the produced clusters: 1/(balance*k) 
%            (e.g. balance=2 means that a cluster may be at least half the balanced cluster size 1/k.
%            if not defined, perfect balance is considered (i.e. balance=1).
%  csep - separation threshold, if not defined csep = sqrt(d) is set
%  info - a cell with data structs having the parameters for each cluster type
%  probtype - the probability to generate a cluster of each different cluster structure
%            (1) gaussian (2) student t (3) uniform rectangle (4) uniform oval
%  rotate - [probRot num_min num_max theta_min theta_max] each cluster will be rotated with probability probRot.
%         If it is to be rotated, then it will applied "num in num_min...num_max" times by selecting each time a random 
%         dimension to rotate. The theta of rotation will be selected from range [theta_min, theta_max].
%  rndseed - specific random seed
%returns
%  X - data set (n x d);
%  C - the instance labels
%  Types - type of clusters generated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argyris Kalogeratos, April 2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GAUSSIAN = 1; 
STUDENT_T = 2; 
UNIFORM_RECTANGLE = 3; 
UNIFORM_OVAL = 4;

type_num = 4;
Types = zeros(type_num,1);


%-- parse input arguments --%
[found, balance, varargin] = Parsepar(varargin, 'balance');
if (~found), balance = 1; end

[found, csep, varargin] = Parsepar(varargin, 'csep');
if (~found), csep = sqrt(d); end

[do_rotate, rotate, varargin] = Parsepar(varargin, 'rotate');
if (isempty(rotate)), do_rotate = 0; end
if (~do_rotate), theta = 0; end

[found, probType, varargin] = Parsepar(varargin, 'probType');
if (isempty(probType)), probType = ones(type_num,1) ./ type_num; end

[found, rndseed] = Parsepar(varargin, 'rndseed');
if (isempty(rndseed)), rndseed = sum(100*clock); end
rand('state', rndseed);  randn('state', rndseed);


rot = struct;
if (do_rotate)
     rot.prob = rotate(1);      rot.num_min = rotate(2);   rot.num_max = rotate(3);
     rot.theta_min = rotate(4); rot.theta_max = rotate(5);
else rot.prob = 0;
end

if (sum(probType) ~= 1)
    probType = cumsum(probType ./ sum(probType));
    disp('Warning: normalizing probability vector for cluster types!');
else
    probType = cumsum(probType);
end


C = zeros(n,1);

% mixing weights
if (balance == 1)
    W = ones(k,1) ./ k;
else
    while (1)
        W = rand(k,1); 
        W = W / sum(W);
        if all(W > 1/(balance*k)), break; end
    end
end

% create well-separated clusters
while (1)
    M = randn(k,d)*k;
    if (k > 1),
        D = sqrt(Sqdist(M',M'));
        if min(D(find(D))) >= csep*sqrt(d), break; end
    else
        break
    end
end

% create clusters
X = []; N = 0;
for j = 1:k
    nj = ceil(n*W(j));

    type = find(rand <= probType, 1);
    Types(type) = Types(type) + 1;

    switch (type)
    case GAUSSIAN
        fprintf('Gaussian, \n');
        U = rand(d,d)-0.5; 
        U = sqrtm(inv(U*U')) * U;
        U = diag(rand(1,d)+0.2) * U;
        XCluster = randn(nj,d) * U;

    case STUDENT_T
        min_df = 1; max_df = 3;
        [XCluster, df] = Student_t_cluster (nj, d, [], min_df, max_df);
        fprintf('Student t, '); disp(df);

    case UNIFORM_RECTANGLE
        scale = 0.5 + rand*2;
        dimvar = scale*(1+(2*rand(d,1)));
        dimvar = dimvar - min(dimvar);
        dimvar(dimvar==0) = 1;

        rotate = set_rotation (rot, d);
        XCluster = Ellipse_circle_rectangle(nj, d, 'rectangle', 'side_ratio', 0, 'density_factor', 1, 'dim_var', dimvar, 'rotate', rotate);
        fprintf('Uniform rectangle, '); disp(dimvar');

    case UNIFORM_OVAL
        scale = 0.5 + rand*2;
        dimvar = scale*(1+(2*rand(d,1)));
        dimvar = dimvar - min(dimvar);
        dimvar(dimvar==0) = 1;

        rotate = set_rotation (rot, d);
        XCluster = Ellipse_circle_rectangle(nj, d, 'ball', 'side_ratio', 0, 'density_factor', 1, 'dim_var', dimvar, 'rotate', rotate);
        fprintf('Uniform oval, '); disp(dimvar');
    end
    
    XCluster = repmat(M(j,:),nj,1) + XCluster;
    X = [X; XCluster];
    C(N+1 : N+nj) = j;
    N = N+nj;
end

% standardize to zero mean and unit variance
n = size(X,1);
X = X - repmat(mean(X),n,1);
X = X ./ repmat(std(X),n,1);

end


function rotate = set_rotation (rot, d)
    rotate = [];
    
    if (d > 1)
        if (rand < rot.prob) % do rotate
            if (rot.num_max == rot.num_min),
                 num = rot.num_min;
            else num = rot.num_min - 1 + randsample(rot.num_max - rot.num_min, 1);
            end

            for i=1:num,
                dperm = randperm(d);
                rotate = [rotate  dperm(1:2)  rot.theta_min + rand*(rot.theta_max - rot.theta_min)];
            end
        end
    end
end

