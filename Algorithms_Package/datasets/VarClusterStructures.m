function [X,Cx] = VarClusterStructures(n,k,d,c)
%VarClusterStructures - generator of a c-separated mixture of various 
%                       cluster structures (shapes)
%
%[X,T] = VarClusterStructures(n,k,d,c) 
%  n - sample size
%  k - number of components
%  d - dimension
%  c - separation threshold
%returns
%  X - data set (n x d);


% Argyris Kalogeratos, April 2012.

rand('state',sum(100*clock));
Cx = zeros(n,1);

% mixing weights
while 1
  W = rand(k,1); 
  W = W / sum(W);
  if all(W > 1/(2*k))
    break;
  end
end

% create well-separated clusters
while (1)
    M = randn(k,d)*k;
    if (k > 1), 
        D = sqrt(sqdist(M',M'));
        if min(D(find(D))) >= c*sqrt(d)
            break;
        end
    end
end

% create Gaussian clusters, training X and test T set
X = [];
T = [];
Nx = 0; Nt = 0;
for j = 1:k
  nj = ceil(n*W(j));
  U = rand(d,d)-0.5; 
  U = sqrtm(inv(U*U')) * U;
  U = diag(rand(1,d)+0.2) * U;
  Xj = randn(nj,d) * U;
  X = [X; repmat(M(j,:),nj,1) + Xj];
  Cx(Nx+1 : Nx+nj) = j;
  Nx = Nx+nj;
end

% standardize to zero mean and unit variance
n = size(X,1);
X = X - repmat(mean(X),n,1);
X = X ./ repmat(std(X),n,1);
plot(X(:,1),X(:,2), '.')