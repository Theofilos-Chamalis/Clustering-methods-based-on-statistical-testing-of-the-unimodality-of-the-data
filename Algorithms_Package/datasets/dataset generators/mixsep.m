function [X,T] = mixsep(n,k,d,c)
%mixsep - generator of a c-separated mixture
%
%[X,T] = mixsep(n,k,d,c) 
%  n - sample size
%  k - number of components
%  d - dimension
%  c - separation threshold
%returns
%  X - training set (n x d); samples from the mixture
%  T - test set (n/2 x d) 

% Nikos Vlassis, 2000

rand('state',sum(100*clock));

% mixing weights
while 1
  W = rand(k,1); 
  W = W / sum(W);
  if all(W > 1/(2*k))
    break;
  end
end

% create well-separated clusters
while 1
  M = randn(k,d)*k;
  D = sqrt(sqdist(M',M'));
  if min(D(find(D))) >= c*sqrt(d)
    break;
  end
end

% create Gaussian clusters, training X and test T set
X = [];
T = [];
for j = 1:k
  nj = ceil(n*W(j));
  U = rand(d,d)-0.5; 
  U = sqrtm(inv(U*U')) * U;
  U = diag(rand(1,d)+0.2) * U;
  Xj = randn(nj,d) * U;
  X = [X; repmat(M(j,:),nj,1) + Xj];
  Tj = randn(ceil(nj/2),d) * U;
  T = [T; repmat(M(j,:),ceil(nj/2),1) + Tj];
end

% standardize to zero mean and unit variance
n = size(X,1);
m = size(T,1);
meanX = mean(X);
X = X - repmat(meanX,n,1);
T = T - repmat(meanX,m,1);
stdX = std(X);
X = X ./ repmat(stdX,n,1);
T = T ./ repmat(stdX,m,1);
