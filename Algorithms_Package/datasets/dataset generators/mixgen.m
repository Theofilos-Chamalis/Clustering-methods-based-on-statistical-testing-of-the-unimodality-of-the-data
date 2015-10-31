function X = mixgen(n,d,k) %[X,logl] = mixgen(n,d,k)
%mixgen(n,d,k) - generate n d-dim points forming k clusters

% Nikos Vlassis, 2000

rand('state',sum(100*clock));

while 1
  W = rand(k,1); 
  W = W / sum(W);
% if all(W > 0.1)
  if all(W > 0.03)
    break;
  end
end

M = 4 * (rand(k,d) - .5);
R = zeros(k,d^2);
X = [];
for j = 1:k
  nj = ceil(n*W(j));
  U = rand(d,d)-0.5; 
  U = sqrtm(inv(U*U'))*U;
  Xj = randn(nj,d) * diag(0.1*rand(1,d)+0.02) * U;
  X = [X; repmat(M(j,:),nj,1) + Xj];
  Rj = chol(cov(Xj,1));
  R(j,:) = Rj(:)';
end

%logl = mean(log(em_gauss(X,M,R)*W));

% standarization of data
X = X - repmat(mean(X), size(X, 1), 1);
X = X ./ repmat(std(X), size(X, 1), 1);
