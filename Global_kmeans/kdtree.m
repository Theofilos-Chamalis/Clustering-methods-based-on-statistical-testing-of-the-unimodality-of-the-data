function C = kdtree(X,I,C,b)
%kdtree - recursive building of a k-dimensional tree
%
%C = kdtree(X,I,C,b)
%  X - (n x d) matrix of input data 
%  I - (n x 1) indeces to points X 
%  C - (m x d) dynamic list containing the centroids of the terminal nodes
%  b - bucket size
%returns
%  C - (m x d) updated C

% Nikos Vlassis, 2001, http://www.science.uva.nl/~vlassis 
% based on (Sproull, Algorithmica 6, 1991)

XI = X(I,:);
ni = size(XI,1);

if ni==1
  C = [C; XI];

elseif (ni>0) & (ni <= b) 		% terminal node
  
  C = [C; mean(XI)];   

elseif ni>0 

  % project XI on their 1st pc 
  M = mean(XI);
  S = cov(XI);
  [W,var] = eig(S);[var,m]=max(diag(var));
  W=W(:,m);
  XIproj = (XI - ones(ni,1)*M) * W;

  % split projection dimension perpendicularly on the mean
  cutval = 0;
  Ileft = I(find(XIproj <= cutval));
  Iright = I(find(XIproj > cutval));

  % recursively call kd-tree on subpartitions
  C = kdtree(X,Ileft,C,b);
  C = kdtree(X,Iright,C,b);

end
