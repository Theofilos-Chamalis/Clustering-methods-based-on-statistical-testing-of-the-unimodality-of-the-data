function d = sqdist(a,b)
% sqdist - computes pairwise squared Euclidean distances between points

% original version by Roland Bunschoten, 1999

aa = sum(a.*a); bb = sum(b.*b); ab = a'*b; 
d = abs(repmat(aa',[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
