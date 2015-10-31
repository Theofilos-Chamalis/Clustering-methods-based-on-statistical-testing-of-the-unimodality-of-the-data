clear all;
close all;
warning off;

load X1;
%k-means
global_kmeans(X,[],5,0,0,1);
