%------------
% Initialization for Clustering Algorithms based on Unimodality Tests package previously
% called Dip-means package. Run bistest.m for a demonstration.
%------------
% Copyright (C) 2012-2015, Chamalis Theofilos, Kalogeratos Argyris.
%------------

clear;
fclose('all');
clc;

% copyright/copyleft info
fprintf('Clustering Algorithms based on Unimodality Tests package v.0.2. Copyright (C) 2012-2015 Chamalis Theofilos, Kalogeratos Argyris.\n'); 
fprintf('This is free software distributed under the GNU General Public License; for details see LICENSE.txt.\n\n');

% add folders and their subfolders to the path
addpath(genpath('Algorithms_Package'));
addpath(genpath('Experiment_report'));
addpath(genpath('Global_kmeans'));
addpath(genpath('Plots2D_Of_Datasets'));
addpath(genpath('PG_means'));
addpath(genpath('External_Bootstrap'));
addpath(genpath('wgPlot'));
            
%%% open important files [optional]
% files to run frameworks' algorithms
edit('bistest.m');

% main files that call the clustering algorithms
edit('bisect_kmeans_default.m');
edit('bisect_agglodip.m');
edit('bisect_agglopdip.m');

% files for Hartigan's test
%edit('HartigansDipSignifTest_no_boot.m');
%edit('HartigansDipSignifTest.m');
%edit('HartigansDipTest.m');
%edit('hartigansdiptestdemo.m');

% some of the functions used
%edit('test_projected_unimodality.m');
%edit('test_unimodal_cluster.m');
%edit('test_unimodal_cluster2.m');
%edit('data_Proj');
%edit('MergeSmallClusters.m');
%edit('UnimodalityChecking.m');
%edit('plotClusterResults');
%edit('compute_unifpdfboot');