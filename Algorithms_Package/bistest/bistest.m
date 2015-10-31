%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Clustering Algorithms based on Unimodality Tests package demonstration.
%------------
% Copyright (C) 2012-2015, Chamalis Theofilos, Kalogeratos Argyris.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear the matrices from the data of possible previous runs
% or all data in workspace
%clear ('X','C');
%clc;

% remove all warnings if any, found at older matlab versions
w = warning ('off','all');

% find all windows of type figure (if any), and close them
delete(findall(0,'Type','figure'))

% define the RNG seed
%rseed = sum(100*clock);    
rseed = 10;
% rand('state', rseed);  
% randn('state', rseed);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create or load data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%LOADING OF DATASETS
%--------------------

% Pendigits (UCI) (0,4,7 digits)
%--------------------
%load('PendigitsTestSet047');

% Synthetic Various Cluster Structures
%--------------------
%(n,k,d,[Gauss,Stud-T,UnifRect,UnifOval],separ~dflt=sqrt(d))
%[X,C] = Varclusterstructures(6000,20,32,'probtype',[0.4,0.2,0.2,0.2]);

% Iris dataset
%-------------------
%load('iris.mat'); % X and C are loaded

% Combo Setting 2d dataset (Includes 7 clusters of several distributions, ideal for demonstration)
%-------------------
load('combo_setting.mat'); % X and C are loaded    

% X10 dataset
%-------------------
%load('X10');

% uncomment the line below to hide the true clusters' members 
%clear C;


if (exist('C', 'var'))
     k = length(unique(C));
     real_k = k;
else
    real_k = -1;  % the ground truth labels are not available
end

[N,d] = size(X);  
DATASIZE = N;

% set the number of initial clusters that agglomerative algorithms
% begin with
if(real_k~=-1)
    numberOfInitialClusters = 3*real_k;
else
    numberOfInitialClusters = 9;
end

% set the minimum number of elements a cluster (agglomerative) may have
%-------------------
smallestCl = 6;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% select which methods to test (the names for each id are indicated in the next line)
methods = [1,2,3,4,5,6];

num_methods = length(methods);
method_names = {'   Agglopdip', '   Pdip-means', '   Agglodip', '   Dip-means', 'X-means', 'G-means', '',}; % the last one is kernel k-means or k-means guided by the dip-test criterion (dip-means^*)
method_name = method_names(methods);

split_struct = cell(num_methods,1);

split_struct{1} = struct;
split_struct{1}.pval_threshold    = 0.00;    % the probability <= to which the cluster must split
split_struct{1}.exhaustive_search = 1;      % whether Hartigan test must be done for all cluster objects or to stop when first prob=0 is found
split_struct{1}.voting            = 0;           % the spliting criterion is based on a voting (0<voting<=1), or on the worst indication from the objects of the cluster (voting=0)
split_struct{1}.nboot             = 1000;      % number of normal distribution samples to test each object with Hartigan test
split_struct{1}.overall_distr     = 0;          % in pdip-means we use individual object-viewers so it is set to 0

split_struct{2} = struct;
split_struct{3} = struct;
split_struct{4} = struct;
split_struct{2} = split_struct{1};
split_struct{3} = split_struct{1};
split_struct{4} = split_struct{1};

split_struct{6} = struct;
split_struct{6}.pval_threshold    = 0.999;  % the a value (alpha value)

 split_struct{7} = split_struct{3}; 

 split_trials = 10;                                   % times to try a split (we use 10 for more stable results)
 merge_trials = split_trials;

% choose whether to solve problem in kernel space (it is enabled if a non-negative kernel is provided)
kernelSpace = 0.0;    
if (kernelSpace > 0), Kernel = RBFkernel (X, kernelSpace);
                      method_names{7} = 'kernel dip-means';
else                  Kernel = [];
                      method_names{7} = 'dip-means^*';
end

result = zeros(max(methods), 7);
j = 1;

% override manually the methods that you wish to run
methods = [2];

for m=methods,
    if (m == 1)
        fprintf('\n+++++++++++++++++++++++++++++++++AGGLOPDIP++++++++++++++++++++++++++++++++\n');
        tic;
        [R, sumer, R_ref, sumer_ref] = bisect_agglopdip(X, numberOfInitialClusters, 'merge_struct', split_struct{m}, 'merge_trials', merge_trials, 'mergeSELECT', 6, 'splitMODE', 0, 'refineMODE', 1, 'smallest_cluster', smallestCl,'attempts', 1, 'rndseed', 0+rseed);
        toc;
        disp(' '); 
%         % % the next two variables can be used in all of the below methods to
%         % % compare their resulting clusters to the real ones (if we have them) 
%         R_ref_plot = R_ref;
%         k_plot = length(unique(R_ref_plot));
% 
%         %%% plot the resulting clustering found by any algorithm and the true clusters (if they exist) and pause.
%         delete(findall(0,'Type','figure'))
%         plotClusterResults;
%         pause;
    elseif (m == 2)
        fprintf('+++++++++++++++++++++++++++++++++PDIP-MEANS+++++++++++++++++++++++++++++++++++++++++\n\n');
        tic;
        [R, sumer, R_ref, sumer_ref] = bisect_kmeans_default(X, 'split_struct', split_struct{m}, 'split_trials', split_trials, 'splitSELECT', 6, 'splitMODE', 0, 'refineMODE', 1, 'smallest_cluster', smallestCl, 'attempts', 1, 'rndseed', 0+rseed);
        toc;
        disp(' ');      
%         % % the next two variables can be used in all of the below methods to
%         % % compare their resulting clusters to the real ones (if we have them) 
%         R_ref_plot = R_ref;
%         k_plot = length(unique(R_ref_plot));
% 
%         %%% plot the resulting clustering found by any algorithm and the true clusters (if they exist) and pause.
%         delete(findall(0,'Type','figure'))
%         plotClusterResults;
%         pause;
     elseif (m == 3)
        fprintf('\n+++++++++++++++++++++++++++++++++AGGLODIP+++++++++++++++++++++++++++++++++++\n');
        tic;
        [R, sumer, R_ref, sumer_ref] = bisect_agglodip(X, numberOfInitialClusters,'merge_struct', split_struct{m}, 'merge_trials', merge_trials, 'mergeSELECT', 6, 'splitMODE', 0, 'refineMODE', 0, 'smallest_cluster', smallestCl,'attempts', 1, 'rndseed', 0+rseed);
        toc;
        disp(' ');      
%         % % the next two variables can be used in all of the below methods to
%         % % compare their resulting clusters to the real ones (if we have them) 
%         R_ref_plot = R_ref;
%         k_plot = length(unique(R_ref_plot));
%
%         %%% plot the resulting clustering found by any algorithm and the true clusters (if they exist) and pause.
%         delete(findall(0,'Type','figure'))
%         plotClusterResults;
%         pause;
    elseif (m == 4)
        fprintf('+++++++++++++++++++++++++++++++++DIP-MEANS++++++++++++++++++++++++++++++++++++++++++\n\n');
        tic;
        [R, sumer, R_ref, sumer_ref] = bisect_kmeans_default(X, 'split_struct', split_struct{m}, 'split_trials', split_trials, 'splitSELECT', 3, 'splitMODE', 0, 'refineMODE', 1, 'smallest_cluster', smallestCl, 'attempts', 1, 'rndseed', 0+rseed);
        toc;
        disp(' ');
%         % % the next two variables can be used in all of the below methods to
%         % % compare their resulting clusters to the real ones (if we have them) 
%         R_ref_plot = R_ref;
%         k_plot = length(unique(R_ref_plot));
% 
%         %%% plot the resulting clustering found by any algorithm and the true clusters (if they exist) and pause.
%         delete(findall(0,'Type','figure'))
%         plotClusterResults;
%         pause;
    elseif(m == 5)
        fprintf('+++++++++++++++++++++++++++++++++X-MEANS++++++++++++++++++++++++++++++++++++++++++++\n\n');
        tic;
        [R, sumer, R_ref, sumer_ref] = bisect_kmeans_default(X, 'split_struct', split_struct{m}, 'split_trials', split_trials, 'splitSELECT', 4, 'splitMODE', 0, 'refineMODE', 1, 'smallest_cluster', smallestCl, 'attempts', 1, 'rndseed', 0+rseed);
        toc;
        disp(' '); 
%         % % the next two variables can be used in all of the below methods to
%         % % compare their resulting clusters to the real ones (if we have them) 
%         R_ref_plot = R_ref;
%         k_plot = length(unique(R_ref_plot));
% 
%         %%% plot the resulting clustering found by any algorithm and the true clusters (if they exist) and pause.
%         delete(findall(0,'Type','figure'))
%         plotClusterResults;
%         pause;
    else
        fprintf('+++++++++++++++++++++++++++++++++G-MEANS++++++++++++++++++++++++++++++++++++++++++++\n\n');
        tic;
        [R, sumer, R_ref, sumer_ref] = bisect_kmeans_default(X, 'split_struct', split_struct{m}, 'split_trials', split_trials, 'splitSELECT', 5, 'splitMODE', 0, 'refineMODE', 1, 'smallest_cluster', smallestCl, 'attempts', 1, 'rndseed', 0+rseed);
        toc;
        disp(' ');  
%         % % the next two variables can be used in all of the below methods to
%         % % compare their resulting clusters to the real ones (if we have them) 
%         R_ref_plot = R_ref;
%         k_plot = length(unique(R_ref_plot));
% 
%         %%% plot the resulting clustering found by any algorithm and the true clusters (if they exist) and pause.
%         delete(findall(0,'Type','figure'))
%         plotClusterResults;
%         pause;
    end   
     
    k = length(unique(R_ref));
    result(m, 1:3) = [k, sumer, sumer_ref];
    
    % use of clustering quality metrics
    if (real_k > 0)
     %[pq, RI, ARI, conf_matrix, conf_matrix_probC, conf_matrix_probR] =  partition_quality(C,R_ref);
        [~, RI, ARI, ~, ~, ~] =  partition_quality(C,R_ref);
        VI = varinfo(C,R_ref);
        result(m, 4:6) = [RI, ARI, VI];
    end       
    j = j+1;
end
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('----------------------------------------------------------------------------------------------------------------------------------------\n\t\t\t\tClustering results for real_k = %g\n----------------------------------------------------------------------------------------------------------------------------------------\n', real_k);
if (real_k > 0)
    for m=methods,
        fprintf('%g. %10s -- k: %3g, RI: %1.5f, ARI: %1.5f, VI: %1.5f, error: %5.5f\n', m, method_names{m}, result(m,1), result(m,4), result(m,5), result(m,6), result(m,3));
    end
else % real_k == -1: the ground truth labels are not available
    for m=methods,
        fprintf('%g. %10s -- k: %3g, error: %5.5f, (supervised measures N/A)\n', m, method_names{m}, result(m,1), result(m,3));
    end
end


% %% plot the resulting clustering found by any algorithm and the true data
% %% clusters (if they exist) at the end to visually compare them
% delete(findall(0,'Type','figure'))
% plotClusterResults;

fprintf('RNG seed used: %f\n\n\n', rseed);