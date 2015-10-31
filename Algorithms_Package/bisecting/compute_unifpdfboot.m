%------------
% This function generates boot_dips variable 
% and stores it in an external file named
% unifpdfbootext.mat for later use by 
% test_projected_unimodality.m to speed up
% the pdip and pdip reverse algorithms.
%------------
% Input
%	MaxN: The maximum number of points a cluster can have
%	nboot: Bootstrap sample size (default is 1000)   
% Output
%   	None.
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
%------------


function [] = compute_unifpdfboot(MaxN,nboot)
    % This function generates boot_dips variable
    % and stores it in an external file named
    % unifpdfbootext.mat. It is used offline.
    % Default MaxN value is 3000 and nboot 1000
    tic;
    boot_dip = zeros(MaxN,nboot);
    parfor j = 4:MaxN
        for i=1:nboot
                   unifpdfboot = sort(unifrnd(0,1,1,j));
                   unif_dip    = HartigansDipTest(unifpdfboot);
                   boot_dip_temp(i) = unif_dip;
        end;
        boot_dip_temp = sort(boot_dip_temp);
        boot_dip(j,:) = boot_dip_temp;
    end
    boot_dips = boot_dip;
    toc;
    save('unifpdfbootext.mat','boot_dips');
end
