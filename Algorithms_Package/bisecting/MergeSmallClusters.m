function [gIdx_ref_Init,c_all,initclNo] = MergeSmallClusters(gIdx_ref_Init,c_all,smallest_cluster,X)

    repz = 1;
    az = unique(gIdx_ref_Init);
    out = [az,histc(gIdx_ref_Init(:),az)];
    smallcl=[0];
    
    for i1=1:length(out)
        if(out(i1,2)<smallest_cluster)
            if(smallcl(1)==0)
                smallcl(1)=i1;
            else
                smallcl=[smallcl;i1];
            end        
        end
    end
    
    while(smallcl(1,1)~=0)
        % Check if there are any small clusters 
        if(smallcl(1,1)~=0)
            smallcl(:,2)=0;
            smallclsize = size(smallcl);

            for i1=1:smallclsize(1,1)
                    q = pdist2(c_all(smallcl(i1,1),:),c_all);
                    q(smallcl(i1)) = 999;
                    [~,smalldist] = min(q);
                    smallcl(i1,2) = smalldist;
                    smallcl2(i1,1) = smallcl(i1,2);
                    smallcl2(i1,2) = smallcl(i1,1);
            end

            % Remove the duplicates from smallcl
             [La,Lb]=ismember(smallcl,smallcl2,'rows');
             samepairs = [];

             for i1=1:length(Lb)
                 if(Lb(i1)~=0)
                     if(Lb(i1)<=Lb(Lb(i1)))
                         samepairs = [samepairs; Lb(i1),Lb(Lb(i1))];
                     end
                 end
             end

            if(~isempty(samepairs))
               samepairs=sort(samepairs,1,'descend');
               samepairs(:,1) = samepairs(:,2);
               samepairs(:,2) = [];
               for i1=1:length(samepairs)
                smallcl(samepairs(i1),:)=[];
               end
               smallcl = flipud(smallcl);
            else
                smallcl = flipud(smallcl);
            end

            smallclsize = [];
            smallclsize = size(smallcl);

            for i1=1:smallclsize(1,1)
                for i2=1:length(gIdx_ref_Init)
                    if(gIdx_ref_Init(i2) == smallcl(i1,1))
                        gIdx_ref_Init(i2) = smallcl(i1,2);
                        c_all(smallcl(i1,1)) = c_all(smallcl(i1,2));
                    end
                end
            end

            for i1=1:smallclsize(1,1)
                for i2=1:length(gIdx_ref_Init)
                    if(gIdx_ref_Init(i2)>smallcl(i1,1))
                        gIdx_ref_Init(i2) = gIdx_ref_Init(i2) - 1;
                    end
                end
                c_all(smallcl(i1),:)=[];
            end
        end
        
        [gIdx_ref_Init, c_all, ~, ~, ~] = ark_kmeans(X', [], c_all', 20, -1, 1, 0, 0, 0, 0);
        c_all = c_all';
        
        %Fill again smallcl variable to check for small clusters
        clear az out smallcl;

        az = unique(gIdx_ref_Init);
        out = [az,histc(gIdx_ref_Init(:),az)];
        smallcl=[0];

        for i1=1:length(out)
            if(out(i1,2)<smallest_cluster)
                if(smallcl(1)==0)
                    smallcl(1)=i1;
                else
                    smallcl=[smallcl;i1];
                end        
            end
        end
        
%         fprintf('\nThis is Repz number : %d \n',repz);
%         fprintf('\nThis what smallcl has : \n');
%         smallcl
        initclNo = length(unique(gIdx_ref_Init));
%         fprintf('\nThe number of clusters is : %d \n',initclNo);
        repz = repz+1;
    end
	
	initclNo = length(unique(gIdx_ref_Init));
	fprintf('\n\nAfter merging the small clusters, the number of clusters is now: %d\n\n',initclNo);
	
end
