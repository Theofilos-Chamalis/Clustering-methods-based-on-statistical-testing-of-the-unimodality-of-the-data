function [X, df_set] = student_t_cluster (n, d, df_set, mindf, maxdf)

if (~isempty(df_set)) 
    if (numel(df_set) > d)
        disp('error');
    else
        X = trnd(repmat(df_set, n, 1), n, d);
    end
else
    df_set = mindf + rand(1,d) * (maxdf - mindf);
    X = trnd(repmat(df_set, n, 1), n, d);
end
plot(X(:,1),X(:,2), '.')