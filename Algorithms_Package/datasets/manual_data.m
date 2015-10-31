function X = manual_data(X, n_iter, iter_var, make_fig) 

if isempty(make_fig)
    figure; plot(X(:,1), X(:,2), '.b'); hold on;
end

while (1) %for i=1:ceil(n_new/n_iter),
    [x,y] = ginput(1);
    new_points = [[x y]; repmat([x-iter_var/2 y-iter_var/2], n_iter-1, 1) + iter_var * rand(n_iter-1,2)];
    X = [X; new_points];
    plot(new_points(:,1), new_points(:,2), '.g');
    
    fprintf('Size of dataset: %g\n', size(X,1));
    reply = input('Do you want more? y/n [y]: ', 's');
    if (~isempty(reply) && reply == 'n')
        break
    end
end