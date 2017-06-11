fprintf('Optimal n is %d \n', sum(sum(x_opt)));
all_n = unique(history_all(1,:));
for n = all_n
    extract = history_all(:, history_all(1,:) == n);
    [betas, index] = sort(extract(7,:), 'ascend');
    SLs = extract(6,index);
    index_thres = sum(SLs<=0.8);
    f_beta = extract(3,index);
    diff_f_beta = diff(f_beta);
    if ~isempty(diff_f_beta)
        mono = all(diff_f_beta>0);
        fprintf('n = %d, f(beta)^n is monotonic?: %d \n', n, mono);
        if mono == 0
            betas
            f_beta
        end
        if index_thres > 1
            figure;
            hold on;
            plot(betas(1:end-2), f_beta(1:end-2));
            title(['Fix n = ', num2str(n), ', Bisection beta. When beta is between point ',  num2str(index_thres),  ' and ' , num2str(index_thres+1), ' (~= ', num2str(betas(index_thres)), ' to ',  num2str(betas(index_thres+1)), ' ), SL ~= 0.8']);
            xlabel('beta');
            ylabel('f^n(beta)');
            hold off;
        end
    end
end