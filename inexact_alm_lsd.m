function [A,E,tol_iter] = inexact_alm_lsd(D,graph)
addpath PROPACK;
[n,p] = size(D); d=min(n,p);

lambda = 1 / sqrt(n);

% initialize
Y = D;
norm_two = lansvd(Y, 1, 'L');
norm_inf = norm( Y(:), inf) / lambda;
dual_norm = max(norm_two, norm_inf);
Y = Y / dual_norm;

mu = 12.5/norm_two; % this one can be tuned

rho = 1.5; 

% Y = zeros(size(D));
E = zeros(size(D));
A = zeros(size(D));

tol_out = 1e-7; %1e-7

converged = false;
iter_out = 0;
sv = 10;
tol_iter = 0;

while ~converged
    iter_out = iter_out + 1;

    G = D - E + Y/mu;
    if choosvd(d, sv) == 1
        [U Si V] = lansvd(G, sv, 'L');
    else
        [U Si V] = svd(G, 'econ');
    end
    diagS = diag(Si);
%    svp = length(find(diagS > 1/mu));
    diagS = diagS(1:sv);
    svn = length(find(diagS > 1/mu));
    svp = svn;
    ratio = diagS(1:end-1)./diagS(2:end);
    [max_ratio, max_idx] = max(ratio);
    if max_ratio > 2
        svp = min(svn, max_idx);
    end
    if svp < sv
       sv = min(svp + 1, d);
    else
       sv = min(svp + round(0.05*d), d);
    end  
%      
    A = U(:, 1:svp) * diag(diagS(1:svp) - 1/mu) * V(:, 1:svp)';
    
    G2 = D - A + Y/mu;
    E = solve_lstruct(G2,lambda/mu,graph);
    Z = D - A - E;
    Y = Y + mu*Z;
    mu = min(mu*rho, mu * 1e7);
    err_out = norm(Z,'fro')/norm(D,'fro');


   % disp(['Iteration' num2str(iter_out) ' #svd ' num2str(tol_iter) ' r(A) ' num2str(svn)...
   %     ' |E|_0 ' num2str(length(find(abs(E)>0)))...
   %     ' stopCriterion ' num2str(err_out)]);
    disp(['Iteration' num2str(iter_out) ' r(A) ' num2str(svp)...
        ' |E|_0 ' num2str(length(find(abs(E)>0)))...
        ' stopCriterion ' num2str(err_out)]);
    if err_out < tol_out
        converged = true;
    end    
end

function [E] = solve_lstruct(W,lambda,graph)
n = size(W,2);
E = W;
%for i=1:n
    %E(:,i) = solve_l2(W(:,i),lambda,graph);
E = solve_ls(W,lambda,graph);
%end

function [x] = solve_ls(w,lambda,graph)
% min lambda |x|_2 + |x-w|_2^2
% graph paramters
graph_param.regul='graph';
graph_param.lambda= lambda; % regularization parameter
graph_param.num_threads=-1; % all cores (-1 by default)
graph_param.verbose=true;   % verbosity, false by default
graph_param.pos=false;       % can be used with all the other regularizations
graph_param.intercept=false; % can be used with all the other regularizations
x = mexProximalGraph(w,graph,graph_param);