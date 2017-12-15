%%
%*************************************************************************************
%----rpca_lsd main function------
%--framework£º
% 1£ºUsing DR method to solve L£¨A£©, low rank part (background)
% 2£ºUsing ProxFlow method to solve S£¨E£©, sturcture sparsity part (foreground)
%*************************************************************************************
% The lansvd function can be downloaded here
% http://soi.stanford.edu/~rmunk/PROPACK/

function [A,E,tol_iter] = rpca_lsd(D,graph)
addpath PROPACK;
[n,p] = size(D); d=min(n,p);
% parameter setting
%kappa = 1.1/100;%can be tuned
kappa = 1.1; 
%tau = 0.61;%can be tuned
tau = 1;
lambda = tau*kappa;
eta = (1-tau)*kappa;
mu = 30/norm(sign(D)); 
rho = 1.1; 
alpha = 1;
beta = .2;
 
Energy = 0;

Y = zeros(size(D));% init  Lagrange multipliers 
E = zeros(size(D));% init   
A = D;

 %init threshold of convergence 
tol_inner1 = 1e-4; 
tol_inner2 = 1e-6;  
tol_out = 1e-3; 
err_out = 1;

%init maximun iterations 
MaxIter_out = 500;
MaxIter_inner1 = 1;
MaxIter_inner2 = 20;

iter_out = 0;
sv = 10;
tol_iter = 0;

while iter_out < MaxIter_out && err_out > tol_out
    iter_out = iter_out + 1;

    Ak = A; Ek = E;
    iter_inner1 = 0;
    err_inner1 = 1;

    while iter_inner1 < MaxIter_inner1 && err_inner1 > tol_inner1
        iter_inner1 = iter_inner1 + 1;

        G = D - Ek + Y/mu;
        Akk = G;
        Ahk = zeros(size(Akk));
        err_inner2 = 1;
        iter_inner2 = 0;

        while iter_inner2 < MaxIter_inner2 && err_inner2 > tol_inner2
            iter_inner2 = iter_inner2 + 1;
            if choosvd(d, sv) == 1
                [U Si V] = lansvd(Akk, sv, 'L');
            else
                [U Si V] = svd(Akk, 'econ');
            end
            diagS = diag(Si);
            diagS = diagS(1:sv);
            svn = length(find(diagS > beta));
            svp = svn;
            ratio = diagS(1:end-1)./diagS(2:end);
            [max_ratio, max_idx] = max(ratio);
            if max_ratio > 2
                svp = min(svn, max_idx);
            end  
            if svp < sv
                sv = min(svp + 1, d);
            else
                sv = min(svp + 10, d);
            end
     
            Ahk = U(:, 1:svp) * diag(diagS(1:svp) - beta) * V(:, 1:svp)';

            B = 2*Ahk - Akk + mu*beta*G;
            ns = norm(B);
            %B = solve_l1l2(B/(1 + mu*beta),beta*eta/(1+beta*mu));
            B = bsxfun(@times,B/(1 + mu*beta),max(0,1 - beta*eta./ns));
            Akk = Akk + alpha*(B - Ahk);
            err_inner2 = alpha*norm(B-Ahk,'fro');

            tol_iter = tol_iter + 1;

        end
        G = D - Ahk + Y/mu;
       % ns = norm(G);
       % Ep = bsxfun(@times, G, max(0,1-lambda/mu./ns));
       % Ep = solve_lstruct(G,0.01*lambda/mu,graph);
        Ep = solve_lstruct(G,0.01*lambda/mu,graph);
        err_inner1 = max(norm(Ek-Ep,'fro'),norm(Ak-Ahk,'fro'));
        Ek = Ep;
        Ak = Ahk;

    end

    A = Ak; E = Ek;
    err_out = norm(D-A-E,'fro')/norm(D,'fro');


    Y = Y + mu*(D - A - E);
    mu = rho*mu;
   % disp(['Iteration' num2str(iter_out) ' #svd ' num2str(tol_iter) ' r(A) ' num2str(svn)...
   %     ' |E|_0 ' num2str(length(find(abs(E)>0)))...
   %     ' stopCriterion ' num2str(err_out)]);
    disp(['Iteration' num2str(iter_out) ' #svd ' num2str(tol_iter) ' r(A) ' num2str(rank(A))...
        ' |E|_0 ' num2str(length(find(abs(E)>0)))...
        ' stopCriterion ' num2str(err_out)]);
     %revise
    Eold = Energy;
    Energy = norm(A,2)+norm(E,2);
    disp(['Iter:    ' num2str(iter_out) '   E:  ' num2str(Energy)]);
    if (abs(Eold-Energy)/(Eold+eps)<1e-8)
        break
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
% graph parameter
graph_param.regul='graph';
graph_param.lambda= lambda; % regularization parameter
graph_param.num_threads=-1; % all cores (-1 by default)
graph_param.verbose=true;   % verbosity, false by default
graph_param.pos=false;       % can be used with all the other regularizations
graph_param.intercept=false; % can be used with all the other regularizations
x = mexProximalGraph(w,graph,graph_param);

