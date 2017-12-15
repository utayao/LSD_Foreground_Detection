
function Mask = ForegroundMask(S,M,L,idx)
    S = abs(S);
    S_back_temp = S<0.5*max(S(:)); 
    S_diff = abs(M-L).*S_back_temp;
    id = find(S_diff>0);
    mu_s = mean(S_diff(id));
    sigma_s = std(S_diff(id));
    th = mu_s + 2*sigma_s;
    Mask = S>th;
end

