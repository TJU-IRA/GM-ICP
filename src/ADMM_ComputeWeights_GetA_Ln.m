function [A_hat, b_hat, xx, SortWeights] = ADMM_ComputeWeights_GetA_Ln(p1,p2,n1,n2, alpha1, u1, u2,u3, PP_norm)
n = n2; 
d = p1 - p2;
dist = dot(d,n)';
b = abs(dist);

w = sqrt(sum((p1-p2).^2))';
SortWeights = sort(w); 
kernelWidth = alpha1* (u1^PP_norm);
xx = length(find(SortWeights <(u3*3*sqrt(3))) )/ size(w,1);
weights1 = exp(-b.^PP_norm/kernelWidth);
weightss = weights1.^(1/PP_norm);
%% 
c = cross(p1,n);
A = [c' n'];
WW =repmat(weightss, 1, 6);
A_hat = A.*WW;
b_hat = dist.*weightss;
end

