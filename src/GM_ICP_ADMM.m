function [T,count] = GM_ICP_ADMM(SP,TP,SN,TN,alpha_, PP_norm)

T = eye(4);
Btree = KDTreeSearcher(TP');
[idx, dist] = knnsearch(Btree,TP','k',7);

dist = dist(:,2:7);
u3 = median(median(dist,2))/(3*sqrt(3)); 

x_nearest = zeros(size(idx, 1), 1);
for i = 1:size(idx, 1)
    k_dist = zeros(7,1);
    for j  = 1:7
        TP1 = TP(:, idx(i,1));
        TP2 = TP(:, idx(i, j));
        NP1 = TN(:, idx(i,1));
        k_dist( j) = abs(dot(TP1 - TP2, NP1));
    end
    x_nearest(i) = median(k_dist(2:7));
end
u2 = median(x_nearest)/6;
[p1,p2,n1,n2,r] = match_points(SP,TP,SN,TN,Btree);

u1 = 3*median(r);   
stop1 = 0; count = 0;

while(~stop1) 
    for i=1:100
        [A_hat, b_hat, xx, sortw]= ADMM_ComputeWeights_GetA_Ln(p1,p2,n1,n2, alpha_ ,u1, u2, u3, PP_norm);   
        if PP_norm <= 1.01
            T0 = ADMM_incrments_solver_small1(A_hat,b_hat, PP_norm);
        else
            if abs(PP_norm - 1.5) < 1e-2
                T0 = ADMM_incrments_solver_equal15(A_hat,b_hat, PP_norm);
            end
            if abs(PP_norm - 2.0) < 1e-2
                T0 = ADMM_incrments_solver_equal2(A_hat,b_hat);
            end
        end
   
        T = T0*T;
        p12 = T*[SP;ones(1,size(SP,2))]; p1 = p12(1:3,:);
        n1 = T(1:3,1:3)*SN;
        [p1,p2,n1,n2,~] = match_points(p1,TP,n1,TN,Btree);
        stop2 = norm(T0-eye(4));
        if stop2 < 1e-5
            break;
        end
    end
    if abs(u1-u2)<1e-6
        stop1 = 1;
    end
    if xx > 0.30
         u1_ = median(sortw(1:ceil(size(sortw,1)*xx)))/2;
         u1 = min([u1_, u1/4]);
         xx
    else
         u1 = u1/4;
    end
    count = count+i;
    if u1<u2
        u1 = u2;
    end
end

