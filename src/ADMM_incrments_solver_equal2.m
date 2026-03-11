function T = ADMM_incrments_solver_equal2(A_hat, b_hat)
% 计算运行矩阵 p = 2.0的情况
A_pinv = -pinv(A_hat);
x_t= A_pinv*b_hat;
rot = x_t(1:3); trans = x_t(4:6);
rotangle = norm(rot);
TR = rotation_matrix(rotangle, rot);
T = [TR(1:3,1:3) trans; 0 0 0 1];
end
