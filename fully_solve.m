function [w, A, H, f, L] = fully_solve(X)
n = size(X, 1);
m = nchoosek(n, 2);
w = [];
[H, U] = get_complete_matrices(X);
A = abs(U);
o = optimoptions(@quadprog, 'Algorithm', 'active-set', 'MaxIter', 400);
[w, f, flag, output, lambda] = quadprog(H, zeros(m, 1), -A, -ones(n, 1), [], [] , zeros(m, 1), [] , w, o);
L = U*diag(w)*U';
end
