function [H, U] = get_complete_matrices( X )
[n, d] = size(X);
m = nchoosek(n, 2);
M = sparse(d*n, m);
U = sparse(n, m);
% Despite MATLAB obnoxious warning, it's faster to index sparse matrix
% directly
% M = zeros(d*n, m);
% U = zeros(n, m);
komplete = [ones(1, n-1); -speye(n-1)];
from = 1;
to = n-1;
for v = 1:n
	U(v:n, from:to) = komplete(1:n-(v-1), 1:n-v);
	from = to + 1;
	to = from + n - (v+1) - 1;
end
% U = sparse(U);
T = U'*X;
for k=1:d
	first_row = 1 + (k-1)*n;
	last_row = n + (k-1)*n;
	Yk = spdiags(T(:,k), 0, m, m);
	M(first_row:last_row, :) = U*Yk;
end
% M = sparse(M);
H=M'*M;
end