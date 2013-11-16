function [w, A, H, f, L, s] = fully_solve(X, kind, mew)
n = size(X, 1);
m = nchoosek(n, 2);
w = [];
[H, U] = get_complete_matrices(X);
f = sparse(m, 1);
A = abs(U);
A_constraint = -A;
b = -ones(n, 1);
lower_bound = zeros(m, 1);
if strcmpi(kind, 'soft')
	H = [H+mew*(A'*A) -mew*A'; -mew*A mew*eye(n)];
	% remove the factor 2 as MATLAB optimize 1/2 x'*H*x + f'*x
	f = -mew*[ones(1, n)*A -ones(1, n)]';
	A_constraint = [];
	b = [];
	lower_bound = zeros(m+n, 1);
end
if strmatch(version('-release'), '2013')
	o = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'MaxIter', 500, 'Display', 'off', 'TolFun', 1e-15);
else
	o = optimset('Algorithm', 'interior-point-convex', 'MaxIter', 500, 'Display', 'off', 'TolFun', 1e-15);
end
[w,~,~,~,lambda] = quadprog(H, f, A_constraint, b, [], [] , lower_bound, [] , w, o);
if strcmpi(kind, 'soft'); der =H*w + f ;s = w(m+1:end); w = w(1:m);
	save('out', 'der', 'lambda')
find(der(1:m)<-1e-6)
sprintf('%f\t%f', norm(H(1:m,1:m)*w)^2, mew*norm(ones(n,1) - A*w +s)^2)
end
L = U*diag(w)*U';
end
