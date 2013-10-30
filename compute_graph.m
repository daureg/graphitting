function [w, Aw, L] = compute_graph(X, kind, mu)
[n, d] = size(X);
m = nchoosek(n, 2);
M = sparse(d*n, m);
U = sparse(n, m);
w = [];
MAX_ITER = 50;
nb_iter = 0;

% At the end, we cannot have more than $\frac{d+1}{n}$ edges according to
% Theorem 3.1. But we must start with only a small subset of them. So we first
% select randomly 7\% of them (for no specific reason but cinematographic
% one). Actually, it may be more sensitive to use some kind of heuristic like
% nearest neighbors at this stage to be more efficient later.

edges = randi(m, 1, 0.07*(d+1)*n);

while (numel(edges) > 0 && nb_iter < MAX_ITER)
	% Remove the multiple of $n$ to avoid self loop (and the first one in case
	% it's 1).
	edges = edges(mod(edges, n) != 0)(2:end);
	% Then we update the corresponding element $(i,j)=e$ of $U$ with
	% respectively $1$ and $-1$.
	vertex_j = rem(edges, n);
	vertex_i = div(edges, n) + 1;
	positive = bsxfun (@(x,y) sub2ind(size(U), x, y), vertex_i, edges);
	negative = bsxfun (@(x,y) sub2ind(size(U), x, y), vertex_j, edges);
	U(positive) = 1;
	U(negative) = -1;
	A = abs(U);
	assert(sum(A)/2 <= (d+1)*n, 'there are too many edges');

	T = U'*X; % $y^{(k)} = U^Tx^{k}$ is thus the $k$th column of $T$
	% TODO: use parfor
	for k=1:d
		first_row = 1 + (k-1)*n;
		last_row = n + (k-1)*n;
		Yk = spdiags(T, [k], m, m);
		M(first_row:last_row, :) = U*Yk;
	end
	% Now that we have built our matrices, we can solve the minimization problem
	% TODO use SDPT3, although the documentation is quite intimidating:
	% http://www.math.nus.edu.sg/~mattohkc/sdpt3/guide4-0-draft.pdf
	if strcmpi(kind, 'hard')
		[w, f, flag, output, lambda] = quadprog(M'*M, sparse(m, 1), ...
																-A, -ones(m, 1), [], [], [], [], w);
		z = lambda.lower;
		derivative = 2*M'*M*w - A'*z;
	else
		% According to the paper, we want to solve
		% $ \min_{w,s} ||Mw||^2 + \mu||\mathbf{1} - Aw - s|| $
	 	% subject to $w,s\geq 0$, but I don't see how to formulate that for 
		% quadprog or lsqnonneg (see http://math.stackexchange.com/q/545280)
		error(strcat(kind, ' is not yet implemented'));
	end
	[val, may_be_added]=sort(derivative(find(derivative<0)));
	% The paper says: we add to our quadratic program the edges with the smallest
	% $\frac{d\Lambda}{d w_{i,j}}$ values, which I think mean not all. For now,
	% let's take half of them. TODO take only the one below average or look at
	% diff.
	edges = may_be_added(1:end/2);
	nb_iter = nb_iter + 1;
end
Aw = A*w;
L = U'*spdiags (w, [0], m, m)*U;
