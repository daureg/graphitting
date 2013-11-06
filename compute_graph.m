function [w, A, H, f] = compute_graph(X, kind, mu)
% All stories have a hero and our is not different. So meet $X$, a handsome and
% brave set of $n$ $d$-dimensional vectors.
[n, d] = size(X);
m = nchoosek(n, 2);
% She is in love with the equally beautiful $U$, a matrix containing all the
% $m$ possible edges between $X$'s nodes.
U = sparse(n, m);
% Yet for now, $U$, like many boys of his age, is quite empty. Therefore $X$,
% must fill him with the weights in $w$. But to be fair, she has no feasible
% solution to propose so far.
w = [];
% Fortunately, she will be helped by some friends, although they shall be
% presented later, as they are, with all due respect, mainly calculations'
% artifice (and as such, they don't have a clue about how to start).
M = sparse(d*n, m);
MAX_ITER = 50;
nb_iter = 0;
bin_upper = n*(0:n-1) - cumsum(0:n-1);
considered_last_time = [];
[HK, UK] = get_complete_matrices(X);
AK = abs(UK);

% One day, an old man told $X$ that she could for instance start with this
% random small subset: a third of $(d+1)n$ edges, as this was indeed common
% knowledge, provided by a book called Theorem 3.1, that $U$ cannot contains
% more. At this point, the astute reader may wonder why we do not use a more
% sensible initial choice like bind each node with its closest neighbor. Well,
% I am telling the story so we do like this. But feel free to contribute!
edges = randi(m, 1, floor(0.33*(d+1)*n));
% shamelessly cheating!
load('edges.mat');

% And thus begin the quest of $X$, until she can not add more edges to $U$ or
% until she get fed up and realize that organizing illegal fights of turtles is
% much more exciting than finding love.
while (numel(edges) > 0 && nb_iter < MAX_ITER)
	% $U$ was also quite stubborn and felt that linear indexing of edges
	% was not doing justice to his amazing 2D abilities. Therefore $X$ has
	% to resort to her cunning to convert them. The edges $1$ through $n-1$
	% were from $i=1$, those from $n$ to $n+(n-2)-1=2n-3$ started at $i=2$
	% and so on.  It turns out that bin\_upper has memorized all these
	% upper bounds so finding $i$s was simply a matter of finding the
	% maximum possible bounds.
	vertex_i = arrayfun(@(x) find(x' <= bin_upper, 1, 'first'), edges) - 1;
	% Then $j$ follows
	vertex_j = vertex_i + edges - bin_upper(vertex_i);
	% and the pairs $(i,j)$ could be converted into $U$ indexes
	positive = bsxfun (@(x,y) sub2ind(size(U), x, y), vertex_i, edges);
	negative = bsxfun (@(x,y) sub2ind(size(U), x, y), vertex_j, edges);
	% in order to represent the newly selected edges.
	U(positive) = 1;
	U(negative) = -1;
	A = abs(U);
	assert(sum(A(:))/2 <= (d+1)*n, 'there are too many edges');

	% To assess their compatibility, the tradition was to compute the
	% Frobenius norm of $X$ times the graph Laplacian. Both find this
	% method awkward and they choose to reformulate it as an Euclidean
	% distance. But doing so require the help of a friend: $M$.
	T = U'*X; % $y^{(k)} = U^Tx^{k}$ is thus the $k$th column of $T$
	% First $X$ mixes her columns with $U$, producing $d$ new vectors of
	% length $n$: $y^{(k)}$.
	for k=1:d
		first_row = 1 + (k-1)*n;
		last_row = n + (k-1)*n;
		% These new vectors were soon promoted as diagonal matrices and
		% filled $M$ from top to bottom (although it would have been
		% faster to do it in parallel).
		Yk = spdiags(T(:,k), [0], m, m);
		M(first_row:last_row, :) = U*Yk;
	end
	H = M'*M;

	% Having done all this preparatory work, $X$ could finally go see an
	% oracle living in the mountain, the so called quadprog, and ask him to
	% set $w$ optimally according to $M$. (Actually, she had also heard of
	% another one, SDPT3, potentially faster and able to deal with sparse
	% matrix instead of converting M'*M to a full one and taking $2n^4$
	% bytes of memory. But she had to ask her question in a slightly
	% different language:
	% http://www.math.nus.edu.sg/~mattohkc/sdpt3/guide4-0-draft.pdf
	if strcmpi(kind, 'hard')
		% In one method, she had to ensured that the weighted sum of
		% degree was at least one for each node, or in the language of
		% the oracle: $-Aw \leq -\bm{1}$. But this was only possible
		% for the nodes that were part of at least one edge, that is
		% for the non zero rows of $A$).
		tic;
		y=sdpvar(m,1);
		C=[y>=0, A*y-ones(n,1)>=0];
		O=y'*H*y;
		os=sdpsettings('solver', 'sdpt3', 'savesolveroutput', 1);
		sol=solvesdp(C,O,os);
		w = double(y);
		toc
		f=-mean(sol.solveroutput.obj);
		% do only one iteration since we need a way to compute the derivative
		break;
		o = optimoptions(@quadprog, 'Algorithm', 'active-set', 'MaxIter', 500);
		[w, f, flag, output, lambda] = quadprog(H, sparse(m, 1), -A, -(sum(A, 2)>0), [], [], zeros(m,1), [], w, o);
		toc
		z = lambda.ineqlin;
		derivative = 2*HK*w + AK'*z;
		f
		save('out', 'H', 'w', 'A', 'lambda', 'derivative');
	else
		% There was another method were a portion $\alpha$ of the nodes
		% were allowed to have degree less than one. But she still
		% has to think about to formulate
		% $ \min_{w,s} ||Mw||^2 + \mu||\bm{1} - Aw - s||^2 $
	 	% for quadprog or lsqnonneg (http://math.stackexchange.com/q/545280)
		error(strcat(kind, ' is not yet implemented'));
	end
	% Because the new $(w, z)$ were supposed to be feasible solution,
	% $\frac{d\Lambda}{d w}$ has to be positive. Therefore, she finds the
	% edges where it was not the case to add them in the next step.
	% [val, may_be_added]=sort(derivative(find(derivative<0)));
	% This was badly erroneous
	may_be_added = find(derivative<0)'
	% perform cheap regularization, namely remove weirdly large value
	% (find another way to do it in general)
	w(w>1e6) = 0;
	% Of course maybe there was nothing to do. Or more concerning, the
	% oracle was rambling and returned a solution that yields the same set
	% of edges to add as previously, in which case there was no point in
	% continuing any further.
	if ((isempty(may_be_added)) || (length(may_be_added) == length(considered_last_time) && all(may_be_added == considered_last_time)))
		break;
	end
	% She decide to add only half of them but probably there were other
	% ways of doing it (like adding the "smallest one" ?)
	% edges = may_be_added(1:max(1, floor(end/2)))';
	edges = may_be_added;
	% this is quite arguable
	w(may_be_added) = mean(w);
	considered_last_time = may_be_added;
	nb_iter = nb_iter + 1;
end
% When $X$ finds the perfect weights for her graph (and hopefully not because
% she just give up), she have to fill some paperwork like computing weighted
% degree and Laplacian to make their union official.
Aw = A*w;
W = spdiags (w, [0], m, m);
L = U*W*U';
end
