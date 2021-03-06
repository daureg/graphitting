function [w, A, H, f] = compute_graph(X, kind, mew)
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
w = zeros(m, 1);
% Fortunately, she will be helped by some friends, although they shall be
% presented later, as they are, with all due respect, mainly calculations'
% artifice (and as such, they don't have a clue about how to start).
M = sparse(d*n, m);
MAX_ITER = 70;
nb_iter = 1;
bin_upper = n*(0:n-1) - cumsum(0:n-1);
considered_last_time = [];
[HK, UK] = get_complete_matrices(X);
AK = abs(UK);
if strcmpi(kind, 'soft')
	HK = [HK+mew*(AK'*AK) -mew*AK'; -mew*AK mew*eye(n)];
	fK = -mew*[ones(1, n)*AK -ones(1, n)]';
end
all_edges = 1:m;

% One day, an old man told $X$ that she could for instance start with this
% random small subset: a third of $(d+1)n$ edges, as this was indeed common
% knowledge, provided by a book called Theorem 3.1, that $U$ cannot contains
% more. At this point, the astute reader may wonder why we do not use a more
% sensible initial choice like bind each node with its closest neighbor. Well,
% I am telling the story so we do like this. But feel free to contribute!
% load('edges.mat');
% sol = ed;
MAX_EDGES = (d+1)*n;
edges = randi(m, 1, floor(0.07*MAX_EDGES));
% maintain the list of edges currently in the graph
current_edges = unique(edges);
% what proportion of edges are allowed above the upper bound
overflow_factor = 1.1;
% what proportion of <0 derivates edges we add at most
adding_factor = 1;
% edges=sol;
to_remove = [];
% shamelessly cheating!
% load('edges.mat');
% % edges = sort([18 edges]);
% tmp = randperm(numel(edges));
% cedges = sort(edges(tmp(1:end-3)));
% save cedges cedges
% edges = cedges;

increment = logspace(-1, 0, 10);%.1:.05:1;
% And thus begin the quest of $X$, until she can not add more edges to $U$ or
% until she get fed up and realize that organizing illegal fights of turtles is
% much more exciting than finding love.
while (numel(edges) > 0 && nb_iter <= MAX_ITER)
	% $U$ was also quite stubborn and felt that linear indexing of edges
	% was not doing justice to his amazing 2D abilities. Therefore $X$ has
	% to resort to her cunning to convert them. The edges $1$ through $n-1$
	% were from $i=1$, those from $n$ to $n+(n-2)-1=2n-3$ started at $i=2$
	% and so on.  It turns out that bin\_upper has memorized all these
	% upper bounds so finding $i$s was simply a matter of finding the
	% maximum possible bounds.
	[positive, negative] = from_edges_to_index(edges, bin_upper, size(U));
	U(positive) = 1;
	U(negative) = -1;
	fprintf('%d: remove %d edges because of small weigth\n', nb_iter, length(to_remove))
	[is, js] = from_edges_to_index(to_remove, bin_upper, size(U));
	U([is js]) = 0;
	A = abs(U);
	i = 1;
	wnz = sum(w>1e-6)/m;
	real = find(sum(A, 1) > 1);
	% sprintf('%d, %d', nnz(A)/2, length(current_edges))
	disp(setdiff(real, current_edges))
	disp(setdiff(current_edges, real))
	disp(find(real ~= current_edges))
	assert(nnz(A)/2 == length(current_edges), 'c_e does not agree with A')
	assert(nnz(A)/2 <= ceil(overflow_factor*MAX_EDGES), 'check your budget');
	% while (nnz(A)/2 > MAX_EDGES*overflow_factor && i<10)
	% 	% I should have commentated while writing because even only one half
	% 	% hour after, I'm already not sure what I meant. But the main purpose
	% 	% is too keep removing small weights until we are below the maximum
	% 	% number of edges. One clever way of avoiding that would probably be
	% 	% not add that many edges in the first place, but only the more
	% 	% negative one.
	% 	% warning('There were too many edges at iteration %d so I removed the %.2f%% smallest of them (I''m not heartless, but science must prevail!)', nb_iter,100*wnz+increment(i))
	% 	to_remove = all_edges(w<quantile(w(w>1e-6), wnz+increment(i)/100));
	% 	w(to_remove) = 0;
	% 	[is, js] = from_edges_to_index(to_remove, bin_upper, size(U));
	% 	U([is js]) = 0;
	% 	A = abs(U);
	% 	i = i+1;
	% end
	% if i > 1
	% 	fprintf('removal\n')
	% end
	% assert(i<10, 'failed to remove enough edge')


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
		Yk = spdiags(T(:,k), 0, m, m);
		M(first_row:last_row, :) = U*Yk;
	end
	H = 2.*(M'*M);
	if strcmpi(kind, 'soft')
		H = [H+mew*(A'*A) -mew*A'; -mew*A mew*eye(n)];
		% remove the factor 2 as MATLAB optimize 1/2 x'*H*x + f'*x
		f = -mew*[ones(1, n)*A -ones(1, n)]';
	end

	% Having done all this preparatory work, $X$ could finally go see an
	% oracle living in the mountain, the so called quadprog, and ask him to
	% set $w$ optimally according to $M$. (Actually, she had also heard of
	% another one, SDPT3, potentially faster and able to deal with sparse
	% matrix instead of converting M'*M to a full one and taking $2n^4$
	% bytes of memory. But she had to ask her question in a slightly
	% different language:
	% http://www.math.nus.edu.sg/~mattohkc/sdpt3/guide4-0-draft.pdf
	if (strmatch('2013', version('-release')))
		o = optimoptions(@quadprog, 'Algorithm', 'interior-point-convex', 'MaxIter', 300, 'Display', 'off', 'TolFun', 1e-11);
	else
		o = optimset('Algorithm', 'interior-point-convex', 'MaxIter', 500, 'Display', 'off');
	end
	if strcmpi(kind, 'hard')
		% In one method, she had to ensured that the weighted sum of
		% degree was at least one for each node, or in the language of
		% the oracle: $-Aw \leq -\bm{1}$. But this was only possible
		% for the nodes that were part of at least one edge, that is
		% for the non zero rows of $A$).
		% y=sdpvar(m,1);
		% C=[y>=0, A*y-ones(n,1)>=0];
		% O=y'*(H)*y;
		% os=sdpsettings('solver', 'sdpt3', 'savesolveroutput', 1, 'usex0', 1, 'verbose', 0);
		% tic;
		% sol=solvesdp(C,O,os);
		% w = double(y);
		% toc
		% f=-mean(sol.solveroutput.obj);
		% do only one iteration since we need a way to compute the derivative
		% break;
		% nb_iter
		% nnz(A)/2
		[w, f, flag, output, lambda] = quadprog(H, sparse(m, 1), -A, -(sum(A, 2)>0), [], [], zeros(m,1), [], w, o);
		% [w, f, flag, output, lambda] = quadprog(H, sparse(m, 1), -A, -ones(n, 1), [], [], zeros(m,1), [], w, o);
		z = lambda.ineqlin;
		derivative = 2*HK*w - AK'*z;
		save('out.mat', 'lambda', 'derivative');
		% break;
	else
		% There was another method were a portion $\alpha$ of the nodes
		% were allowed to have degree less than one. But she still
		% has to think about to formulate
		% $ \min_{w,s} ||Mw||^2 + \mu||\bm{1} - Aw + s||^2 $
		% for quadprog or lsqnonneg (http://math.stackexchange.com/q/545280)
		% After carefull thinking, she gets at least a formulation for quadprog:
		[w, f, flag, output, lambda] = quadprog(H, f, [], [], [], [], zeros(m+n,1), [], w, o);
		derivative = HK*w + fK + lambda.lower;
		s = w(m+1:end);
		w = w(1:m);
		alt = HK(1:m,1:m)*w - mew*AK'*(s + ones(n, 1)) + lambda.lower(1:m);
		derivative = derivative(1:m);
		% norm(alt - derivative)
		% sum(derivative(sol)<-1e-5)
		save('out.mat', 'lambda', 'derivative');
		% if nb_iter==15
		% 	break;
		% end
		% error(strcat(kind, ' is not yet implemented'));
	end
	% Because the new $(w, z)$ were supposed to be feasible solution,
	% $\frac{d\Lambda}{d w}$ has to be positive. Therefore, she finds the
	% edges where it was not the case to add them in the next step.
	[~, idx] = sort(derivative);
	nb_neg = sum(derivative<0);
	if (nb_neg > 0)
		may_be_added = setdiff(idx(1:nb_neg), current_edges)';
	else
		%TODO: actually, maybe we just need another round with some 0
		%edges removed?
		may_be_added = [];
	end
	% This was badly erroneous
	% [val, may\_be\_added]=sort(derivative(find(derivative<0)));

	may_be_added = find(derivative<0)';
	% perform cheap regularization, namely clamp oddly large values
	% (find another way to do it in general)
	w(w>20) = 1;
	% Of course maybe there was nothing to do. Or more concerning, the
	% oracle was rambling and returned a solution that yields the same set
	% of edges to add as previously, in which case there was no point in
	% continuing any further.
	% She first remove edges that were probably numerical zero
	real = find(sum(A, 1) > 1);
	assert(all(real == current_edges), 'c_e diverge from real')
	to_remove = find(w<1e-7)';
	w(to_remove) = 0;
	to_remove = intersect(real, to_remove);
	fprintf('%d: I will remove %d of the current edges because their weigths are negative or too small\n', nb_iter, length(to_remove))
	current_edges = setdiff(current_edges, to_remove);
	if ((isempty(may_be_added)) || (length(may_be_added) == length(considered_last_time) && all(may_be_added == considered_last_time)))
		break;
	end
	% Then she decides to add only half of them but probably there were other
	% ways of doing it (like adding the "smallest one" ?)
	budget = ceil(MAX_EDGES*overflow_factor) - length(current_edges);
	fprintf('%d: graph has %d edges (max = %d), so budget = %d\n', nb_iter, length(current_edges), ceil(MAX_EDGES*overflow_factor), budget)
	allowed = floor(min(budget, length(may_be_added))*adding_factor);
	fprintf('%d: %d candidates (and a_f=%.1f) so allowed = %d\n', nb_iter, length(may_be_added), adding_factor, allowed)
	% here we add the more negative ones, farthest from 0
	edges = may_be_added(1:allowed);
	% here we add the negative ones closest to 0
	% edges = may_be_added(end:-1:end-allowed+1);
	current_edges = union(current_edges, edges);
	% edges = may\_be\_added(1:max(1, floor(end/2)))';
	% edges = may_be_added;
	% edges = sol;
	% sprintf('%.3f\t%.3f', mean(derivative(setdiff(may_be_added, sol))), mean(derivative(sol)))
	% this is quite arguable
	w(edges) = mean(w);
	w(edges) = w(edges) + normrnd(0, mean(w), length(edges), 1);
	considered_last_time = may_be_added;
	nb_iter = nb_iter + 1;
	fprintf('iter: %d\n', nb_iter)
end
% When $X$ finds the perfect weights for her graph (and hopefully not because
% she just give up), she have to fill some paperwork like computing weighted
% degree and Laplacian to make their union official.
% TODO: When it will work, use this as output argument
if strcmpi(kind, 'soft')
	sprintf('%f\t%f', norm(HK(1:m,1:m)*w)^2, mew*norm(ones(n,1) - A*w +s)^2)
	Aw = A;
else
	Aw = A*w;
end
W = spdiags (w, [0], m, m);
L = U*W*U';
sprintf('use %d out %d possible iterations', nb_iter-1, MAX_ITER)
end
