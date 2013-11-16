function [w, Aw, L, report] = compute_alpha_graph(X, alpha, tol, cheat)
[n, ~] = size(X);
m = nchoosek(n, 2);
mew = normrnd(1, .2);
% tau0 = 1.5;
nb_iter = 1;
MAX_ITER = 13;
% Set $\lambda$ so that $\tau\geq1$ with equality at MAX\_ITER.
% lambda = (tau0 - 1)/MAX\_ITER;
can_improve = true;
report = zeros(MAX_ITER, 4);
while (can_improve)
	if cheat
		[w, A, H, f, L] = fully_solve(X, 'soft', mew);
	else
		[w, A, L] = compute_graph(X, 'soft', mew);
	end
	tmp = max(zeros(n, 1), ones(n, 1) - A(:,1:m)*w);
	alpha_bar = tmp'*tmp/n;
	% Maybe we don't need this complication and keep a fixed $\tau=\tau_0$.
	% tau = tau0/(1 + nb\_iter*lambda);
	delta = (alpha_bar - alpha)/alpha;
	% if (alpha\_bar < alpha)
		% We want to increase $\bar{\alpha}$ so we need decrease $\mu$,
		% which in turn require $\tau\leq 1$.
	% 	tau = 1/tau;
	% end
	% It is supposed to correspond to: "we then adjust $\mu$ up or down
	% proportionally to how far $\frac{\eta(w)}{n}=\bar{\alpha}$ is from the
	% desired value of $\alpha$."
	rep = [alpha_bar abs(alpha_bar - alpha)/alpha mew];
	% actually, it's a bit weird to reduce $\mu$ when $\delta$ is
	% negative, because it happens when $\bar{\alpha} < \alpha$, meaning
	% that the weights satisfy the constraints even "more" that what we
	% are looking for.
	mew = ((abs(delta)+1)^sign(delta))*mew;
	rep = [rep mew];
	report(nb_iter, :) = rep;
	nb_iter = nb_iter + 1;
	can_improve = (abs(alpha_bar - alpha) > tol) && (nb_iter <= MAX_ITER);
end
Aw = A(:,1:m)*w;
% sprintf('\ta_bar\trel\ttau\tmu_n\tmu_n+1')
% report
end
