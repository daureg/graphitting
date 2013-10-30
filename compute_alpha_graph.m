function [w, Aw, L] = compute_alpha_graph(X, alpha, tol)
[n, d] = size(X);
m = nchoosek(n, 2);
mu = 5*rand();
tau0 = 1.5;
MAX_ITER = 50;
% Set $\lambda$ so that $\tau\geq1$ with equality at MAX\_ITER.
lambda = (tau0 - 1)/MAX_ITER;
can_improve = true;
while (can_improve)
	w, Aw, L = compute_graph(X, 'soft', mu);
	tmp = max(zeros(m, 1), ones(m, 1) - A*w);
	alpha_bar = tmp'*tmp/n;
	% Maybe we don't need this complication and keep a fixed $\tau=tau_0$.
	tau = tau0/(1 + nb_iter*lambda);
	if (alpha_bar < alpha)
		% We want to increase $\bar{\alpha}$ so we need decrease $\mu$, which in
		% turn require $\tau\leq 1$
		tau = 1/tau;
	end
	% It is supposed to correspond to: "we then adjust $\mu$ up or down
	% proportionally to how far $\frac{\eta(w)}{n}=\bar{\alpha}$ is from the
	% desired value of $\alpha$."
	mu = tau*abs(alpha_bar - alpha)/alpha
	nb_iter = nb_iter + 1;
	can_improve = abs(alpha_bar - alpha) < tol && nb_iter < MAX_ITER;
end
end
