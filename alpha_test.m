load small.mat
tic;
[wf, Af] = compute_graph(X, 'hard', 0);
toc
res=[]
tic;
for alpha=logspace(-1, -5, 10)
	alpha
	[wr, Aw, Lr, rep] = compute_alpha_graph(X, alpha, alpha/10, 1);
	res = [res norm(wf-wr)/norm(wf)];
end
toc
disp('as alpha goes to 0, so must the relative difference in norm between the hard solution and the alpha one')
res
alpha = 0.1;
% tic;
% [wr, Awr, Lr, repr] = compute_alpha_graph(X, alpha, alpha/10, 1);
% toc
% tic;
% [ww, Aww, Lw, repw] = compute_alpha_graph(X, alpha, alpha/10, 0);
% toc
% edges = find(wr>1e-5);
% % disp('the relative difference in norm between the full solution and the iterative one should be close to zero')
% norm(wr-ww)/norm(wr)

mew=1e4;
[wr, Ar] = fully_solve(X, 'soft', mew);
[ww, Aw] = compute_graph(X, 'soft', mew);
% norm(wr-ww)/norm(wr)
