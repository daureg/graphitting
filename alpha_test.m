alpha = .05;
load small.mat
tic;
[wr, Ar, Lr, rep] = compute_alpha_graph(X, alpha, .0005);
toc
