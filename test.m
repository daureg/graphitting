load small.mat
tic;
w = fully_solve(X);
toc
tic;
[wi, ~, ~] = compute_graph(X, 'hard', 0);
toc
