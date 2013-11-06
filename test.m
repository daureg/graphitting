load small.mat
tic;
[wr, Ar, Hr, fr, Lr] = fully_solve(X);
toc
tic;
[wmw, Aw, Hw, fw] = compute_graph(X, 'hard', 0);
toc
