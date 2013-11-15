% addpath(genpath('yalmip'))
% addpath('SDPT3-4.0')
% startup();
load small.mat
tic;
[wr, Ar, Hr, fr, Lr] = fully_solve(X, 'hard', 0);
toc
tic;
[ww, Aw, Hw, fw] = compute_graph(X, 'hard', 0);
toc
% derivative=(2*HK*ww)';
% load cedges.mat
% edges = find(wr>.0001);
% m=105;
% excluded=setdiff(1:m, edges);
% % d are the 3 edges that we have remove from the 32 correct ones.
% d=setdiff(edges,cedges)
% % t=ww; t(t>10)=0; der=(2*HK*t)';
% load out
% % so we'd like the derivative to be largely negative on those one, to indicate
% % us that we should add them next.
% [v,i]=sort(derivative);
% v(1:5)'
% setdiff(i(1:5)', cedges)
edges = find(wr>1e-4);
[wr(edges) ww(edges)]
norm(wr-ww)
