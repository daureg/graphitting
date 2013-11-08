function [positive, negative] = from_edges_to_index(edges, bin_upper, sizeU)
	vertex_i = arrayfun(@(x) find(x' <= bin_upper, 1, 'first'), edges) - 1;
	% Then $j$ follows
	vertex_j = vertex_i + edges - bin_upper(vertex_i);
	% and the pairs $(i,j)$ could be converted into $U$ indexes
	positive = bsxfun (@(x,y) sub2ind(sizeU, x, y), vertex_i, edges);
	negative = bsxfun (@(x,y) sub2ind(sizeU, x, y), vertex_j, edges);
end
