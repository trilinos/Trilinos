%Prototype function to generate unsmoothed P from aggregates
function [P, Nullspace, KeepNodes] = createP(A, KeepNodesFine)
	N = size(A, 1);
	Nkeep = length(KeepNodesFine);
	OK_IDX = setdiff(1:N, KeepNodesFine);
	A_OK = A(OK_IDX,OK_IDX);
	h = muelu('setup', A_OK, 'max levels', 2, 'coarse: max size', int32(size(A, 1) / 3));
	P_OK = muelu('get', h, 1, 'P');
	muelu('cleanup', h);
	% NOTE To Ray: This probably won't work for multiple PDEs/NODE
	P = sparse(N, size(P_OK, 2) + Nkeep);
	P(OK_IDX, 1:size(P_OK, 2)) = P_OK;
	P(KeepNodesFine, size(P_OK, 2) + 1:end) = speye(Nkeep, Nkeep);
	% NOTE: Totally Faking nullspace
	Nullspace = ones(size(P_OK, 2), 1);
	KeepNodes=size(P_OK, 2) + [1:Nkeep];
end 
