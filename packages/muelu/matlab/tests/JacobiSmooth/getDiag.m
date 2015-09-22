%Use getDiag() as "setup" function for MatlabSmoother
function D = getDiag(A)
	D = sparse(diag(diag(A))); %first call generates vector of diagonal, second call puts that in a matrix
end
