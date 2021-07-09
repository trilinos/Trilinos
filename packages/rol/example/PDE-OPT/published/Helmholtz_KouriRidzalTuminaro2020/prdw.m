function [y] = prdw( x, usr_par)
%
%   compute and/or apply a perturbed Rees-Dollar-Wathen KKT preconditioner
%

  global GLB_PROB;

  nu = GLB_PROB.nu;
  nz = GLB_PROB.nz;

  beta = GLB_PROB.beta;
  K = GLB_PROB.K;
  M = GLB_PROB.M;
  C = GLB_PROB.C;
  R = GLB_PROB.R;
  J = GLB_PROB.J;
  A1 = K;

  A2 = ctranspose(A1);

  perturb = min(1e-4, beta);

  y = zeros(size(x));

  if (GLB_PROB.reuse_factors)
    if (~GLB_PROB.factors_computed)
      fprintf("  Computing factors (to be reused) ...\n  ");
      tic
      [GLB_PROB.LC, GLB_PROB.UC, GLB_PROB.PC, GLB_PROB.QC]     = lu(C+perturb*M);
      [GLB_PROB.LR, GLB_PROB.UR, GLB_PROB.PR, GLB_PROB.QR]     = lu(R);
      [GLB_PROB.LA1, GLB_PROB.UA1, GLB_PROB.PA1, GLB_PROB.QA1] = lu(A1);
      %[GLB_PROB.LA2, GLB_PROB.UA2, GLB_PROB.PA2, GLB_PROB.QA2] = lu(A2);
      toc
      GLB_PROB.factors_computed = true;
    end
    y(1:nu)         = GLB_PROB.QC * (GLB_PROB.UC \ (GLB_PROB.LC \ (GLB_PROB.PC * x(1:nu))));
    y(nu+1:nu+nz)   = (1/beta) * ( GLB_PROB.QR * (GLB_PROB.UR \ (GLB_PROB.LR \ (GLB_PROB.PR * x(nu+1:nu+nz)))) );
    tmp             = (C+perturb*M) * ( GLB_PROB.QA1 * (GLB_PROB.UA1 \ (GLB_PROB.LA1 \ (GLB_PROB.PA1 * x(nu+nz+1:end)))) );
    y(nu+nz+1:end)  = GLB_PROB.PA1' * (GLB_PROB.LA1' \ (GLB_PROB.UA1' \ (GLB_PROB.QA1' * tmp)));
  else
    y(1:nu)         = (C+perturb*M) \ x(1:nu);
    y(nu+1:nu+nz)   = (1/beta) * (R \ x(nu+1:nu+nz));
    y(nu+nz+1:end)  = A2 \ ((C+perturb*M) * (A1 \ x(nu+nz+1:end)));
  end

  GLB_PROB.ct = GLB_PROB.ct + 1;

end % function prdw
