function [y] = imshift( x, usr_par)
%
%   compute and/or apply a complex-shift KKT preconditioner
%

  global GLB_PROB;

  nu = GLB_PROB.nu;
  nz = GLB_PROB.nz;

  beta = GLB_PROB.beta;
  r = sqrt(-1)/sqrt(beta);
  K = GLB_PROB.K;
  M = GLB_PROB.M;
  C = GLB_PROB.C;
  R = GLB_PROB.R;
  J = GLB_PROB.J;
  A1 = K - r*M;

  A2 = ctranspose(A1);

  y = zeros(size(x));

  if (GLB_PROB.reuse_factors)
    if (~GLB_PROB.factors_computed)
      fprintf("  Computing factors (to be reused) ...\n  ");
      tic
      [GLB_PROB.LC, GLB_PROB.UC, GLB_PROB.PC, GLB_PROB.QC]     = lu(C);
      [GLB_PROB.LR, GLB_PROB.UR, GLB_PROB.PR, GLB_PROB.QR]     = lu(R);
      [GLB_PROB.LA1, GLB_PROB.UA1, GLB_PROB.PA1, GLB_PROB.QA1] = lu(A1);
      %[GLB_PROB.LA2, GLB_PROB.UA2, GLB_PROB.PA2, GLB_PROB.QA2] = lu(A2);
      toc
      GLB_PROB.factors_computed = true;
    end
    y(1:nu)         = GLB_PROB.QC * (GLB_PROB.UC \ (GLB_PROB.LC \ (GLB_PROB.PC * x(1:nu))));
    y(nu+1:nu+nz)   = (1/beta) * ( GLB_PROB.QR * (GLB_PROB.UR \ (GLB_PROB.LR \ (GLB_PROB.PR * x(nu+1:nu+nz)))) );
    tmp             = M * ( GLB_PROB.QA1 * (GLB_PROB.UA1 \ (GLB_PROB.LA1 \ (GLB_PROB.PA1 * x(nu+nz+1:end)))) );
    y(nu+nz+1:end)  = GLB_PROB.PA1' * (GLB_PROB.LA1' \ (GLB_PROB.UA1' \ (GLB_PROB.QA1' * tmp)));
    %y(nu+nz+1:end)  = GLB_PROB.QA2 * (GLB_PROB.UA2 \ (GLB_PROB.LA2 \ (GLB_PROB.PA2 * tmp)));
  else
    y(1:nu)         = C \ x(1:nu);
    y(nu+1:nu+nz)   = (1/beta) * (R \ x(nu+1:nu+nz));
    y(nu+nz+1:end)  = A2 \ (M * (A1 \ x(nu+nz+1:end)));
  end


  GLB_PROB.ct = GLB_PROB.ct + 1;

end % function imshift
