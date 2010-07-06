
#include "NOX_Epetra_LinearSystem_SGJacobi.hpp"

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "EpetraExt_BlockVector.h"

NOX::Epetra::LinearSystemSGJacobi::
LinearSystemSGJacobi(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams_, 
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& detsolve_,
  const Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
  const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
  const Teuchos::RCP<Epetra_Operator>& J,
  const NOX::Epetra::Vector& cloneVector,
  const NOX::Epetra::Vector& detcloneVector,
  const Teuchos::RCP<NOX::Epetra::Scaling> s):
  detsolve(detsolve_),
  Cijk(Cijk_),
  jacInterfacePtr(iJac),
  jacPtr(J),
  leftHandSide(Teuchos::rcp(new Epetra_Vector(cloneVector.getEpetraVector()))),
  rightHandSide(Teuchos::rcp(new Epetra_Vector(cloneVector.getEpetraVector()))),
  detvec(Teuchos::rcp(new Epetra_Vector(detcloneVector.getEpetraVector()))),
  scaling(s),
  timer(cloneVector.getEpetraVector().Comm()),
  utils(printingParams)
{
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(jacPtr, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
}

NOX::Epetra::LinearSystemSGJacobi::
~LinearSystemSGJacobi()
{
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyJacobian(const NOX::Epetra::Vector& input, 
      		     NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(true);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  jacPtr->SetUseTranspose(false);

  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyJacobianInverse(Teuchos::ParameterList &params, 
      			    const NOX::Epetra::Vector &input, 
      			    NOX::Epetra::Vector &result)
{
  double startTime = timer.WallTime();
 
    int max_iter = params.get("Max Iterations",100 );
    double sg_tol = params.get("Tolerance", 1e-12);

//  *rightHandSide = input.getEpetraVector();

   //const Stokhos::VectorOrthogPoly<Epetra_Operator>& sg_J_poly_ref = *sg_J_poly;
 
    Teuchos::RCP<Epetra_Vector> sg_y =
     Teuchos::rcp(new Epetra_Vector(rightHandSide->Map()));

    Teuchos::RCP<Epetra_Vector> sg_df =
     Teuchos::rcp(new Epetra_Vector(rightHandSide->Map()));

    std::vector< Teuchos::RCP< Epetra_Vector> > sg_dx_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_f_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_df_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_kx_vec_all ;

    // Extract blocks
    EpetraExt::BlockVector sg_dx_block(View, detvec->Map(), result.getEpetraVector());
    EpetraExt::BlockVector sg_f_block(View, detvec->Map(), input.getEpetraVector());
    EpetraExt::BlockVector sg_df_block(View, detvec->Map(), *sg_df);

    int sz = sg_J_poly->basis()->size();
    int num_KL = sg_J_poly->basis()->dimension();
    const Teuchos::Array<double>& norms = sg_J_poly->basis()->norm_squared();
    for (int i=0; i<sz; i++) {
      sg_dx_vec_all.push_back(sg_dx_block.GetBlock(i));
      sg_f_vec_all.push_back(sg_f_block.GetBlock(i));
      sg_df_vec_all.push_back(sg_df_block.GetBlock(i));
      sg_kx_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(detvec->Map())));
    }

//  Teuchos::RCP< Epetra_Vector> kx ;
    Teuchos::RCP<Epetra_Vector> k0dx =
      Teuchos::rcp(new Epetra_Vector(detvec->Map()));
    Teuchos::RCP<Epetra_Vector> dx =
      Teuchos::rcp(new Epetra_Vector(detvec->Map()));
    Teuchos::RCP<Epetra_Vector> df =
      Teuchos::rcp(new Epetra_Vector(detvec->Map()));
 
    // Print initial residual norm
    double norm_f,norm_df;
    sg_f_block.Norm2(&norm_f);
    stokhos_op->Apply(sg_dx_block,*(sg_y));
    sg_df->Update(1.0,*sg_y,-1.0,sg_f_block,0.0);
    sg_df->Norm2(&norm_df);

    std::cout << "\nInitial residual norm = " << norm_f << std::endl;

    detsolve->setJacobianOperatorForSolve(sg_J_poly->getCoeffPtr(0));
int iter = 0;
//for (int iter=0;iter<2;iter++){
while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++;
     // Extract blocks

     // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
  //  ordinal_type Cijk_size = Cijk.size();
    for (int i=0; i<sz; i++) {
      (sg_df_vec_all[i])->Update(1.0, *sg_f_vec_all[i], 0.0);
    }
    for (int k=1; k<num_KL+1; k++) {
      int nj = Cijk->num_j(k);
      const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
      for (int jj=0; jj<nj; jj++) {
        int j = j_indices[jj];
        (*sg_J_poly)[k].Apply(*(sg_dx_vec_all[j]),*(sg_kx_vec_all[j]));
      }
      for (int jj=0; jj<nj; jj++) {
        int j = j_indices[jj];
        const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
        const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
        int ni = i_indices.size();
        for (int ii=0; ii<ni; ii++) {
          int i = i_indices[ii];
          double c = cijk_values[ii];  // C(i,j,k)
          sg_df_vec_all[i]->Update(-1.0*c/norms[i],*(sg_kx_vec_all[j]),1.0);
        }
      }
    } //End of k loop

    for(int i=0; i<sz; i++) {
      NOX::Epetra::Vector nox_df(sg_df_vec_all[i], NOX::Epetra::Vector::CreateView);
      NOX::Epetra::Vector nox_dx(sg_dx_vec_all[i], NOX::Epetra::Vector::CreateView);

     (*sg_J_poly)[0].Apply(*(sg_dx_vec_all[i]),*(sg_kx_vec_all[i]));
  //    nox_df.print(std::cout);
      nox_dx.init(0.0);
      // Solve linear system
      {
       TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");
       detsolve->applyJacobianInverse(params.sublist("Deterministic Krylov Solver"), nox_df, nox_dx);
      }
     (sg_df_vec_all[i])->Update(-1.0,*(sg_kx_vec_all[i]),1.0);
    }
//    stokhos_op->Apply(sg_dx_block,*(sg_y));
//    sg_df->Update(1.0,*sg_y,-1.0,sg_f_block,0.0);
    sg_df->Norm2(&norm_df);
    std::cout << "rel residual norm at iteration "<< iter <<" is " << norm_df/norm_f << std::endl;
  } //End of iter loop 

  //result.getEpetraVector() = *leftHandSide;

  double endTime = timer.WallTime();
  if (utils.isPrintType(Utils::LinearSolverDetails))
    utils.out() << "\n       Time required for one linear solve : " 
         << (endTime - startTime) << " (sec.)" << endl;;

  //return status;
  return true;
}

bool NOX::Epetra::LinearSystemSGJacobi::
applyRightPreconditioning(bool useTranspose,
      			      Teuchos::ParameterList& params, 
      			      const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  return false;
}

Teuchos::RCP<NOX::Epetra::Scaling> NOX::Epetra::LinearSystemSGJacobi::
getScaling()
{
  return scaling;
}

void NOX::Epetra::LinearSystemSGJacobi::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
  scaling = s;
  return;
}

bool NOX::Epetra::LinearSystemSGJacobi::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *jacPtr);
  sg_J_poly = stokhos_op->getOperatorBlocks();
  detsolve->setJacobianOperatorForSolve(sg_J_poly->getCoeffPtr(0));
  return success;
}

bool NOX::Epetra::LinearSystemSGJacobi::
createPreconditioner(const NOX::Epetra::Vector& x, 
      			    Teuchos::ParameterList& p,
      			    bool recomputeGraph) const
{
  std::cout << "calling createPreconditioner..." << std::endl;
  EpetraExt::BlockVector sg_x_block(View, detvec->Map(), x.getEpetraVector());
  bool success = detsolve->createPreconditioner(*(sg_x_block.GetBlock(0)), 
                                                p.sublist("Deterministic Krylov Solver"), 
                                                recomputeGraph);
  return success;
}

bool NOX::Epetra::LinearSystemSGJacobi::
destroyPreconditioner() const
{
  return detsolve->destroyPreconditioner();
}

bool NOX::Epetra::LinearSystemSGJacobi::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
      		Teuchos::ParameterList& linearSolverParams) const
{  
  EpetraExt::BlockVector sg_x_block(View, detvec->Map(), x.getEpetraVector());
  bool success = detsolve->recomputePreconditioner(*(sg_x_block.GetBlock(0)), 
                                                   linearSolverParams.sublist("Deterministic Krylov Solver"));
  return success;

}

NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemSGJacobi::
getPreconditionerPolicy(bool advanceReuseCounter)
{
  return detsolve->getPreconditionerPolicy(advanceReuseCounter);
} 

bool NOX::Epetra::LinearSystemSGJacobi::
isPreconditionerConstructed() const
{
  return detsolve->isPreconditionerConstructed();
}

bool NOX::Epetra::LinearSystemSGJacobi::
hasPreconditioner() const
{
  return detsolve->hasPreconditioner();
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getJacobianOperator() const
{
  return jacPtr;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getJacobianOperator()
{
  return jacPtr;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGJacobi::
getGeneratedPrecOperator()
{
  return Teuchos::null;
}

void NOX::Epetra::LinearSystemSGJacobi::
setJacobianOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(jacPtr, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
  return;
}

void NOX::Epetra::LinearSystemSGJacobi::
setPrecOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  return;
}
