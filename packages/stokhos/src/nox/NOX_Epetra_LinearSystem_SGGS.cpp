
#include "NOX_Epetra_LinearSystem_SGGS.hpp"

#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "EpetraExt_BlockVector.h"
//bool full_expansion = true;

NOX::Epetra::LinearSystemSGGS::
LinearSystemSGGS(
  Teuchos::ParameterList& printingParams, 
  Teuchos::ParameterList& linearSolverParams, 
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

NOX::Epetra::LinearSystemSGGS::
~LinearSystemSGGS()
{
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobian(const NOX::Epetra::Vector& input, 
      		     NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobianTranspose(const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(true);
  int status = jacPtr->Apply(input.getEpetraVector(), 
				  result.getEpetraVector());
  jacPtr->SetUseTranspose(false);

  return (status == 0);
}

bool NOX::Epetra::LinearSystemSGGS::
applyJacobianInverse(Teuchos::ParameterList &params, 
      			    const NOX::Epetra::Vector &input, 
      			    NOX::Epetra::Vector &result)
{
  double startTime = timer.WallTime();

    int max_iter = params.get("Max Iterations",100 );
    double sg_tol = params.get("Tolerance", 1e-12);
    bool MatVecTable = params.get("Save MatVec Table", true);

//  *rightHandSide = input.getEpetraVector();

   //const Stokhos::VectorOrthogPoly<Epetra_Operator>& sg_J_poly_ref = *sg_J_poly;
 
    Teuchos::RCP<Epetra_Vector> sg_y =
     Teuchos::rcp(new Epetra_Vector(rightHandSide->Map()));

    Teuchos::RCP<Epetra_Vector> sg_df =
     Teuchos::rcp(new Epetra_Vector(rightHandSide->Map()));

    std::vector< Teuchos::RCP< Epetra_Vector> > sg_dx_vec_all ;
    std::vector< Teuchos::RCP< Epetra_Vector> > sg_dxold_vec_all ;
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
//      sg_df_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(detvec->Map())));
      sg_kx_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(detvec->Map())));
      sg_dxold_vec_all.push_back(Teuchos::rcp(new Epetra_Vector(detvec->Map())));
    }

//  Teuchos::RCP< Epetra_Vector> kx ;
    Teuchos::RCP<Epetra_Vector> kx =
      Teuchos::rcp(new Epetra_Vector(detvec->Map()));
    Teuchos::RCP<Epetra_Vector> dx =
      Teuchos::rcp(new Epetra_Vector(detvec->Map()));
    Teuchos::RCP<Epetra_Vector> df =
      Teuchos::rcp(new Epetra_Vector(detvec->Map()));
 
    // Print initial residual norm
    double norm_f,norm_df;
    std::vector<double> norm_df_k,norm_f_k;
    norm_df_k.resize(sz); 
    norm_f_k.resize(sz); 
    for (int i=0;i<sz;i++) {
      norm_df_k[i] = 1.0;
      norm_f_k[i] = 1.0;
    }
    sg_f_block.Norm2(&norm_f);
    stokhos_op->Apply(sg_dx_block,*(sg_y));
    sg_df->Update(1.0,*sg_y,-1.0,sg_f_block,0.0);
    sg_df->Norm2(&norm_df);

    std::cout << "\nInitial residual norm = " << norm_f << std::endl;

    // Set deterministic solver K0
    detsolve->setJacobianOperatorForSolve(sg_J_poly->getCoeffPtr(0));

// 2D array to store mat-vec results
//int expn;
std::vector< std::vector< Teuchos::RCP< Epetra_Vector> > > Kx_table(sz);
for (int i=0;i<sz;i++) {
/*  if(full_expansion){
  Kx_table[i].resize(sz);
  expn =sz;
  }
//  if(!full_expansion){
  else {
  Kx_table[i].resize(num_KL+1);
  expn = num_KL+1;
  }*/
  Kx_table[i].resize(num_KL+1);
  for (int j=0;j<num_KL+1;j++) {
    Kx_table[i][j] = Teuchos::rcp(new Epetra_Vector(detvec->Map()));
  }
}

int iter = 0;
//for (int iter=0;iter<28;iter++){
while (((norm_df/norm_f)>sg_tol) && (iter<max_iter)) {
    TEUCHOS_FUNC_TIME_MONITOR("Total global solve Time");
    iter++;

/*       for(int i=0; i<sz; i++) {
         (*sg_J_poly)[0].Apply(*(sg_dx_vec_all[i]),*(sg_kx_vec_all[i]));
       }
*/
    // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
  //  ordinal_type Cijk_size = Cijk.size();
    for (int k=0; k<sz; k++) {
//      if ((norm_df_k[k]/norm_f_k[k])>1e-12) {
//      df->Update(1.0, *sg_f_vec_all[k], 0.0);
      (sg_df_vec_all[k])->Update(1.0, *sg_f_vec_all[k], 0.0);
      int nl = Cijk->num_values(k);
      for (int l=0; l<nl; l++) {
        int i,j;
        double c;
        Cijk->value(k,l,i,j,c); 
        if (i!=0) {
         if (!MatVecTable) {
          (*sg_J_poly)[i].Apply(*(sg_dx_vec_all[j]),*(Kx_table[j][i]));
         }
         sg_df_vec_all[k]->Update(-1.0*c/norms[k],*(Kx_table[j][i]),1.0);      
        }
      }

      NOX::Epetra::Vector nox_df(sg_df_vec_all[k], NOX::Epetra::Vector::CreateView);
      NOX::Epetra::Vector nox_dx(sg_dx_vec_all[k], NOX::Epetra::Vector::CreateView);

      (*sg_J_poly)[0].Apply(*(sg_dx_vec_all[k]),*(sg_kx_vec_all[k]));
      //std::cout << "nox_df" << nox_df << std::endl;
//      nox_df.print(std::cout);
      nox_dx.init(0.0);
      // Solve linear system
     {
      TEUCHOS_FUNC_TIME_MONITOR("Total deterministic solve Time");
      detsolve->applyJacobianInverse(params.sublist("Deterministic Krylov Solver"), nox_df, nox_dx);
//     Teuchos::RCP<const Epetra_Operator> M = detsolve->getGeneratedPrecOperator();
//     M->ApplyInverse(*(sg_df_vec_all[k]), *(sg_dx_vec_all[k]));
     }
      (sg_dxold_vec_all[k])->Update(1.0, *sg_dx_vec_all[k], -1.0);
      (sg_dx_vec_all[k])->Norm2(&norm_f_k[k]);
      (sg_df_vec_all[k])->Update(-1.0,*(sg_kx_vec_all[k]),1.0);
//      (sg_df_vec_all[k])->Norm2(&norm_df_k[k]);
      (sg_dxold_vec_all[k])->Norm2(&norm_df_k[k]);
      (sg_dxold_vec_all[k])->Update(1.0, *sg_dx_vec_all[k], 0.0);
//      (sg_f_vec_all[0])->Norm2(&norm_f);

      if (MatVecTable) {
       for(int i=1;i<num_KL+1;i++) {
	(*sg_J_poly)[i].Apply(*(sg_dx_vec_all[k]),*(Kx_table[k][i]));
       }
      }
//     } //End of if((norm_df/norm_f)>1e-12) loop
    } //End of k loop

  /*  for(int i=0; i<sz; i++) {
      (sg_df_vec_all[i])->Update(-1.0,*(sg_kx_vec_all[i]),1.0);
    }
*/
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

bool NOX::Epetra::LinearSystemSGGS::
applyRightPreconditioning(bool useTranspose,
      			      Teuchos::ParameterList& params, 
      			      const NOX::Epetra::Vector& input, 
      			      NOX::Epetra::Vector& result) const
{
  return false;
}

Teuchos::RCP<NOX::Epetra::Scaling> NOX::Epetra::LinearSystemSGGS::
getScaling()
{
  return scaling;
}

void NOX::Epetra::LinearSystemSGGS::
resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
  scaling = s;
  return;
}

bool NOX::Epetra::LinearSystemSGGS::
computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), 
						  *jacPtr);
  sg_J_poly = stokhos_op->getOperatorBlocks();
  detsolve->setJacobianOperatorForSolve(sg_J_poly->getCoeffPtr(0));
  return success;
}

bool NOX::Epetra::LinearSystemSGGS::
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

bool NOX::Epetra::LinearSystemSGGS::
destroyPreconditioner() const
{
  return detsolve->destroyPreconditioner();
}

bool NOX::Epetra::LinearSystemSGGS::
recomputePreconditioner(const NOX::Epetra::Vector& x, 
      		Teuchos::ParameterList& linearSolverParams) const
{  
  EpetraExt::BlockVector sg_x_block(View, detvec->Map(), x.getEpetraVector());
  bool success = detsolve->recomputePreconditioner(*(sg_x_block.GetBlock(0)), 
                                                   linearSolverParams.sublist("Deterministic Krylov Solver"));
  return success;

}

NOX::Epetra::LinearSystem::PreconditionerReusePolicyType 
NOX::Epetra::LinearSystemSGGS::
getPreconditionerPolicy(bool advanceReuseCounter)
{
  return detsolve->getPreconditionerPolicy(advanceReuseCounter);
} 

bool NOX::Epetra::LinearSystemSGGS::
isPreconditionerConstructed() const
{
  return detsolve->isPreconditionerConstructed();
}

bool NOX::Epetra::LinearSystemSGGS::
hasPreconditioner() const
{
  return detsolve->hasPreconditioner();
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getJacobianOperator() const
{
  return jacPtr;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getJacobianOperator()
{
  return jacPtr;
}

Teuchos::RCP<const Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

Teuchos::RCP<Epetra_Operator> NOX::Epetra::LinearSystemSGGS::
getGeneratedPrecOperator()
{
  return Teuchos::null;
}

void NOX::Epetra::LinearSystemSGGS::
setJacobianOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(jacPtr, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
  return;
}

void NOX::Epetra::LinearSystemSGGS::
setPrecOperatorForSolve(const 
      	 Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  return;
}
