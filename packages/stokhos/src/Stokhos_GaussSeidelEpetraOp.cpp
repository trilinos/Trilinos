// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Epetra_config.h"
#include "EpetraExt_BlockMultiVector.h"
#include "Stokhos_GaussSeidelEpetraOp.hpp"
#include "EpetraExt_BlockVector.h"

Stokhos::GaussSeidelEpetraOp::GaussSeidelEpetraOp(
  const Teuchos::RCP<const Epetra_Map>& base_map_,
  const Teuchos::RCP<const Epetra_Map>& sg_map_,
  unsigned int num_blocks_,
  Teuchos::ParameterList& linearSolverParams, 
  const Teuchos::RCP<NOX::Epetra::LinearSystem>& detsolve_,
  const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
  const Teuchos::RCP<Epetra_Operator>& J):
  label("Stokhos Gauss-Seidel Preconditioner"),
  base_map(base_map_),
  sg_map(sg_map_),
  useTranspose(false),
  num_blocks(num_blocks_),
  detsolve(detsolve_),
  Cijk(Cijk_),
  jacPtr(J),
  detvec(Teuchos::rcp(new Epetra_Vector(*base_map))),
  params(linearSolverParams)
{
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(jacPtr, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
}

Stokhos::GaussSeidelEpetraOp::~GaussSeidelEpetraOp()
{
}

void
Stokhos::GaussSeidelEpetraOp::setOperatorAndConstructPreconditioner(
  const Teuchos::RCP<Epetra_Operator>& J, const Epetra_Vector& x) 
{
  jacPtr = J;
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(jacPtr, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
  EpetraExt::BlockVector sg_x_block(View, detvec->Map(), x);
  detsolve->setJacobianOperatorForSolve(sg_J_poly->getCoeffPtr(0));
  detsolve->createPreconditioner(*(sg_x_block.GetBlock(0)),
				 params.sublist("Deterministic Krylov Solver"),
				 false);
}

int 
Stokhos::GaussSeidelEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;

  return 0;
}

int 
Stokhos::GaussSeidelEpetraOp::Apply(const Epetra_MultiVector& Input, 
			     Epetra_MultiVector& Result) const
{
  return stokhos_op->Apply(Input,Result);
}

int 
Stokhos::GaussSeidelEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
				    Epetra_MultiVector& Result) const
{
    int max_iter = params.get("Max Iterations",100 );
    double sg_tol = params.get("Tolerance", 1e-12);
    bool MatVecTable = params.get("Save MatVec Table", true);

    // We have to be careful if Input and Result are the same vector.
    // If this is the case, the only possible solution is to make a copy
    const Epetra_MultiVector *input = &Input;
    bool made_copy = false;
    if (Input.Values() == Result.Values()) {
      input = new Epetra_MultiVector(Input);
      made_copy = true;
    } 
 
    Teuchos::RCP<Epetra_Vector> sg_y =
     Teuchos::rcp(new Epetra_Vector(*sg_map));

    Teuchos::RCP<Epetra_Vector> sg_df =
     Teuchos::rcp(new Epetra_Vector(*sg_map));

    Teuchos::Array< Teuchos::RCP< Epetra_Vector> > sg_dx_vec_all ;
    Teuchos::Array< Teuchos::RCP< Epetra_Vector> > sg_dxold_vec_all ;
    Teuchos::Array< Teuchos::RCP< Epetra_Vector> > sg_f_vec_all ;
    Teuchos::Array< Teuchos::RCP< Epetra_Vector> > sg_df_vec_all ;
    Teuchos::Array< Teuchos::RCP< Epetra_Vector> > sg_kx_vec_all ;

    // Extract blocks
    EpetraExt::BlockVector sg_dx_block(View, detvec->Map(), *(Result(0)));
    EpetraExt::BlockVector sg_f_block(View, detvec->Map(), *((*input)(0)));
    EpetraExt::BlockVector sg_df_block(View, detvec->Map(), *sg_df);
    sg_dx_block.PutScalar(0.0);

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

//    std::cout << "\nInitial residual norm = " << norm_f << std::endl;

    // Set deterministic solver K0
    //detsolve->setJacobianOperatorForSolve(sg_J_poly->getCoeffPtr(0));

// 2D array to store mat-vec results
//int expn;
std::vector< std::vector< Teuchos::RCP< Epetra_Vector> > > Kx_table;
if (MatVecTable) {
 Kx_table.resize(sz);
 for (int i=0;i<sz;i++) {
   Kx_table[i].resize(sg_J_poly->size());
   for (int j=0;j<sg_J_poly->size();j++) {
     Kx_table[i][j] = Teuchos::rcp(new Epetra_Vector(detvec->Map()));
   }
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
        if (i!=0 && i<=num_KL+1) {
         if (MatVecTable) {
          sg_df_vec_all[k]->Update(-1.0*c/norms[k],*(Kx_table[j][i]),1.0);      
         }
         else {
          (*sg_J_poly)[i].Apply(*(sg_dx_vec_all[j]),*kx);
          sg_df_vec_all[k]->Update(-1.0*c/norms[k],*kx,1.0);      
         }  
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
      // Teuchos::RCP<const Epetra_Operator> M = detsolve->getGeneratedPrecOperator();
      // M->ApplyInverse(*(sg_df_vec_all[k]), *(sg_dx_vec_all[k]));
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
//    std::cout << "rel residual norm at iteration "<< iter <<" is " << norm_df/norm_f << std::endl;
  } //End of iter loop 

  //result.getEpetraVector() = *leftHandSide;


  if (made_copy)
    delete input;

  //return status;
  return 0; 
}

double 
Stokhos::GaussSeidelEpetraOp::NormInf() const
{
  return stokhos_op->NormInf();
}


const char* 
Stokhos::GaussSeidelEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::GaussSeidelEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::GaussSeidelEpetraOp::HasNormInf() const
{
  return stokhos_op->HasNormInf();
}

const Epetra_Comm & 
Stokhos::GaussSeidelEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::GaussSeidelEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::GaussSeidelEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}
