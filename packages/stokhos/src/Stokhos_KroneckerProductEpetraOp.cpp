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
#include "Stokhos_KroneckerProductEpetraOp.hpp"

#include "Ifpack.h"
//#include <Stokhos_IfpackPreconditionerFactory.hpp>

Stokhos::KroneckerProductEpetraOp::KroneckerProductEpetraOp(
   const Teuchos::RCP<const Epetra_Map>& base_map_,
   const Teuchos::RCP<const Epetra_Map>& sg_map_,
   unsigned int num_blocks_,
   const Teuchos::RCP<Epetra_Operator>& mean_op_, 
   const Teuchos::RCP<const Stokhos::Sparse3Tensor<int,double> >& Cijk_,
   const Teuchos::RCP<Epetra_Operator>& J):
  label("Stokhos Kronecker Product Preconditioner"),
  base_map(base_map_),
  sg_map(sg_map_),
  mean_op(mean_op_),
  useTranspose(false),
  num_blocks(num_blocks_),
  Cijk(Cijk_)
{
  stokhos_op = Teuchos::rcp_dynamic_cast<Stokhos::MatrixFreeEpetraOp>(J, true);
  sg_J_poly = stokhos_op->getOperatorBlocks();
  //Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis =
    // sg_J_poly->basis();
  basis = sg_J_poly->basis();
  // Number of stochastic rows
  int num_rows = basis->size();

  //Teuchos::RCP<Epetra_CrsMatrix> matA0 =
    ///  Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(0), true);
  matA0 =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(0), true);

  // Replicated local map
  Epetra_LocalMap G_map(num_rows, 0, matA0->Comm());

  // Graph to be created
  Teuchos::RCP<Epetra_CrsGraph> graph =
    Teuchos::rcp(new Epetra_CrsGraph(Copy, G_map, 0));

  // Loop over Cijk entries including a non-zero in the graph at
   // indices (i,j) if there is any k for which Cijk is non-zero
//  int Cijk_size = Cijk->size();
  //for (int k=0; k<Cijk_size; k++) {
  for (int k=0; k<sg_J_poly->size(); k++) { 
    int nj = Cijk->num_j(k);
    const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
    for (int jj=0; jj<nj; jj++) {
      int j = j_indices[jj];
      const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
      int ni = i_indices.size();
      for (int ii=0; ii<ni; ii++) {
        int i = i_indices[ii];
        graph->InsertGlobalIndices(i, 1, &j);
      }
    }
  }

  // Sort, remove redundencies, transform to local, ...
  graph->FillComplete();

  // Construct G matrices, G_{ij} = \sum tr(A'B)/tr(A'A)*<psi_alpha,psi_i,psi_j>.
  //Teuchos::RCP<Epetra_CrsMatrix> G = 
    //  Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
  G = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));

  if (mean_op != Teuchos::null)
    label = std::string("Stokhos Kronecker Product Preconditioner:\n") + 
      std::string("		***** ") + 
      std::string(mean_op->Label());
}

Stokhos::KroneckerProductEpetraOp::~KroneckerProductEpetraOp()
{
}

void
Stokhos::KroneckerProductEpetraOp::setMeanOperator(const Teuchos::RCP<Epetra_Operator>& op)
{
  mean_op = op;  
  label = std::string("Stokhos Kronecker Product Preconditioner:\n") + 
    std::string("		***** ") + 
    std::string(mean_op->Label());
  
  int num_rows = basis->size();
  // Replicated local map
   Epetra_LocalMap G_map(num_rows, 0, matA0->Comm());

   double traceAB0 = MatrixTrace(matA0, matA0);
   int NumMyElements = G_map.NumMyElements();
   int * MyGlobalElements = G_map.MyGlobalElements();
   //for (int k=0; k<Cijk_size; k++) {
   for (int k=0; k<sg_J_poly->size(); k++) {
     Teuchos::RCP<Epetra_CrsMatrix> matA_k =
       Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(k), true);
     double traceAB = MatrixTrace(matA_k, matA0);
     int nj = Cijk->num_j(k);
     const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
     for( int irow=0 ; irow<NumMyElements; ++irow ) {
     for (int jj=0; jj<nj; jj++) {
       int j = j_indices[jj];
       if (irow==j) { 
       const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
       const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
       int ni = i_indices.size();
       double *Values = new double[ni];
       int *Indices = new int[ni];
       for (int ii=0; ii<ni; ii++) {
         //int i = i_indices[ii];
         //double c = cijk_values[ii];  // C(i,j,k)
         Indices[ii] = i_indices[ii];
         Values[ii] = cijk_values[ii]*traceAB/traceAB0;  // C(i,j,k)
       }
       G->SumIntoGlobalValues(MyGlobalElements[irow], ni, Values, Indices); 
     }       
     }
     }
   }
 G->FillComplete();
 // ************Construct preconditioner*************
  //
 // ParameterList ifpackList;
  Teuchos::RCP<Teuchos::ParameterList> ifpackList =
      Teuchos::rcp(new Teuchos::ParameterList);

  // allocates an IFPACK factory. No data is associated
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  std::string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.

  //Teuchos::RCP<Ifpack_Preconditioner> Prec = 
    // Teuchos::rcp( Factory.Create(PrecType, &*G, OverlapLevel) );
  Prec = Teuchos::rcp( Factory.Create(PrecType, &*G, OverlapLevel) );

  assert(Prec != Teuchos::null);
  // specify parameters for ILU
  ifpackList->set("fact: drop tolerance", 1e-9);
  ifpackList->set("fact: ilut level-of-fill", 1.0);
  // the combine mode is on the following:
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  // Their meaning is as defined in file Epetra_CombineMode.h
  ifpackList->set("schwarz: combine mode", "Add");
  // sets the parameters
  //IFPACK_CHK_ERR(Prec->SetParameters(*ifpackList));
  Prec->SetParameters(*ifpackList);

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  //IFPACK_CHK_ERR(Prec->Initialize());
  Prec->Initialize();

  // Builds the preconditioners, by looking for the values of
  // the matrix.
  //IFPACK_CHK_ERR(Prec->Compute());
  Prec->Compute();

}

Teuchos::RCP<Epetra_Operator>
Stokhos::KroneckerProductEpetraOp::getMeanOperator()
{
  return mean_op;
}

int 
Stokhos::KroneckerProductEpetraOp::SetUseTranspose(bool UseTranspose) 
{
  useTranspose = UseTranspose;
  mean_op->SetUseTranspose(useTranspose);

  return 0;
}

int 
Stokhos::KroneckerProductEpetraOp::Apply(const Epetra_MultiVector& Input, 
			     Epetra_MultiVector& Result) const
{
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (unsigned int i=0; i<num_blocks; i++) {
    mean_op->Apply(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

  return 0;
}

int 
Stokhos::KroneckerProductEpetraOp::ApplyInverse(const Epetra_MultiVector& Input, 
				    Epetra_MultiVector& Result) const
{
  TEUCHOS_FUNC_TIME_MONITOR("Total Kronecker Product Prec Time");
/*  // We have to be careful if Input and Result are the same vector.
  // If this is the case, the only possible solution is to make a copy
  const Epetra_MultiVector *input = &Input;
  bool made_copy = false;
  if (Input.Values() == Result.Values()) {
    input = new Epetra_MultiVector(Input);
    made_copy = true;
  }
*/
  EpetraExt::BlockMultiVector sg_input(View, *base_map, Input);
  EpetraExt::BlockMultiVector sg_result(View, *base_map, Result);
  for (unsigned int i=0; i<num_blocks; i++) {
    mean_op->ApplyInverse(*(sg_input.GetBlock(i)), *(sg_result.GetBlock(i)));
  }

/*  Teuchos::RCP<Epetra_CrsMatrix> matA0 =
      Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(0), true);
  double traceAB0 = MatrixTrace(matA0, matA0);
  //std::cout << "print traceAB0 = " << traceAB0 << std::endl;
  //std::cout << "print matA0 = " << *matA0 << std::endl;
  //Epetra_Comm comm = matA0->Comm();  
  Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> > basis =
    sg_J_poly->basis();
//  Teuchos::RCP<Epetra_CrsGraph> graph = Stokhos::sparse3Tensor2CrsGraph(basis,Cijk,comm);
  // Number of stochastic rows
    int num_rows = basis->size();

    // Replicated local map
    Epetra_LocalMap G_map(num_rows, 0, matA0->Comm());

    // Graph to be created
    Teuchos::RCP<Epetra_CrsGraph> graph =
      Teuchos::rcp(new Epetra_CrsGraph(Copy, G_map, 0));

   // Loop over Cijk entries including a non-zero in the graph at
    // indices (i,j) if there is any k for which Cijk is non-zero
    int Cijk_size = Cijk->size();
    //for (int k=0; k<Cijk_size; k++) {
    for (int k=0; k<sg_J_poly->size(); k++) {
      int nj = Cijk->num_j(k);
      const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
      for (int jj=0; jj<nj; jj++) {
        int j = j_indices[jj];
        const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
        int ni = i_indices.size();
        for (int ii=0; ii<ni; ii++) {
          int i = i_indices[ii];
          graph->InsertGlobalIndices(i, 1, &j);
        }
      }
    }

    // Sort, remove redundencies, transform to local, ...
    graph->FillComplete();

   // Construct G matrices, G_alpha_{ij} = <psi_alpha,psi_i,psi_j>.
   Teuchos::RCP<Epetra_CrsMatrix> G = 
      Teuchos::rcp(new Epetra_CrsMatrix(Copy, *graph));
   //double *Values = new double[basis->size()];
  // int *Indices = new int[basis->size()];
   int NumMyElements = G_map.NumMyElements();
   int * MyGlobalElements = G_map.MyGlobalElements();
   //for (int k=0; k<Cijk_size; k++) {
   for (int k=0; k<sg_J_poly->size(); k++) {
     Teuchos::RCP<Epetra_CrsMatrix> matA_k =
       Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>((*sg_J_poly).getCoeffPtr(k), true);
     double traceAB = MatrixTrace(matA_k, matA0);
     int nj = Cijk->num_j(k);
     const Teuchos::Array<int>& j_indices = Cijk->Jindices(k);
     for( int irow=0 ; irow<NumMyElements; ++irow ) {
     for (int jj=0; jj<nj; jj++) {
       int j = j_indices[jj];
       if (irow==j) { 
       const Teuchos::Array<double>& cijk_values = Cijk->values(k,jj);
       const Teuchos::Array<int>& i_indices = Cijk->Iindices(k,jj);
       int ni = i_indices.size();
       double *Values = new double[ni];
       int *Indices = new int[ni];
       for (int ii=0; ii<ni; ii++) {
         //int i = i_indices[ii];
         //double c = cijk_values[ii];  // C(i,j,k)
         Indices[ii] = i_indices[ii];
         Values[ii] = cijk_values[ii]*traceAB/traceAB0;  // C(i,j,k)
       }
       G->SumIntoGlobalValues(MyGlobalElements[irow], ni, Values, Indices); 
     }       
     }
     }
   }
 G->FillComplete();
 //std::cout << "print G" << *G << std::endl;
*/  //std::cout << "sg_result(990) = "<< (*sg_result.GetBlock(0))[0][100] << std::endl;
/* // Compute ILU preconditioner of G
 Teuchos::RCP<Teuchos::ParameterList> precParams =
      Teuchos::rcp(new Teuchos::ParameterList);
 precParams.set("Ifpack Preconditioner", "ILU");
 precParams.set("Overlap", 0);
// Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType, &*G, OverlapLevel) );
 */
/*
 // ************Construct preconditioner*************
  //
//  ParameterList ifpackList;
  Teuchos::RCP<Teuchos::ParameterList> ifpackList =
      Teuchos::rcp(new Teuchos::ParameterList);

  // allocates an IFPACK factory. No data is associated
  // to this object (only method Create()).
  Ifpack Factory;

  // create the preconditioner. For valid PrecType values,
  // please check the documentation
  std::string PrecType = "ILU"; // incomplete LU
  int OverlapLevel = 1; // must be >= 0. If Comm.NumProc() == 1,
                        // it is ignored.

  Teuchos::RCP<Ifpack_Preconditioner> Prec = Teuchos::rcp( Factory.Create(PrecType, &*G, OverlapLevel) );

  assert(Prec != Teuchos::null);
  // specify parameters for ILU
  ifpackList->set("fact: drop tolerance", 1e-9);
  ifpackList->set("fact: ilut level-of-fill", 1.0);
  // the combine mode is on the following:
  // "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
  // Their meaning is as defined in file Epetra_CombineMode.h
  ifpackList->set("schwarz: combine mode", "Add");
  // sets the parameters
  IFPACK_CHK_ERR(Prec->SetParameters(*ifpackList));

  // initialize the preconditioner. At this point the matrix must
  // have been FillComplete()'d, but actual values are ignored.
  IFPACK_CHK_ERR(Prec->Initialize());

  // Builds the preconditioners, by looking for the values of
  // the matrix.
  IFPACK_CHK_ERR(Prec->Compute());
*/
  //int vecLen = sg_input.GetBlock(0)->GlobalLength(); // Global length of vector.
  // Replicated local map
  int num_rows = basis->size();
  Epetra_LocalMap G_map(num_rows, 0, matA0->Comm());
  int NumMyElements = G_map.NumMyElements();
  int * MyGlobalElements = G_map.MyGlobalElements();

  int vecLen = sg_input.GetBlock(0)->MyLength(); // Global length of vector.
  Teuchos::RCP<Epetra_MultiVector> result_MVT
     = Teuchos::rcp(new Epetra_MultiVector(G_map, vecLen));
  for( int irow=0 ; irow<NumMyElements; ++irow ) {
    for (int icol=0;icol<vecLen;icol++) {
      result_MVT->ReplaceGlobalValue (MyGlobalElements[irow], icol, (*sg_result.GetBlock(MyGlobalElements[irow]))[0][icol]);
    }
  }
  Teuchos::RCP<Epetra_MultiVector> result_MVT_new
     = Teuchos::rcp(new Epetra_MultiVector(G_map, vecLen));
  Prec->ApplyInverse(*result_MVT, *result_MVT_new);
 
  for( int irow=0 ; irow<NumMyElements; ++irow ) {
    for (int icol=0;icol<vecLen;icol++) {
      (sg_result.GetBlock(MyGlobalElements[irow]))->ReplaceGlobalValue (icol, 0, (*result_MVT_new)[icol][MyGlobalElements[irow]]);
    }
  }
//  std::cout << "print multi_vec = " << (*result_MVT_new)[100][100] << endl;
  return 0;
}

double 
Stokhos::KroneckerProductEpetraOp::NormInf() const
{
  return mean_op->NormInf();
}


const char* 
Stokhos::KroneckerProductEpetraOp::Label () const
{
  return const_cast<char*>(label.c_str());
}
  
bool 
Stokhos::KroneckerProductEpetraOp::UseTranspose() const
{
  return useTranspose;
}

bool 
Stokhos::KroneckerProductEpetraOp::HasNormInf() const
{
  return true;
}

const Epetra_Comm & 
Stokhos::KroneckerProductEpetraOp::Comm() const
{
  return base_map->Comm();
}
const Epetra_Map& 
Stokhos::KroneckerProductEpetraOp::OperatorDomainMap() const
{
  return *sg_map;
}

const Epetra_Map& 
Stokhos::KroneckerProductEpetraOp::OperatorRangeMap() const
{
  return *sg_map;
}

//double Stokhos::KroneckerProductEpetraOp::MatrixTrace(Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::RCP<Epetra_CrsMatrix> B) const {
//double Stokhos::KroneckerProductEpetraOp::MatrixTrace(Epetra_Operator& A,Epetra_Operator& B) const {
double Stokhos::KroneckerProductEpetraOp::MatrixTrace(Teuchos::RCP<Epetra_CrsMatrix>& matA,Teuchos::RCP<Epetra_CrsMatrix>& matB) const {
//Teuchos::RCP<Epetra_CrsMatrix> matA =
  //    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A, true);
//Teuchos::RCP<Epetra_CrsMatrix> matB =
  //    Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(B, true);
int n = matA->NumMyRows(); // # of rows on the processor.
double traceAB = 0;
for (int i=0;i<n;i++) {
  int m = matA->NumMyEntries(i); // # of non zero entries in row i.
  for (int j=0;j<m;j++) { 
    traceAB += ((*matA)[i][j])*((*matB)[i][j]);
  }
}
return traceAB;
}

