#ifndef __PB_EpetraLSCHelpers_hpp__
#define __PB_EpetraLSCHelpers_hpp__

// stl includes
#include <string>
#include <vector>

// Epetra includes
#include "Epetra_CrsMatrix.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

namespace PB {

namespace Epetra {

// Navier-Stokes specific support functions
//////////////////////////////////////////////////////////

// This function computes the constants needed for the 
// LSC-Stabilized preconditioner developed in 
//    Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
//    for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
//
// arguments (in):
//    A:  2x2 block matrix [ F Bt ; B C ] (see below for details)
//    Qu: Velocity mass matrix
// arguments (out):
//    gamma: gamma parameter (see EHSST2007 4.28)
//    alpha: alpha parameter (see EHSST2007 4.29)
//
void computeLSCConstants(const Epetra_Operator & A,const Epetra_RowMatrix & Qu,
                         double & gamma,double & alpha);

// This function computes the constants needed for the 
// LSC-Stabilized preconditioner developed in 
//
//    Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
//    for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
//
// Arguments that are Epetra_CrsMatrix require the specific functionality of these 
// data structures. (In particular EpetraExt::MatrixMatrix and diagonal extraction)
//
// arguments (in):
//    F:  Convection operator on discrete velocity
//    Bt: Gradient operator on discrete pressure
//    B:  Divergence operator on discrete velocity
//    C:  Stabilization operator
//    Qu: Velocity mass matrix
// arguments (out):
//    gamma: gamma parameter (see EHSST2007 4.28)
//    alpha: alpha parameter (see EHSST2007 4.29)
//
void computeLSCConstants(const Epetra_CrsMatrix & F, const Epetra_CrsMatrix & Bt,
                         const Epetra_CrsMatrix & B, const Epetra_CrsMatrix & C,
                         const Epetra_RowMatrix & Qu, double & gamma,double & alpha);

// This routine will construct the LSC matrices needed for both stabilized and
// stable discretizations. The difference will depend on if the A operator has
// a stabilization matrix or not. If it does not the out_aD parameter will not be
// set. For a reference on the LSC preconditioners see
//
//    Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
//    for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
//
// arguments (in):
//    in_A:  2x2 block matrix [ F Bt ; B C ] (see below for details)
//    in_Qu: Velocity mass matrix
// arguments (out):
//    out_BQBtmgC: Matrix (ML/IFPACK-ready) whose approximate inverse is needed by 
//                 preconditioner (BQBtmgC = B*inv(Q)*Bt - gamma*C ... See EHSST2007)
//    out_aiD:     vector representation of alpha*inv(D) (D = diag(B*inv(diag(F))*Bt-C)
//                 this will be set to null for a case stable discretization
//    
void buildLSCOperators(const Epetra_Operator & in_A,const Epetra_RowMatrix & in_Qu,
                       Teuchos::RCP<const Epetra_CrsMatrix> & out_BQBtmgC,Teuchos::RCP<const Epetra_Vector> & out_aiD);

void buildLSCOperators(const Epetra_Operator & in_A,Teuchos::RCP<const Epetra_CrsMatrix> & out_BQBtmgC,
                       Teuchos::RCP<const Epetra_Vector> & out_aiD);

// function to be used internally, handles both stable and stabilized discretizations
void buildLSCOperators(const Teuchos::RCP<const Epetra_CrsMatrix> & F,const Teuchos::RCP<const Epetra_CrsMatrix> & Bt,
                       const Teuchos::RCP<const Epetra_CrsMatrix> & B,const Teuchos::RCP<const Epetra_CrsMatrix> & C,
                       const Teuchos::RCP<const Epetra_RowMatrix> & in_Qu,
                       Teuchos::RCP<const Epetra_CrsMatrix> & out_BQBtmgC,Teuchos::RCP<const Epetra_Vector> & out_aiD);

// This routine will construct the LSC matrices needed for stabilized
// discretizations.  For a reference on the LSC preconditioners see
//
//    Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
//    for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
//
// arguments (in):
//    F:  Convection operator on discrete velocity
//    Bt: Gradient operator on discrete pressure
//    B:  Divergence operator on discrete velocity
//    C:  Stabilization operator
//    in_Qu: Velocity mass matrix
// arguments (out):
//    out_BQBtmgC: Matrix (ML/IFPACK-ready) whose approximate inverse is needed by 
//                 preconditioner (BQBtmgC = B*inv(Q)*Bt - gamma*C ... See EHSST2007)
//    out_aiD:     vector representation of alpha*inv(D) (D = diag(B*inv(diag(F))*Bt-C)
//    
inline void buildLSCOperators(const Epetra_CrsMatrix & F, const Epetra_CrsMatrix & Bt,
                       const Epetra_CrsMatrix & B, const Epetra_CrsMatrix & C,
                       const Epetra_RowMatrix & in_Qu,
                       Teuchos::RCP<const Epetra_CrsMatrix> & out_BQBtmgC,Teuchos::RCP<const Epetra_Vector> & out_aiD)
{
   buildLSCOperators(Teuchos::rcpFromRef(F),Teuchos::rcpFromRef(Bt),
                     Teuchos::rcpFromRef(B),Teuchos::rcpFromRef(C),Teuchos::rcpFromRef(in_Qu),out_BQBtmgC,out_aiD);
}



// convenience function for breaking a Thyra matrix (possibly blocked)
// composed of Epetra_CrsMatrix objects.  This is primarily intended to
// be used internally by the library.
//
// Preconditions: i) the matrix contains only Epetra_CrsMatrix or Thyra::zero
//                   objects.  Anything else will be treated as a Thyra::zero
//                   if failOnNonZero==false, other wise it wail throw a 
//                   casting exception. 
//               ii) blocks is an empty vector
//                 
// 
// arguments(in):
//    A:             Thyra matrix
//    failOnNonZero: throw an exception if a matrix is not CRS or Thyra::zero
// arguments(out):
//    blocks: vector of blocks of A cast to Epetra_CrsMatrix objects
// returns:
//    stl::pair<int,int> giving block dimension of matrix (this will be 
//    (1,1) for unblocked matrices)
//
std::pair<int,int> thyraMatrixToCrsVector(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A,bool failOnNonZero,
                                          std::vector<Teuchos::RCP<const Epetra_CrsMatrix> > & blocks);


/** \brief Build a vector of the dirchlet row indicies. 
  *
  * Build a vector of the dirchlet row indicies. That is, record the global
  * index of any row that is all zeros except for $1$ on the diagonal.
  *
  * \param[in] mat Matrix to be examined
  * \param[in,out] indicies Output list of indicies corresponding to dirchlet rows.
  */
void identityRowIndicies(const Epetra_CrsMatrix & mat,std::vector<int> & indicies);

} // end namespace Epetra
} // end namespace PB

#endif
