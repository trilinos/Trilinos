/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_ml_preconditioner.H
 *
 * \class MOERTEL::Mortar_ML_Preconditioner
 *
 * \brief 
 *
 * \date Last update do Doxygen: 20-March-06
 *
 */
#ifndef MOERTEL_ML_PRECONDITIONER_H
#define MOERTEL_ML_PRECONDITIONER_H

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_Operator.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_RowMatrix.h>
#ifdef EPETRA_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "ml_common.h"
#include "ml_include.h"
#include "ml_epetra_utils.h"
#include "ml_epetra.h"
#include "ml_epetra_operator.h"
#include "ml_MultiLevelPreconditioner.h"
#include "ml_agg_genP.h"
#include "MLAPI_Error.h"
#include "MLAPI_CompObject.h"
#include "MLAPI_TimeObject.h"
#include "MLAPI_Operator.h"
#include "MLAPI_Operator_Utils.h"
#include "MLAPI_MultiVector.h"
#include "MLAPI_InverseOperator.h"
#include "MLAPI_Expressions.h"
#include "MLAPI_BaseOperator.h"
#include "MLAPI_Workspace.h"
#include "MLAPI_Aggregation.h"
#include "MLAPI_Eig.h"

#include "mrtr_utils.H"

/*!
\brief MOERTEL: namespace of the Moertel package

The Moertel package depends on \ref Epetra, \ref EpetraExt, \ref Teuchos,
\ref Amesos, \ref ML and \ref AztecOO:<br>
Use at least the following lines in the configure of Trilinos:<br>
\code
--enable-moertel 
--enable-epetra 
--enable-epetraext
--enable-teuchos 
--enable-ml
--enable-aztecoo --enable-aztecoo-teuchos 
--enable-amesos
\endcode

*/
namespace MOERTEL
{
/*!
\class Mortar_ML_Preconditioner

\brief 


\author Glen Hansen (Glen.Hansen@inl.gov)

*/
class  Mortar_ML_Preconditioner : public virtual Epetra_Operator
{
public:

  // @{ \name Constructors and destructors
  

  /*!
  \brief Construct the mortar modified ml preocnditioner
  
   Construct the mortar modified ml preconditioner using ML's MLAPI
   
   
   \param Atilde (In) : The mortar modified spd system of equations
   \param A (In) : The original (uncoupled) system matrix
   \param WT (In) : Mortar scaling object 
   \param B (In) : Mortar constraints equations 
   \param params (In) : ML parameters 
   \param constructnow (In) : If true, the preconditioner is constructed immediately 
  */
  explicit Mortar_ML_Preconditioner(Teuchos::RCP<Epetra_CrsMatrix> Atilde,
                                    Teuchos::RCP<Epetra_CrsMatrix> A,
                                    Teuchos::RCP<Epetra_CrsMatrix> WT,
                                    Teuchos::RCP<Epetra_CrsMatrix> B,
                                    Teuchos::RCP<Epetra_Map>       Annmap,
                                    Teuchos::ParameterList& mlparams,
                                    bool constructnow = true);

  
  /*!
  \brief Destroys an instance of this class
  
  Destructor
  
  */
  virtual ~Mortar_ML_Preconditioner() {}

  //@}
  // @{ \name Query methods

  /*!
  \brief Returns the label of this class.
  */
  const char* Label() const {return(&(label_[0]));}

  /*!
  \brief get Comm of this class
  
  Derived from Epetra_Operator, returns ref to the Epetra_Comm of this class 
  
  */
  const Epetra_Comm & Comm() const {return(Atilde_->Comm());}

  /*!
  \brief Get fine level OperatorDomainMap
  
  Derived from Epetra_Operator, get fine level OperatorDomainMap
  
  */
  const Epetra_Map & OperatorDomainMap() const {return Atilde_->OperatorDomainMap();}

  /*!
  \brief Get fine level OperatorRangeMap
  
  Derived from Epetra_Operator, get fine level OperatorRangeMap
  
  */
  const Epetra_Map & OperatorRangeMap() const {return Atilde_->OperatorRangeMap();}
  
  //@}
  // @{ \name Apply methods


  /*!
  \brief Compute the preconditioner
  
   Compute the mortar modified ml preconditioner
  */
  bool Compute();  

  /*!
  \brief ApplyInverse the preconditioner
  
   ApplyInverse the preconditioner. Method is derived from Epetra_Operator. 
   
   
   \param X   (In) : Epetra_MultiVector matching the fine level map of this
                     preconditioner
   \param Y (Out)  : Epetra_MultiVector containing the result on output
  */
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const ;

  // @{ \name Not implemented methods derived from Epetra_Operator

  /*!
  \brief not implemented
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const 
  { cout << "Mortar_ML_Preconditioner::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";throw -1;return(-1);}
    
  /*!
  \brief not implemented
  */
  int SetUseTranspose(bool UseTranspose) { cout << "Mortar_ML_Preconditioner::SetUseTranspose not impl.\n";return -1;}

  /*!
  \brief not implemented
  */
  double NormInf() const {cout << "Mortar_ML_Preconditioner::NormInf not impl.\n";
                          return(-1.0);}


  /*!
  \brief not implemented
  */
  bool UseTranspose() const {cout << "Mortar_ML_Preconditioner::UseTranspose not impl.\n"; return false;}

  /*!
  \brief not implemented
  */
  bool HasNormInf() const {cout << "Mortar_ML_Preconditioner::HasNormInf not impl.\n"; return false;}

  //@}

private:
  // don't want copy-ctor and = operator
  Mortar_ML_Preconditioner(Mortar_ML_Preconditioner& old);
  Mortar_ML_Preconditioner operator = (const Mortar_ML_Preconditioner& old);
  
  // the v-cycle
#if 1 // old working version
  int MultiLevelSA(const MLAPI::MultiVector& b_f, MLAPI::MultiVector& x_f, int level) const;
#endif
#if 0 // experimental version I
  int MultiLevelSA(const MultiVector& b1_f,
                   const MultiVector& b2_f,
                   MultiVector& x1_f, 
                   MultiVector& x2_f, 
                   int level) const;
#endif
#if 0 // experimental version II
  int MultiLevelSA(const MultiVector& b_f,MultiVector& x_f, int level) const;
#endif

  // some helpers
  inline const MLAPI::Operator& P(const int i) const { return (mlapiPmod_[i]); }
  inline const MLAPI::Operator& R(const int i) const { return (mlapiRmod_[i]); }
  inline const MLAPI::Operator& A(const int i) const { return (mlapiAtilde_[i]); }
  inline const MLAPI::Operator& ImBWT(const int i) const { return (mlapiImBWT_[i]); }
  inline const MLAPI::Operator& ImWBT(const int i) const { return (mlapiImWBT_[i]); }
  inline const MLAPI::InverseOperator& S(const int i) const { return (mlapiS_[i]); }

  const char*                    label_;       // the label of this class
  bool                           iscomputed_;  // flag indicating whether preconditioner is ready for use
  Teuchos::ParameterList&                 mlparams_;    // ML parameter list
  int                            maxlevels_;   // number of levels

  Teuchos::RCP<Epetra_CrsMatrix>  Atilde_;      // mortar modified spd system
  Teuchos::RCP<Epetra_CrsMatrix>  A_;           // original uncoupled system
  Teuchos::RCP<Epetra_CrsMatrix>  WT_;          // mortar scaling
  Teuchos::RCP<Epetra_CrsMatrix>  B_;           // mortar constraints

  Teuchos::RCP<Epetra_Map>        Annmap_;      // row map of Ann

  std::vector<MLAPI::Operator>               mlapiImBWT_;  // operator I - B * W^TT on each grid                     
  std::vector<MLAPI::Operator>               mlapiImWBT_;  // operator I - W * B^T on each grid                     
  std::vector<MLAPI::Operator>               mlapiRmod_;   // modified restriction                     
  std::vector<MLAPI::Operator>               mlapiPmod_;   // modified prolongation                    
  std::vector<MLAPI::Operator>               mlapiAtilde_; // operator Atilde on each grid                       
  std::vector<MLAPI::InverseOperator>        mlapiS_;      // smoother on each grid

#if 0 // additional data for experimental version I 
  Teuchos::RCP<Epetra_Map>        Arrmap_;      // row map of Arr
  Teuchos::RCP<Epetra_CrsMatrix>  Atilde11_;
  Teuchos::RCP<Epetra_CrsMatrix>  Ahat11_;
  Teuchos::RCP<Epetra_CrsMatrix>  M_;
  Teuchos::RCP<Epetra_CrsMatrix>  MT_;
  Operator                       mlapiAtilde11_;
  Operator                       mlapiAhat11_;
  Operator                       mlapiM_;
  Operator                       mlapiMT_;
  InverseOperator                G_;
#endif  

#if 0 // additional data for experimental version II
  Teuchos::RCP<Epetra_Map>        Arrmap_;      // row map of Arr
  Teuchos::RCP<Epetra_CrsMatrix>  Asplit_;      // uncoupled version of A
  Teuchos::RCP<Epetra_CrsMatrix>  Atildesplit_; // uncoupled modified version of Atilde
  Teuchos::RCP<Epetra_CrsMatrix>  Ahat_;        
  std::vector<Operator>               mlapiAtildesplit_;  
  std::vector<Operator>               mlapiAhat_;  
  inline const Operator& Ahat(const int i) const { return (mlapiAhat_[i]); }
#endif  

}; // class  ConstrainedPreconditioner : public virtual Epetra_Operator
}  // namespace ConstrainedPreconditioner


#endif // MOERTEL_ML_PRECONDITIONER_H
