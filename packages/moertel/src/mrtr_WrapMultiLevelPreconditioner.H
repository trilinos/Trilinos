/*
# ************************************************************************
#
#               ML: A Multilevel Preconditioner Package
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file mrtr_WrapMulitLevelPreconditioner.H
 *
 * \class MOERTEL::ConstrainedPreconditioner
 *
 * \brief A wrapper for ML's MultiLevelPreconditioner that enforces Mortar constraints
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#ifndef MOERTEL_WRAPPRECONDITIONER_H
#define MOERTEL_WRAPPRECONDITIONER_H

#include <ctime>
#include <iostream>
#include <vector>

// Trilinos includes
#include <Epetra_Operator.h>
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

#include "mrtr_utils.H"

using namespace std;
using namespace Teuchos;

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
\class ConstrainedPreconditioner

\brief A wrapper for ML's MultiLevelPreconditioner that applies mortar constraints
       to current residual and preconditioned search direction on the fine grid


\author Michael W. Gee (mwgee@sandia.gov)

*/
class  ConstrainedPreconditioner : public virtual Epetra_Operator
{
public:

  // @{ \name Constructors and destructors
  

  /*!
  \brief 
  
  */
  explicit ConstrainedPreconditioner(RefCountPtr<ML_Epetra::MultiLevelPreconditioner> mlprec,
                                     RefCountPtr<Epetra_CrsMatrix> I,
                                     RefCountPtr<Epetra_CrsMatrix> WT,
                                     RefCountPtr<Epetra_CrsMatrix> B);

  
  /*!
  \brief Destroys an instance of this class
  
  Destructor
  
  */
  virtual ~ConstrainedPreconditioner() {}

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
  const Epetra_Comm & Comm() const {return(mlprec_->Comm());}

  /*!
  \brief Get fine level OperatorDomainMap
  
  Derived from Epetra_Operator, get fine level OperatorDomainMap
  
  */
  const Epetra_Map & OperatorDomainMap() const {return mlprec_->OperatorDomainMap();}

  /*!
  \brief Get fine level OperatorRangeMap
  
  Derived from Epetra_Operator, get fine level OperatorRangeMap
  
  */
  const Epetra_Map & OperatorRangeMap() const {return mlprec_->OperatorRangeMap();}
  
  //@}
  // @{ \name Apply methods


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
  { cout << "ConstrainedPreconditioner::Apply does not make sense for an ML-Multilevel Operator\n"
         << "Use ApplyInverse\n";throw -1;return(-1);}
    
  /*!
  \brief not implemented
  */
  int SetUseTranspose(bool UseTranspose) { cout << "ConstrainedPreconditioner::SetUseTranspose not impl.\n";return -1;}

  /*!
  \brief not implemented
  */
  double NormInf() const {cout << "ConstrainedPreconditioner::NormInf not impl.\n";
                          return(-1.0);}


  /*!
  \brief not implemented
  */
  bool UseTranspose() const {cout << "ConstrainedPreconditioner::UseTranspose not impl.\n"; return false;}

  /*!
  \brief not implemented
  */
  bool HasNormInf() const {cout << "ConstrainedPreconditioner::HasNormInf not impl.\n"; return false;}

  //@}

private:
  // don't want copy-ctor and = operator
  ConstrainedPreconditioner(ConstrainedPreconditioner& old);
  ConstrainedPreconditioner operator = (const ConstrainedPreconditioner& old);
  
  char*                                             label_;  // the label of this class
  RefCountPtr<ML_Epetra::MultiLevelPreconditioner>  mlprec_; // the ml preconditioner
  RefCountPtr<Epetra_CrsMatrix>                     WT_;
  RefCountPtr<Epetra_CrsMatrix>                     B_;
  RefCountPtr<Epetra_CrsMatrix>                     I_;  



}; // class  ConstrainedPreconditioner : public virtual Epetra_Operator
}  // namespace ConstrainedPreconditioner


#endif // MOERTEL_WRAPPRECONDITIONER_H
