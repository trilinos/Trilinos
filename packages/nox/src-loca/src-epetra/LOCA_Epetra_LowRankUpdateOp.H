// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER
                                                                                
#ifndef LOCA_EPETRA_LOWRANKUPDATEOP_H
#define LOCA_EPETRA_LOWRANKUPDATEOP_H

#include "Teuchos_RefCountPtr.hpp"

#include "Epetra_Operator.h"	// base class
#include "Epetra_LocalMap.h"    // class data element
#include "Epetra_MultiVector.h" // class data element
#include "NOX_Common.H"         // for string data member

// forward declarations
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Epetra {
    
    /*! 
     * \brief An Epetra operator for implementing the operator
     * \f$P = J + U V^T\f$.
     */
    /*!
     * This class implements the Epetra_Operator interface for 
     * \f$P = J + U V^T\f$ where \f$J\f$ is an Epetra_Operator and 
     * \f$U\f$ and \f$V\f$ are Epetra_MultiVectors.  
     */
    class LowRankUpdateOp : public virtual Epetra_Operator {
      
    public:

      //! Constructor 
      /*!
       * \param global_data [in] The global data object
       * \param jacOperator [in] Jacobian operator J
       * \param U_multiVec [in] Multivector representing U
       * \param V_multiVec [in] Multivector representing V
       */
      LowRankUpdateOp(
        const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	const Teuchos::RefCountPtr<Epetra_Operator>& jacOperator, 
	const Teuchos::RefCountPtr<const Epetra_MultiVector>& U_multiVec, 
	const Teuchos::RefCountPtr<const Epetra_MultiVector>& V_multiVec);

      //! Destructor
      virtual ~LowRankUpdateOp();

      //! Set to true if the transpose of the operator is requested
      virtual int SetUseTranspose(bool UseTranspose);

      /*! 
       * \brief Returns the result of a Epetra_Operator applied to a 
       * Epetra_MultiVector Input in Result as described above.
       */
      virtual int Apply(const Epetra_MultiVector& Input, 
			Epetra_MultiVector& Result) const;

      /*!
       * \brief This method does nothing.
       */
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const;
  
      //! Returns an approximate infinity norm of the operator matrix.
      /*!
       * This is defined only if %NormInf() of the underlying operator \f$J\f$
       * is defined and is given by 
       * \f$\|J\|_\infty+\|U\|_\infty\|V\|_\infty\f$.
       */
      virtual double NormInf() const;
    
      //! Returns a character string describing the operator
      virtual const char* Label () const;
  
      //! Returns the current UseTranspose setting.  Always returns false.
      virtual bool UseTranspose() const;

      /*! 
       * \brief Returns true if the \e this object can provide an 
       * approximate Inf-norm, false otherwise.
       */
      virtual bool HasNormInf() const;

      /*! 
       * \brief Returns a reference to the Epetra_Comm communicator 
       * associated with this operator.
       */
      virtual const Epetra_Comm & Comm() const;

      /*!
       * \brief Returns the Epetra_Map object associated with the 
       * domain of this matrix operator.
       */
      virtual const Epetra_Map& OperatorDomainMap () const;

      /*! 
       * \brief Returns the Epetra_Map object associated with the 
       * range of this matrix operator.
       */
      virtual const Epetra_Map& OperatorRangeMap () const;

    private:

      //! Private to prohibit copying
      LowRankUpdateOp(const LowRankUpdateOp&);

      //! Private to prohibit copying
      LowRankUpdateOp& operator=(const LowRankUpdateOp&);

    protected:

      //! Global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Label for operator
      string label;

      //! Local map for generating temporary matrices
      Epetra_LocalMap localMap;
    
      //! Stores operator representing J
      Teuchos::RefCountPtr<Epetra_Operator> J;
  
      //! Stores multivector representing U
      Teuchos::RefCountPtr<const Epetra_MultiVector> U;

      //! Stores multivector representing V
      Teuchos::RefCountPtr<const Epetra_MultiVector> V;

      //! Flag indicating whether to use the transpose
      bool useTranspose;

      //! Temporary matrix
      mutable Teuchos::RefCountPtr<Epetra_MultiVector> tmpMat;

    };
  } // namespace Epetra
} // namespace LOCA

#endif 
