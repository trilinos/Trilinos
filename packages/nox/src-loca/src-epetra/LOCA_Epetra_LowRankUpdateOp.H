// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_EPETRA_LOWRANKUPDATEOP_H
#define LOCA_EPETRA_LOWRANKUPDATEOP_H

#include "Teuchos_RCP.hpp"
#include "Teuchos_LAPACK.hpp"

#include "Epetra_Operator.h"    // base class
#include "Epetra_LocalMap.h"    // class data element
#include "Epetra_MultiVector.h" // class data element
#include "NOX_Common.H"         // for std::string data member

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
        const Teuchos::RCP<LOCA::GlobalData>& global_data,
    const Teuchos::RCP<Epetra_Operator>& jacOperator,
    const Teuchos::RCP<const Epetra_MultiVector>& U_multiVec,
    const Teuchos::RCP<const Epetra_MultiVector>& V_multiVec,
    bool setup_for_solve);

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

      //! Returns a character std::string describing the operator
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
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Label for operator
      std::string label;

      //! Local map for generating temporary matrices
      Epetra_LocalMap localMap;

      //! Stores operator representing J
      Teuchos::RCP<Epetra_Operator> J;

      //! Stores multivector representing U
      Teuchos::RCP<const Epetra_MultiVector> U;

      //! Stores multivector representing V
      Teuchos::RCP<const Epetra_MultiVector> V;

      //! Flag indicating whether to use the transpose
      bool useTranspose;

      //! Temporary matrix
      mutable Teuchos::RCP<Epetra_MultiVector> tmpMat;

      //! Stores J^{-1}*U
      Teuchos::RCP<Epetra_MultiVector> JinvU;

      //! Stores LU factorization of I + V^T*J^{-1}*U
      Teuchos::RCP<Epetra_MultiVector> lu;

      //! Stores pivots for LU factorization
      std::vector<int> ipiv;

      //! Lapack wrappers
      Teuchos::LAPACK<int,double> lapack;

    };
  } // namespace Epetra
} // namespace LOCA

#endif
