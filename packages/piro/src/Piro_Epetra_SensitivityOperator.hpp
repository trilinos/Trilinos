// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef PIRO_EPETRA_SENSITIVITY_OPERATOR_H
#define PIRO_EPETRA_SENSITIVITY_OPERATOR_H

#include "Epetra_Operator.h"
#include "EpetraExt_ModelEvaluator.h"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "NOX_Epetra_Group.H"
#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H"

namespace Piro {
  namespace Epetra {

    /*!
     * \brief Epetra_Operator representing dg/dp
     */
    class SensitivityOperator : public Epetra_Operator {
    public:

      /** \name Constructors/initializers */
      //@{

      /** \brief Takes the number of elements in the discretization . */
      SensitivityOperator(
	const Teuchos::RCP<const Epetra_Map>& g_map,
	const Teuchos::RCP<const Epetra_Map>& p_map);

      void setup(
	const EpetraExt::ModelEvaluator::Derivative& dfdp,
	const EpetraExt::ModelEvaluator::Derivative& dgdx,
	const EpetraExt::ModelEvaluator::Derivative& dgdp,
	Teuchos::RCP<Teuchos::ParameterList>& piroParams,
	const Teuchos::RCP<NOX::Epetra::Group>& grp,
	const Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy>& tls_strategy);

      //@}

      ~SensitivityOperator();


      /** \name Overridden from Epetra_Operator . */
      //@{

      //! Set to true if the transpose of the operator is requested
      virtual int SetUseTranspose(bool UseTranspose);
      
      /*! 
       * \brief Returns the result of a Epetra_Operator applied to a 
       * Epetra_MultiVector Input in Result as described above.
       */
      virtual int Apply(const Epetra_MultiVector& Input, 
			Epetra_MultiVector& Result) const;
      
      /*! 
       * \brief Returns the result of the inverse of the operator applied to a 
       * Epetra_MultiVector Input in Result as described above.
       */
      virtual int ApplyInverse(const Epetra_MultiVector& X, 
			       Epetra_MultiVector& Y) const;
      
      //! Returns an approximate infinity norm of the operator matrix.
      virtual double NormInf() const;
      
      //! Returns a character string describing the operator
      virtual const char* Label () const;
      
      //! Returns the current UseTranspose setting.
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

      //@}

    private:
    
      //! Private to prohibit copying
      SensitivityOperator(const SensitivityOperator&);
      
      //! Private to prohibit copying
      SensitivityOperator& operator=(const SensitivityOperator&);

      int applyForward(const Epetra_MultiVector& Input, 
		       Epetra_MultiVector& Result) const;
      
      int applyAdjoint(const Epetra_MultiVector& Input, 
		       Epetra_MultiVector& Result) const;

    private:

      //! Label for operator
      std::string label;

      //! Response map (range map for forward sensitivities)
      Teuchos::RCP<const Epetra_Map> g_map;

      //! Parameter map (domain map for forward sensitivities)
      Teuchos::RCP<const Epetra_Map> p_map;

      //! Whether to use the transpose
      bool useTranspose;

      //! Stores df/dp as multi-vector or operator
      EpetraExt::ModelEvaluator::Derivative dfdp;

      //! Stores dg/dx as multi-vector or operator
      EpetraExt::ModelEvaluator::Derivative dgdx;

      //! Stores dg/dp as multi-vector or operator
      EpetraExt::ModelEvaluator::Derivative dgdp;

      //! Group from solver for applyJacobianInverse
      Teuchos::RCP<NOX::Epetra::Group> grp;

      //! Transpose solver
      Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy;

      //! Piro parameters for solvers
      Teuchos::RCP<Teuchos::ParameterList> piroParams;

    };

  }
}
#endif
