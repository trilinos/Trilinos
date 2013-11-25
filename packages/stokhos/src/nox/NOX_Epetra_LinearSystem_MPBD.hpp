/*
// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER
*/

#ifndef NOX_EPETRA_LINEARSYSTEMMPBD_H
#define NOX_EPETRA_LINEARSYSTEMMPBD_H

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_NOX

#include "NOX_Common.H"

#include "NOX_Epetra_LinearSystem.H"	// base class
#include "NOX_Utils.H"                  // class data element

#include "Stokhos_ProductContainer.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_BlockDiagonalOperator.hpp"
#include "EpetraExt_BlockVector.h"

namespace NOX {
  namespace Epetra {
    namespace Interface {
      class Required;
      class Jacobian;
      class Preconditioner;
    }
  }
}

namespace NOX {

  namespace Epetra {

    /*! 
     * \brief Concrete implementation of NOX::Epetra::LinearSolver for 
     * multi-point block diagonal solves
     */
    class LinearSystemMPBD : public virtual NOX::Epetra::LinearSystem {
    public:

      //! Constructor
      LinearSystemMPBD(
	Teuchos::ParameterList& printingParams, 
	Teuchos::ParameterList& linearSolverParams, 
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& block_solver,
	const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
	const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
	const Teuchos::RCP<Epetra_Operator>& J,
	const Teuchos::RCP<const Epetra_Map>& base_map,
	const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
	Teuchos::null);

      //! Destructor.
      virtual ~LinearSystemMPBD();

      /*! 
       * \brief Applies Jacobian to the given input vector and puts the 
       * answer in the result.
       */
      virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
				 NOX::Epetra::Vector& result) const;

      /*!  
       * \brief Applies Jacobian-Transpose to the given input vector and puts
       * the answer in the result.
       */
      virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
					  NOX::Epetra::Vector& result) const;
      
      /*!
       * \brief Applies the inverse of the Jacobian matrix to the given
       * input vector and puts the answer in result.
       */
      virtual bool applyJacobianInverse(Teuchos::ParameterList &params, 
					const NOX::Epetra::Vector &input, 
					NOX::Epetra::Vector &result);
      
      //! Returns false
      virtual bool applyRightPreconditioning(bool useTranspose,
					     Teuchos::ParameterList& params, 
					     const NOX::Epetra::Vector& input, 
					     NOX::Epetra::Vector& result) const;
      
      //! Returns supplied scaling object
      virtual Teuchos::RCP<NOX::Epetra::Scaling> getScaling();
      
      //! Reset supplied scaling object
      virtual void resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s);
      
      //! Evaluates the Jacobian based on the solution vector x.
      virtual bool computeJacobian(const NOX::Epetra::Vector& x);
      
      //! Returns false
      virtual bool createPreconditioner(const NOX::Epetra::Vector& x, 
					Teuchos::ParameterList& p,
					bool recomputeGraph) const;
      
      //! Returns false
      virtual bool destroyPreconditioner() const;
      
      //! Returns false
      virtual bool recomputePreconditioner(const NOX::Epetra::Vector& x, 
					   Teuchos::ParameterList& linearSolverParams) const;

      //! Returns PRPT_REUSE;
      virtual PreconditionerReusePolicyType 
      getPreconditionerPolicy(bool advanceReuseCounter=true);
      
      //! Returns false
      virtual bool isPreconditionerConstructed() const;
      
      //! Returns false
      virtual bool hasPreconditioner() const;
      
      //! Returns jacobian operator
      virtual Teuchos::RCP<const Epetra_Operator> 
      getJacobianOperator() const;
      
      //! Returns jacobian operator
      virtual Teuchos::RCP<Epetra_Operator> getJacobianOperator();
      
      //! Returns Teuchos::null
      virtual Teuchos::RCP<const Epetra_Operator> 
      getGeneratedPrecOperator() const;
      
      //! Returns Teuchos::null
      virtual Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator();
      
      //! Resets the jacobian operator
      virtual void setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp);

      //! Does nothing
      virtual void setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp);
      
    protected:

      enum PREC_STRATEGY {
	STANDARD,
	MEAN,
	ON_THE_FLY
      };
      
      //! Pointer to block solver
      Teuchos::RCP<NOX::Epetra::LinearSystem> block_solver;
      
      //! Reference to the user supplied Jacobian interface functions
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterfacePtr;

      //! Number of multi-point blocks
      int num_mp_blocks;
      
      //! Pointer to the Stokhos block diagonal epetra operator.
      mutable Teuchos::RCP<Stokhos::BlockDiagonalOperator> mp_op;
      
      //! Pointer to the block operators
      mutable Teuchos::RCP<const Stokhos::ProductContainer<Epetra_Operator> > block_ops;
      
      //! Stores base map
      Teuchos::RCP<const Epetra_Map> base_map;
      
      //! Scaling object supplied by the user
      Teuchos::RCP<NOX::Epetra::Scaling> scaling;
      
      //! Printing Utilities object
      NOX::Utils utils;

      PREC_STRATEGY precStrategy;

      //! Preconditioner operator for each point
      mutable Teuchos::Array< Teuchos::RCP<const Epetra_Operator> > precs;

      //! x-vector for on-the-fly preconditioner strategy
      mutable Teuchos::RCP<EpetraExt::BlockVector> prec_x;

    };
    
  } // namespace Epetra
} // namespace NOX

#endif

#endif
