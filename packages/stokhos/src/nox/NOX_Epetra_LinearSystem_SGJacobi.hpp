// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef NOX_EPETRA_LINEARSYSTEMSGJACOBI_H
#define NOX_EPETRA_LINEARSYSTEMSGJACOBI_H

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_NOX

#include "NOX_Common.H"

#include "NOX_Epetra_LinearSystem.H"	// base class
#include "NOX_Utils.H"                  // class data element

#include "Stokhos_ParallelData.hpp"
#include "Stokhos_OrthogPolyBasis.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Stokhos_VectorOrthogPoly.hpp"
#include "Stokhos_VectorOrthogPolyTraitsEpetra.hpp"
#include "Stokhos_SGOperator.hpp"
#include "Stokhos_Sparse3Tensor.hpp"

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
     * stochastic Galerkin systems using Gauss-Seidel iterations.
     */
    class LinearSystemSGJacobi : public virtual NOX::Epetra::LinearSystem {

    public:

      //! Constructor
      LinearSystemSGJacobi(
	Teuchos::ParameterList& printingParams, 
	Teuchos::ParameterList& linearSolverParams, 
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& detsolve,
	const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
	const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
	const Teuchos::RCP<const Stokhos::OrthogPolyBasis<int,double> >& sg_basis,
	const Teuchos::RCP<const Stokhos::ParallelData>& sg_parallel_data,
	const Teuchos::RCP<Epetra_Operator>& J,
	const Teuchos::RCP<const Epetra_Map>& base_map,
	const Teuchos::RCP<const Epetra_Map>& sg_map,
	const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
	Teuchos::null);

      //! Destructor.
      virtual ~LinearSystemSGJacobi();

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

      //! Pointer to deterministic solver
      Teuchos::RCP<NOX::Epetra::LinearSystem> det_solver;

      //! Stores Epetra Cijk tensor
      Teuchos::RCP<const Stokhos::EpetraSparse3Tensor> epetraCijk;
   
      //! Reference to the user supplied Jacobian interface functions
      Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterfacePtr;
      
      //! Pointer to the Stokhos matrixfree epetra operator.
      mutable Teuchos::RCP<Stokhos::SGOperator> sg_op;

      //! Pointer to the PCE expansion of Jacobian.
      mutable Teuchos::RCP<Stokhos::EpetraOperatorOrthogPoly > sg_poly;

      //! Stores base map
      Teuchos::RCP<const Epetra_Map> base_map;

      //! Stores SG map
      Teuchos::RCP<const Epetra_Map> sg_map;

      //! Scaling object supplied by the user
      Teuchos::RCP<NOX::Epetra::Scaling> scaling;

      //! Printing Utilities object
      NOX::Utils utils;
      
      //! Temporary vector used in Gauss-Seidel iteration
      mutable Teuchos::RCP<EpetraExt::BlockVector> sg_df_block;
      
      //! Temporary vector used in Gauss-Seidel iteration
      mutable Teuchos::RCP<Epetra_Vector> kx;

      //! SG operator to implement SG mat-vec
      Teuchos::RCP<Stokhos::SGOperator> mat_free_op;
    };
    
  } // namespace Epetra
} // namespace NOX

#endif

#endif
