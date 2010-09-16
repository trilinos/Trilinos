#ifndef NOX_EPETRA_LINEARSYSTEMSGJACOBI_H
#define NOX_EPETRA_LINEARSYSTEMSGJACOBI_H

#include "NOX_Common.H"

#include "NOX_Epetra_LinearSystem.H"	// base class
#include "NOX_Utils.H"                  // class data element
#include "Epetra_Time.h"                // class data element
#include "Epetra_LinearProblem.h"
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
//! Improved version of the Epetra support class.
namespace Epetra {

/*! 

\brief Concrete implementation of NOX::Epetra::LinearSolver for Amesos. This class
has been written taking NOX::Epetra::LinearSystemAztecOO as an example.

The NOX::Epetra::LinearSystemAmesos object provides a way to 
interface an Epetra based application code to the Amesos linear solver. 
This class handles construction the Amesos solver. All options are determined
through parameter lists and the basic constructor. 

This class contains exclusively NOX::Epetra::LinearSolver virtual methods. Scaling of
the linear system is currently <em>not</em> implemented.

<B>Constructing a Linear System</B>

At the time being there is only one constructor that can be used: among other things, 
parameter, the user must provide a NOX printing parameter sublist, 
a NOX linear solver sublist, and a Jacobian. This constructor is analogous to one of the
NOX:Epetra:LinearSystemAmesos ones.

The user can select an Amesos solver type via the parameter <B>Amesos Solver</B> in 
the NOX linear solver sublist. The default solver is "Amesos_Klu".
 */
class LinearSystemSGJacobi : public virtual NOX::Epetra::LinearSystem {

  public:

    //! Constructor with a user supplied Jacobian Operator  
    LinearSystemSGJacobi(
      Teuchos::ParameterList& printingParams, 
      Teuchos::ParameterList& linearSolverParams, 
      const Teuchos::RCP<NOX::Epetra::LinearSystem>& detsolve,
      const Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> >& Cijk, 
      const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq, 
      const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac, 
      const Teuchos::RCP<Epetra_Operator>& J,
      const NOX::Epetra::Vector& cloneVector,
      const NOX::Epetra::Vector& detcloneVector,
      const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject = 
      Teuchos::null);

    //! Destructor.
    virtual ~LinearSystemSGJacobi();

  /*! 
    \brief Applies Jacobian to the given input vector and puts the answer in the result.

    Computes 
    \f[ v = J u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
    and \f$v\f$ is the result vector.  Returns true if successful.
  */
    virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
			       NOX::Epetra::Vector& result) const;

  /*!  
    \brief Applies Jacobian-Transpose to the given input vector and puts the answer in the result.

    Computes 
    \f[ v = J^T u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result vector.  Returns true if successful.
    
  */
    virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
					NOX::Epetra::Vector& result) const;

  /*!
    \brief Applies the inverse of the Jacobian matrix to the given
    input vector and puts the answer in result.
    
    Computes 
    \f[ v = J^{-1} u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
    and \f$v\f$ is the result vector.
    
    The parameter list contains the linear solver options.
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
    virtual void setJacobianOperatorForSolve(const 
		   Teuchos::RCP<const Epetra_Operator>& solveJacOp);

    //! Does nothing
    virtual void setPrecOperatorForSolve(const 
		   Teuchos::RCP<const Epetra_Operator>& solvePrecOp);
    
  protected:

    //! Pointer to deterministic solver
    Teuchos::RCP<NOX::Epetra::LinearSystem> detsolve;
   
   //! Pointer to Cijk
    Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk;
 
    //! Reference to the user supplied Jacobian interface functions
    Teuchos::RCP<NOX::Epetra::Interface::Jacobian> jacInterfacePtr;
    
    //! Pointer to the Jacobian operator.
    mutable Teuchos::RCP<Epetra_Operator> jacPtr;

    //! Pointer to the Stokhos matrixfree epetra operator.
    mutable Teuchos::RCP<Stokhos::SGOperator> stokhos_op;

    //! Pointer to the PCE expansion of Jacobian.
    mutable Teuchos::RCP<Stokhos::VectorOrthogPoly<Epetra_Operator> > sg_J_poly;
   
    //! Pointer to the left-hand side of the Stokhos problem
    Teuchos::RCP<Epetra_Vector> leftHandSide;

    //! Pointer to the right-hand side of the Stokhos problem
    Teuchos::RCP<Epetra_Vector> rightHandSide;

    //! Pointer to the deterministic vector 
    Teuchos::RCP<Epetra_Vector> detvec;

    //! Scaling object supplied by the user
    Teuchos::RCP<NOX::Epetra::Scaling> scaling;

    //! Epetra_Time object
    Epetra_Time timer;

    //! Printing Utilities object
    NOX::Utils utils;
};

} // namespace Epetra
} // namespace NOX

#endif
