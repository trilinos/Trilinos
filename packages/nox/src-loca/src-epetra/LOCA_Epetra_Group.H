// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef LOCA_EPETRA_GROUP_H
#define LOCA_EPETRA_GROUP_H

#include "LOCA_Abstract_Group.H"               // base class
#include "LOCA_Abstract_TransposeSolveGroup.H" // base class
#include "NOX_Epetra_Group.H"                  // base class
#include "LOCA_Parameter_Vector.H"             // class data element
#include "NOX_Common.H"                        // class data element (string)
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "LOCA_Epetra_Interface_TimeDependentMatrixFree.H"

// Forward declares
namespace Teuchos {
  class ParameterList;
}
namespace LOCA {
  namespace Epetra {
    namespace Interface {
      class Required;
    }
    namespace TransposeLinearSystem {
      class AbstractStrategy;
    }
  }
}
namespace EpetraExt {
  class BlockCrsMatrix;
  class BlockVector;
}

namespace LOCA {

  //! %LOCA %Epetra support
  namespace Epetra {

    //! Extension of the NOX::Epetra:New:Group to %LOCA.  
    /*!
      This class is derived both from the NOX::Epetra::Group and
      LOCA::Abstract::Group classes and therefore inherits the implementation
      of the NOX::Abstract::Group interface provided by NOX::Epetra::Group.
      
      This class provides implementations of %LOCA AbstractGroup virtual 
      methods specific to the %Epetra group.  It stores a parameter vector
      for setting/retrieving parameter values 
      (LOCA::Continuation::AbstractGroup), provides a facility for computing
      eigenvalues (LOCA::Continuation::AbstractGroup) using the Anasazi 
      package, and augements the Jacobian matrix for homotopy
      (LOCA::Homotopy::AbstractGroup), It currently has no support for
      storing a mass matrix and therefore cannot be used in Hopf or 
      generalized eigenvalue computations.  Since it is
      derived from the LOCA::Abstract::Group (which is in-turn derived
      from all FiniteDifference groups), this group implicitly uses the
      finite-difference implementations of parameter and second derivatives
      provided by the FiniteDifference groups.  

      The computeF() and computeJacobian() methods of the NOX::Epetra::Group
      parent class are overloaded here.  They both set the entire contents
      of the parameter vector in the problem interface before calling the
      NOX::Epetra::Group computeF() and computeJacobian().
    */

    class Group : public virtual NOX::Epetra::Group, 
		  public virtual LOCA::Abstract::Group,
		  public virtual LOCA::Abstract::TransposeSolveGroup {

    public:
      //! Constructor with NO linear system (VERY LIMITED).
      /*! WARNING: If this constructor is used, then methods that require
       * a Jacobian or preconditioning will not be available.  You will be
       * limited to simple algorithms like nonlinear-CG with no
       * preconditioning.
       */
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess,
	    const LOCA::ParameterVector& p);

      //! Standard Constructor.
      Group(const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	    const LOCA::ParameterVector& p);

      //! Constructor with time-dependent interface and shifted linear system
      Group(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	Teuchos::ParameterList& printingParams,
        const Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependent>& i,
	NOX::Epetra::Vector& initialGuess,
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& shiftedLinSys,
	const LOCA::ParameterVector& p);

      //! Constructor with time-dependent matrix-free interface and shifted linear system
      Group(
	const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	Teuchos::ParameterList& printingParams,
        const Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependentMatrixFree>& i,
	NOX::Epetra::Vector& initialGuess,
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linSys,
	const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& shiftedLinSys,
	const LOCA::ParameterVector& p);


      /*! 
       * \brief Copy constructor. If type is DeepCopy, takes ownership of
       * valid shared Jacobian and shared preconditioning matrix. 
       */
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      //! Assignment operator
      virtual Group& operator=(const Group& source);

      /*! 
       * @name Overloaded NOX::Epetra::Group  methods.
       */
      //@{

      //! Assignment operator.
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator.
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Epetra::Group& source);

      //! Cloning function
      virtual Teuchos::RefCountPtr<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Overloaded computeF()
      /*!
       * Calls LOCA::Epetra::Interface::setParams before evalulating F.
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeF();

      //! Overloaded computeJacobian()
      /*!
       * Calls LOCA::Epetra::Interface::setParams before evalulating J.
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeJacobian();

      //@}

      /*! 
       * @name Implementation of LOCA::Abstract::TransposeSolveGroup methods.
       */
      //@{

      //! Solve Jacobian-tranpose system 
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverse(Teuchos::ParameterList& params, 
				    const NOX::Abstract::Vector& input, 
				    NOX::Abstract::Vector& result) const;

      //! Solve Jacobian-tranpose system with multiple right-hand sides
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup virtual methods.
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set the parameters
      virtual void setParams(const ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);
  
      //! Set parameter indexed by paramID
      virtual void setParam(string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(string paramID) const;

      //! Projects solution to a few scalars for multiparameter continuation
      /*! 
       * This method is called every time a solution is saved by the 
       * multiparameter continuation code MF for later visualization
       * and should project the solution vector down to a few scalars.
       * The array \c px will be preallocated to the proper length
       * given by projectToDrawDimension().
       *
       * The implementation here is to call the corresponding method
       * in the interface.
       */
      virtual void projectToDraw(const NOX::Abstract::Vector& x,
				 double *px) const;

      //! Returns the dimension of the project to draw array
      /*!
       * The implementation here is to call the corresponding method 
       * in the interface.
       */
      virtual int projectToDrawDimension() const;

      //! Compute a scaled dot product
      /*! 
       * The implementation here uses the scaling vector \f$s\f$ if one
       * is supplied:
       * \f[
       *     \sum_{i=1}^n a_i*b_i*s_i*s_i.
       * \f]
       * If the scaling vector is not provided, the standard dot product
       * is used.
       */
      virtual double
      computeScaledDotProduct(const NOX::Abstract::Vector& a,
			      const NOX::Abstract::Vector& b) const;

      //! Call the user interface print() routine, solution vector
      virtual void printSolution(const double conParam) const;

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Abstract::Vector& x,
                                 const double conParam) const;

      //! Scales a vector using scaling vector
      /*! 
       * The implementation here uses the scaling vector \f$s\f$ if one
       * is supplied:
       * \f[
       *     x_i = a_i*s_i.
       * \f]
       * If the scaling vector is not provided, the vector is rescaled by
       * the square root of its length.
       */
      virtual void
      scaleVector(NOX::Abstract::Vector& x) const;

      //@}

      /*! 
       * @name Implementation of LOCA::Homotopy::AbstractGroup virtual methods.
       */
      //@{

     /*! 
       * \brief Replace Jacobian \f$J\f$ by \f$pJ+(1-p)I\f$ where \f$I\f$ is 
       * the identity matrix and \f$p\f$ is a scalar.
       */ 
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double p);

      //@}

      /*! 
       * @name Implementation of LOCA::TimeDependent::AbstractGroup virtual methods.
       */
      //@{

      //! Compute the shifted matrix
      virtual NOX::Abstract::Group::ReturnType
      computeShiftedMatrix(double alpha, double beta);

      //! Multiply the shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const;

      //! Multiply the shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixMultiVector(
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;
      
      /*! 
       * \brief Apply the inverse of the shifted matrix by a multi-vector, as 
       * needed by the shift-and-invert and generalized Cayley transformations.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyShiftedMatrixInverseMultiVector(
			        Teuchos::ParameterList& params, 
				const NOX::Abstract::MultiVector& input,
				NOX::Abstract::MultiVector& result) const;
      //@}

      /*! 
       * @name Implementation of LOCA::Hopf::MooreSpence::AbstractGroup virtual methods.
       */
      //@{

      //! Is  \f$J+i\omega B\f$ valid
      virtual bool isComplex() const;

      //! Compute \f$J+i\omega B\f$
      /*!
       * The argument \b frequency stores \f$\omega\f$.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeComplex(double frequency);
      
      //! Compute \f$(J+i\omega B)(y+iz)\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplex(const NOX::Abstract::Vector& input_real,
		   const NOX::Abstract::Vector& input_imag,
		   NOX::Abstract::Vector& result_real,
		   NOX::Abstract::Vector& result_imag) const;

      //! Compute \f$(J+i\omega B)(y+iz)\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexMultiVector(
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const;

      //! Solve \f$(J+i\omega B)(y+iz) = a+ib\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexInverseMultiVector(
			    Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const;

      //@}

      //! Return the userInterface.
      virtual Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> 
      getUserInterface();

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Epetra::Vector& x,
                                 const double conParam) const;

      //! Sets the scale vector.  
      void setScaleVector(const NOX::Abstract::Vector& s);

      //! Sets the Jacobian operator
      void setJacobianOperatorForSolve(
		 const Teuchos::RefCountPtr<const Epetra_Operator>& op) const;

    protected:

      //! resets the isValid flags to false
      virtual void resetIsValid();

    protected:

      //! Global data
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Printing parameters
      Teuchos::ParameterList& printParams;

      //! Parameter vector
      LOCA::ParameterVector params;

      //! Reference to the user supplied interface functions
      Teuchos::RefCountPtr<LOCA::Epetra::Interface::Required> userInterface;

      //! Interface for shifted-matrix
      Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependent> userInterfaceTime;

      //! Interface for matrix-free shifted-matrix
      Teuchos::RefCountPtr<LOCA::Epetra::Interface::TimeDependentMatrixFree> userInterfaceTimeMF;

      //! Shared shifted linear system
      Teuchos::RefCountPtr< NOX::SharedObject<NOX::Epetra::LinearSystem, 
					      LOCA::Epetra::Group> > shiftedSharedLinearSystem;

      //! Is preconditioner for shifted matrix valid
      mutable bool isValidShiftedPrec;

      //! \f$\alpha\f$ for matrix-free shifted matrix
      double alpha_;

      //! \f$\beta\f$ for matrix-free shifted matrix
      double beta_;

      //! Extra vector needed for intermediate calculations of LOCA routines.
      /*! NOTE: there already is a tmpVectorPtr in the NOX::Epetra::Group.  
       * This is a second temporary vector if that one extra isn't enough.
       */
      Teuchos::RefCountPtr<Epetra_Vector> tmpVectorPtr2;

      //! Stores a pointer to the scale vector
      Teuchos::RefCountPtr<NOX::Abstract::Vector> scaleVecPtr;

      //! Stores transpose linear solver strategy
      mutable Teuchos::RefCountPtr<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy;

      //! Shared complex system
      mutable Teuchos::RefCountPtr< NOX::SharedObject<NOX::Epetra::LinearSystem, LOCA::Epetra::Group> > complexSharedLinearSystem;

      //! Complex matrix
      Teuchos::RefCountPtr<EpetraExt::BlockCrsMatrix> complexMatrix;

      //! Complex vector
      Teuchos::RefCountPtr<EpetraExt::BlockVector> complexVec;

      //! Is complex matrix valid
      bool isValidComplex;

      //! Is complex matrix preconditioner valid
      mutable bool isValidComplexPrec;

    };

  } // namespace Epetra
} // namespace LOCA


#endif
