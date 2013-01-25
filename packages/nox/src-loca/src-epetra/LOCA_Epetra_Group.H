// $Id: LOCA_Epetra_Group.H,v 1.40 2007/06/21 16:22:55 rhoope Exp $ 
// $Source: /space/CVS/Trilinos/packages/nox/src-loca/src-epetra/LOCA_Epetra_Group.H,v $ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LOCA_EPETRA_GROUP_H
#define LOCA_EPETRA_GROUP_H

#include "LOCA_Abstract_Group.H"               // base class
#include "LOCA_Abstract_TransposeSolveGroup.H" // base class
#include "NOX_Epetra_Group.H"                  // base class
#include "LOCA_Parameter_Vector.H"             // class data element
#include "NOX_Common.H"                        // class data element (std::string)
#include "LOCA_Epetra_Interface_TimeDependent.H"
#include "LOCA_Epetra_Interface_FreeEnergy.H"
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
class Epetra_BlockMap;

namespace LOCA {

  //! %LOCA %Epetra support
  namespace Epetra {

    //! Extension of the NOX::Epetra::Group to %LOCA.  
    /*!
     * This class extends the NOX::Epetra::Group to %LOCA enabling continuation
     * and bifurcation capabilities using %Epetra.  It is derived from
     * the NOX::Epetra::Group (basic %Epetra support), the
     * LOCA::Abstract::Group (brings in all %LOCA abstract base classes),
     * and the LOCA::Abstract::TransposeSolveGroup (for 
     * applyJacobianTransposeInverse() methods).  It stores a parameter vector
     * for setting/retrieving parameter values and overloads the computeF() 
     * and computeJacobian() methods of the NOX::Epetra::Group
     * parent class to set the entire contents
     * of the parameter vector in the problem interface before calling the
     * NOX::Epetra::Group computeF() and computeJacobian().  

     * Since it is derived from the LOCA::Abstract::Group (which is in turn 
     * derived from all FiniteDifference groups), it uses the finite-difference
     * implementations for all parameter derivatives and second derivatives.
     * However this behavior can be modified by calling the setDerivUtils()
     * method of the LOCA::MultiContinuation::FiniteDifferenceGroup parent
     * class.
     *
     * This class provides complete support for all continuation and 
     * bifurcation methods including shift-invert and Cayley methods
     * for computing eigenvalues and Hopf bifurcations.  However this support
     * is only enabled by calling the appropriate constructor described below.
     */
    class Group : public virtual NOX::Epetra::Group, 
		  public virtual LOCA::Abstract::Group,
		  public virtual LOCA::Abstract::TransposeSolveGroup {

    public:
      //! Constructor with NO linear system (VERY LIMITED).
      /*! 
       * WARNING: If this constructor is used, then methods that require
       * a Jacobian or preconditioning will not be available.  You will be
       * limited to simple algorithms like nonlinear-CG with no
       * preconditioning.
       */
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RCP<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess,
	    const LOCA::ParameterVector& p);

      //! Standard Constructor enabling most %LOCA support
      /*! 
       * This is the most commonly used constructor and provides support for
       * all %LOCA algorithms except shift-invert and Cayley transformations
       * and Hopf bifurcations.
       */
      Group(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	    Teuchos::ParameterList& printingParams, 
	    const Teuchos::RCP<LOCA::Epetra::Interface::Required>& i, 
	    NOX::Epetra::Vector& initialGuess, 
	    const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
	    const LOCA::ParameterVector& p);

      //! Constructor with time-dependent interface and shifted linear system
      /*!
       * Use this constructor to enable shift-invert and Cayley transformations
       * or Hopf bifurcations.  It requires another interface to compute the
       * shifted matrix \f$\alpha J + \beta M\f$ where \f$J\f$ is the 
       * Jacobian matrix and \f$M\f$ is the mass matrix, and a linear
       * system object to solve this system.  Setting linSys = shiftedLinSys
       * is a valid option for passing the shifted solver, but this will cause
       * the shifted matrix to overwrite the Jacobian possibly resulting in
       * more matrix fills.
       * See declareSeparateMatrixMemory() method below to assert separate memory.
       */
      Group(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	Teuchos::ParameterList& printingParams,
        const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent>& i,
	NOX::Epetra::Vector& initialGuess,
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& shiftedLinSys,
	const LOCA::ParameterVector& p);

      /*!
       * \brief Constructor with time-dependent matrix-free interface and 
       * shifted linear system
       */
      /*!
       * This constructor may also be used for shift-invert and Cayley
       * transformations, but should be only be used for a matrix-free
       * method for solving the shifted system.
       */
      Group(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	Teuchos::ParameterList& printingParams,
        const Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree>& i,
	NOX::Epetra::Vector& initialGuess,
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& linSys,
	const Teuchos::RCP<NOX::Epetra::LinearSystem>& shiftedLinSys,
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
      virtual Teuchos::RCP<NOX::Abstract::Group>
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
      /*!
       * In addition to all regular linear solver parameters, this method
       * references the following additional parameters:
       * <ul>
       * <li> "Transpose Solver Method" -- [string] 
       *      (default: "Transpose Preconditioner") Method for preconditioning 
       *      the transpose linear system 
       *      (LOCA::Epetra::TransposeLinearSystem::Factory).  Available 
       *      choices are:
       *   <ul>
       *   <li> "Transpose Preconditioner" -- Use the transpose of the 
       *        preconditioner for the original system.
       *   <li> "Left Preconditioning" -- Use the transpose of the 
       *        preconditioner, and apply using left preconditioning.
       *   <li> "Explicit Transpose" -- Form the transpose of the matrix and 
       *        compute the preconditioner.  This method is available only if 
       *        Trilinos is configured with %EpetraExt support 
       *        (--enable-epetraext).
       *   </ul>
       * </ul>
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeInverse(Teuchos::ParameterList& params, 
				    const NOX::Abstract::Vector& input, 
				    NOX::Abstract::Vector& result) const;

      //! Solve Jacobian-tranpose system with multiple right-hand sides
      /*!
       * In addition to all regular linear solver parameters, this method
       * references the following additional parameters:
       * <ul>
       * <li> "Transpose Solver Method" -- [string] 
       *      (default: "Transpose Preconditioner") Method for preconditioning 
       *      the transpose linear system 
       *      (LOCA::Epetra::TransposeLinearSystem::Factory).  Available 
       *      choices are:
       *   <ul>
       *   <li> "Transpose Preconditioner" -- Use the transpose of the 
       *        preconditioner for the original system.
       *   <li> "Left Preconditioning" -- Use the transpose of the 
       *        preconditioner, and apply using left preconditioning.
       *   <li> "Explicit Transpose" -- Form the transpose of the matrix and 
       *        compute the preconditioner.  This method is available only if 
       *        Trilinos is configured with %EpetraExt support 
       *        (--enable-epetraext).
       *   </ul>
       * </ul>
       */
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
      virtual void setParam(std::string paramID, double val);

      //! Return a const reference to the ParameterVector owned by the group. 
      const LOCA::ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(std::string paramID) const;

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.  The implementation here is to call the corresponding
       * method in the interface.
       */
      virtual void 
      preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus);

      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful. The implementation here is to call the corresponding
       * method in the interface.
       */
      virtual void 
      postProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus);

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
       * \brief Replace Jacobian \f$J\f$ by \f$aJ+bI\f$ where \f$I\f$ is 
       * the identity matrix.
       */ 
      virtual NOX::Abstract::Group::ReturnType 
      augmentJacobianForHomotopy(double a, double b);

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

      //! Compute the second shifted matrix (uses different memory then Shifted matrix)
      virtual NOX::Abstract::Group::ReturnType
      computeSecondShiftedMatrix(double alpha, double beta);

      //! Multiply the second shifted matrix by a vector.
      virtual NOX::Abstract::Group::ReturnType
      applySecondShiftedMatrix(const NOX::Abstract::Vector& input,
			 NOX::Abstract::Vector& result) const;

      //! Multiply the second shifted matrix by a multi-vector.
      virtual NOX::Abstract::Group::ReturnType
      applySecondShiftedMatrixMultiVector(
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

      /*! 
       * @name Implementation of LOCA::Hopf::MinimallyAugmented::AbstractGroup virtual methods.
       */
      //@{

      /*!
       * Computes conjugate-tranpose matrix vector product 
       * \f$ (J+i\omega B)^H (x + iy) \f$.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyComplexTranspose(const NOX::Abstract::Vector& input_real,
			    const NOX::Abstract::Vector& input_imag,
			    NOX::Abstract::Vector& result_real,
			    NOX::Abstract::Vector& result_imag) const;

      /*!
       * Computes conjugate-tranpose matrix vector product 
       * \f$ (J+i\omega B)^H (x + iy) \f$.
       */
      virtual NOX::Abstract::Group::ReturnType 
      applyComplexTransposeMultiVector(
			   const NOX::Abstract::MultiVector& input_real,
			   const NOX::Abstract::MultiVector& input_imag,
			   NOX::Abstract::MultiVector& result_real,
			   NOX::Abstract::MultiVector& result_imag) const;

      //! Solve \f$(J+i\omega B)^H (x + iy) = a+ib\f$
      virtual NOX::Abstract::Group::ReturnType
      applyComplexTransposeInverseMultiVector(
			    Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& input_real,
			    const NOX::Abstract::MultiVector& input_imag,
			    NOX::Abstract::MultiVector& result_real,
			    NOX::Abstract::MultiVector& result_imag) const;

      //@}

      //! Method to inject an interface for calucatiuong the free energy
      virtual void setFreeEnergyInterface(
              const Teuchos::RCP<LOCA::Epetra::Interface::FreeEnergy>& iFE);

      //! Method for calling code to guarantee to LOCA that separate matrix
      //  memory was used for linSys and shiftedLinSys (see 3rd constructor)
      void declareSeparateMatrixMemory(bool separateMem=true);

      /*! 
       * @name Implementation of LOCA::PhseTransition::AbstractGroup virtual methods.
       */
      //@{

      virtual double computeFreeEnergy();

      //@}

      //! Return the userInterface.
      virtual Teuchos::RCP<NOX::Epetra::Interface::Required> 
      getUserInterface();

      //! Call the user interface print() routine, any vector
      virtual void printSolution(const NOX::Epetra::Vector& x,
                                 const double conParam) const;

      //! Sets the scale vector.  
      void setScaleVector(const NOX::Abstract::Vector& s);

      //! Sets the Jacobian operator
      void setJacobianOperatorForSolve(
		 const Teuchos::RCP<const Epetra_Operator>& op) const;

      //! Return the Linear System.
      virtual Teuchos::RCP<const NOX::Epetra::LinearSystem> 
      getComplexLinearSystem() const;

      //! Return the Linear System.
      virtual Teuchos::RCP<NOX::Epetra::LinearSystem> 
      getComplexLinearSystem();

      virtual void
      getComplexMaps(Teuchos::RCP<const Epetra_BlockMap>& baseMap,
		     Teuchos::RCP<const Epetra_BlockMap>& globalMap) const;

    protected:

      //! resets the isValid flags to false
      virtual void resetIsValid();

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Printing parameters
      Teuchos::ParameterList& printParams;

      //! Parameter vector
      LOCA::ParameterVector params;

      //! Reference to the user supplied interface functions
      Teuchos::RCP<LOCA::Epetra::Interface::Required> userInterface;

      //! Interface for shifted-matrix
      Teuchos::RCP<LOCA::Epetra::Interface::TimeDependent> userInterfaceTime;

      //! Interface for matrix-free shifted-matrix
      Teuchos::RCP<LOCA::Epetra::Interface::TimeDependentMatrixFree> userInterfaceTimeMF;

      //! Interface for free enerfgy calculation for phase transitions
      Teuchos::RCP<LOCA::Epetra::Interface::FreeEnergy> userInterfaceFreeEnergy;

      //! Shared shifted linear system
      Teuchos::RCP< NOX::SharedObject<NOX::Epetra::LinearSystem, 
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
      Teuchos::RCP<Epetra_Vector> tmpVectorPtr2;

      //! Stores a pointer to the scale vector
      Teuchos::RCP<NOX::Abstract::Vector> scaleVecPtr;

      //! Stores transpose linear solver strategy
      mutable Teuchos::RCP<LOCA::Epetra::TransposeLinearSystem::AbstractStrategy> tls_strategy;

      //! Shared complex system
      mutable Teuchos::RCP< NOX::SharedObject<NOX::Epetra::LinearSystem, LOCA::Epetra::Group> > complexSharedLinearSystem;

      //! Complex matrix
      Teuchos::RCP<EpetraExt::BlockCrsMatrix> complexMatrix;

      //! Complex vector
      Teuchos::RCP<EpetraExt::BlockVector> complexVec;

      //! Is complex matrix valid
      bool isValidComplex;

      //! Is complex matrix preconditioner valid
      mutable bool isValidComplexPrec;

      //! Is complex matrix valid
      bool separateMatrixMemoryDeclared;
    };

  } // namespace Epetra
} // namespace LOCA


#endif
