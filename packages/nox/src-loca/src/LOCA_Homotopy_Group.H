// $Id$
// $Source$

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

#ifndef LOCA_HOMOTOPY_GROUP_H
#define LOCA_HOMOTOPY_GROUP_H

#include "Teuchos_RCP.hpp"

#include "LOCA_MultiContinuation_AbstractGroup.H"    // base class
#include "LOCA_Extended_MultiAbstractGroup.H"        // base class

#include "LOCA_Parameter_Vector.H"  // class data member

namespace LOCA {
  class GlobalData;
  namespace Homotopy {
    class AbstractGroup;
  }
}

namespace LOCA { 

  namespace Homotopy {
  
    /*! 
     * \brief %LOCA's Homotopy Algorithm.  
     */
    /*!
     * The %HomotopyGroup is a concrete implementation of the 
     * LOCA::Continuation::AbstractGroup that modifies the set of nonlinear 
     * equations to be solved to allow for Homotopy to be applied to the 
     * system.  This object should be used in conjunction with the 
     * LOCA::Stepper object to drive the continuation.  This algorithm solves 
     * a system of nonlinear equations supplied by the user (\f$ F(x) \f$) 
     * through continuation.  An artificial parameter \f$ \lambda \f$ is used 
     * to control the continuation.  The idea is to solve a simple equation 
     * starting at \f$ \lambda \f$ = 0 and, using the solution from the 
     * previous step, solve systems of equations that gets progressively 
     * closer to the true system of interest ( at \f$ \lambda \f$ = 1.0 we 
     * recover the original equations \f$ F(x) \f$).  By constraining the 
     * definition of \f$ g(x, \lambda) \f$ and using artificial parameter 
     * contiuation, the continuation branch should be free of multiplicity and 
     * bifurcation phenomena.
     * 
     * The modified system of equations, \f$ g(x, \lambda) \f$, supplied by 
     * the HomotopyGroup is defined as:
     *
     * \f[ g(x, \lambda) = \lambda F(x) + (1.0 - \lambda)(x - a) \f]
     *
     * where \f$x\f$ is the solution vector, \f$ \lambda \f$ is an artificial 
     * parameter, \f$ F(x) \f$ is the set of nonlinear equations the user 
     * supplies, \f$ g(x) \f$ is the corresponding set of 
     * homotopy equations that LOCA will solve, and \f$ a \f$ is a random 
     * vector.
     *
     * This group requires the loca Stepper for continuation
     * from \f$ \lambda \f$ = 0.0 (a simple set of equations to solve) to 
     * \f$ \lambda \f$ = 1.0 (the set of equations requested by the user, 
     * \f$ F(x) \f$).  The Homotopy::Group will generate the Stepper parameter 
     * sublist in the parameter list that is passed in to the constructor.  
     * The user is free to modify this list (it sets default values) before 
     * passing it into the stepper object but should NOT change the starting 
     * and stopping values for the continuation parameter.
     *
     * References:
     * 
     * - ALGORITHM 652 HOMPACK: A Suite of Codes for Globally Convergent 
     * Homotopy Algorithms, Watson, L.T., Billups, S.C, and Morgan, A.P., 
     * ACM Transactions on Mathematical Software, Vol. 13, No. 3, September 
     * 1987, pp281-310.
     */

    class Group : 
      public LOCA::MultiContinuation::AbstractGroup,
      public LOCA::Extended::MultiAbstractGroup {

    public:
  
      /*! 
       * \brief Constructor to set the base group and generate the "%Stepper" 
       * sublist for homotopy continuation.
       */
      /*! 
       * The locaSublist variable is the "LOCA" sublist (of type 
       * Teuchos::ParameterList) that will be used in loca continuation runs.
       *
       * The variables scalarRandomVector and scalarInitialGuess are used to 
       * give some control over the generation of the random vector.  In 
       * certain instances we have seen the random vector force the solution 
       * to a set of variables that are unphysical and could break the 
       * function evaluations (cause them to return nan).  For example, in 
       * heat transfer problems, the temperature could be the dependent 
       * variable.  If the solution vector has an unphysical temperature (
       * the random vector could force the temperature to negative or near 
       * zero values for the solution at \f$ \lambda = 0\f$) then property 
       * evaluations could break.  The random vector can be modified to keep 
       * the values near the initial guess based on values supplied to the 
       * constructor of the HomotopyGroup: 
       * 
       * \f[
       * a = abs(r) * \mbox{scalarRandom} + x_o * \mbox{scalarInitialGuess}
       * \f]
       * 
       * where \f$ r \f$ is the random vector generated by a call to 
       * NOX::Abstract::Vector::random(), \f$ \mbox{scalarRandom} \f$ is a 
       * scalar value, \f$ x_o \f$ is the initial guess to the solution 
       * vector, and \f$ \mbox{scalarInitialGuess} \f$ is a scalar value.  
       * The defualt values force the random vector to be calculated as:
       *
       * \f[
       * a = abs(r)
       * \f]
       *
       * IMPORTANT: For homotopy to work correctly you should not change the 
       * starting and stopping parameter values (0.0 and 1.0 respectively) set 
       * in the "%Stepper" sublist.
       */
      Group(
	 Teuchos::ParameterList& locaSublist,
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Homotopy::AbstractGroup>& g,
	 double scaleRandom = 1.0,
	 double scaleInitialGuess = 0.0);

      //! Constructor with a user supplied random vector.
      Group(
	 Teuchos::ParameterList& locaSublist,
	 const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Homotopy::AbstractGroup>& g,
	 const NOX::Abstract::Vector& randomVector);

      //! Copy constructor.
      Group(const Group& source, NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~Group();

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);
     
      //! Cloning function
      virtual Teuchos::RCP<NOX::Abstract::Group> 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //! Set the solution vector, x, to y.
      virtual void setX(const NOX::Abstract::Vector& y);

      //! Compute this.x = grp.x + step * d.
      virtual void computeX(const NOX::Abstract::Group& g, 
			    const NOX::Abstract::Vector& d,
			    double step);

      //! Compute the homotopy residual $g$ 
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! Compute the Jacobian derivative of the homotopy residual $g$
      virtual NOX::Abstract::Group::ReturnType computeJacobian();

      //! Compute gradient of homotopy residual $g$
      virtual NOX::Abstract::Group::ReturnType computeGradient();

      //! Compute %Newton direction using applyJacobianInverse
      virtual NOX::Abstract::Group::ReturnType 
      computeNewton(Teuchos::ParameterList& params);
  
      //! Computes the homotopy Jacobian vector product
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobian(const NOX::Abstract::Vector& input, 
		    NOX::Abstract::Vector& result) const;

      //! Computes the homotopy Jacobian-transpose vector product
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;

      //! Applies the inverse of the homotopy Jacobian matrix 
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverse(Teuchos::ParameterList& params, 
			   const NOX::Abstract::Vector& input, 
			   NOX::Abstract::Vector& result) const;

      //! Applies Jacobian for homotopy system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
			       NOX::Abstract::MultiVector& result) const;

      //! Applies Jacobian-transpose for homotopy system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianTransposeMultiVector(
				     const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const;

      //! Applies Jacobian inverse for homotopy system
      virtual NOX::Abstract::Group::ReturnType 
      applyJacobianInverseMultiVector(
				    Teuchos::ParameterList& params, 
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;
  
      //! Return \c true if the homotopy residual \f$g\f$ is valid.
      virtual bool isF() const;

      //! Return \c true if the homotopy Jacobian is valid.
      virtual bool isJacobian() const;

      //!  Return \c true if the homotopy gradient is valid.
      virtual bool isGradient() const;

      //! Return \c true if the homotopy Newton direction is valid.
      virtual bool isNewton() const;

      //! Return homotopy solution vector \f$x\f$.  
      virtual const NOX::Abstract::Vector& getX() const;

      //! Return homotopy residual \f$g\f$
      virtual const NOX::Abstract::Vector& getF() const;

      //! Return 2-norm of \f$g\f$. 
      virtual double getNormF() const;

      //! Return homotopy gradient
      virtual const NOX::Abstract::Vector& getGradient() const;

      //! Return homotopy Newton direction.
      virtual const NOX::Abstract::Vector& getNewton() const;

      //! Return homotopy solution vector \f$x\f$.  
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const;

      //! Return homotopy residual \f$g\f$
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const;

      //! Return homotopy gradient
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const;

      //! Return homotopy Newton direction.
      virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const;

      //@}

      /*! 
       * @name Implementation of LOCA::Extended::MultiAbstractGroup 
       * virtual methods 
       */
      //@{

      //! Return underlying group
      virtual 
      Teuchos::RCP<const LOCA::MultiContinuation::AbstractGroup> 
      getUnderlyingGroup() const;
      
	//! Return underlying group
      virtual 
      Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>
      getUnderlyingGroup();
      
      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractGroup
       * virtual methods 
       */
      //@{

      //! Assignment
      virtual void copy(const NOX::Abstract::Group& source);

      //! Set parameters indexed by (integer) paramIDs
      virtual void setParamsMulti(
			  const std::vector<int>& paramIDs, 
			  const NOX::Abstract::MultiVector::DenseMatrix& vals);

      //! Set the parameter vector in the group to p.
      virtual void setParams(const ParameterVector& p);

      //! Set parameter indexed by paramID
      virtual void setParam(int paramID, double val);

      //! Set parameter indexed by paramID
      virtual void setParam(std::string paramID, double val);

      //! Return a const reference to the paramter vector owned by the group.
      virtual const ParameterVector& getParams() const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(int paramID) const;

      //! Return copy of parameter indexed by paramID
      virtual double getParam(std::string paramID) const;

      /*! 
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti(const std::vector<int>& paramIDs, 
		       NOX::Abstract::MultiVector& dfdp, 
		       bool isValidF);

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.
       */
      virtual void 
      preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus);
      
      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.
       */
      virtual void 
      postProcessContinuationStep(
			 LOCA::Abstract::Iterator::StepStatus stepStatus);

      //! Projects solution to a few scalars for multiparameter continuation
      virtual void projectToDraw(const NOX::Abstract::Vector& x,
				 double *px) const;

      //! Returns the dimension of the project to draw array
      virtual int projectToDrawDimension() const;

      /*! 
       * \brief Function to print out solution and continuation 
       * parameter after successful continuation step.
       */
      virtual void printSolution(const double conParam) const;

      /*! 
       * \brief Function to print out solution and continuation 
       * parameter after successful continuation step.
       */
      virtual void printSolution(const NOX::Abstract::Vector& x_,
				 const double conParam) const;

      //@}

    protected:
      
      //! Reset the isValid flags to false.  
      /*! This is called when the solution vector or parameter vector 
	is changed.
      */
      void resetIsValidFlags();

      //! Creates and sets the "Stepper" parameter sublist
      void setStepperParameters(Teuchos::ParameterList& params);

    private:

      //! Prohibit generation and use of operator=()
      Group& operator=(const Group&);

    protected:

      //! Pointer LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Stores the underlying loca group.
      Teuchos::RCP<LOCA::Homotopy::AbstractGroup> grpPtr;   
  
      //! Stores the homotopy residual vector, \f$ g \f$.
      Teuchos::RCP<NOX::Abstract::Vector> gVecPtr;

      //! Stores the random Vector, \f$ a \f$.
      Teuchos::RCP<NOX::Abstract::Vector> randomVecPtr;

      /*! 
       * \brief Stores the homotopy Newton vector, 
       * \f$ \frac{\partial r}{\partial x} \f$.
       */
      Teuchos::RCP<NOX::Abstract::Vector> newtonVecPtr;

      /*! 
       * Stores the homotopy gradient vector if needed, 
       * \f$ \frac{\partial r}{\partial x} \f$.
       */
      Teuchos::RCP<NOX::Abstract::Vector> gradVecPtr;

      //! Is residual vector valid
      bool isValidF;

      //! Is Jacobian matrix valid
      bool isValidJacobian;

      //! Is Newton vector valid
      bool isValidNewton;

      //! Is gradient vector valid
      bool isValidGradient;
      
      //! Copy of the ParameterVector for the underlying grpPtr.
      /*! We copy this and then add the homotopy parameter to the list. */
      LOCA::ParameterVector paramVec;
      
      //! Value of the homotopy continuation parameter.
      /*! Ranges from 0.0 (easy solution) to 1.0 (solution to the system of 
       * interest).
       */
      double conParam;

      //! Continuatioin parameter ID number from the ParameterVector.
      int conParamID;
      
      /*! 
       * \brief Contains the std::string used to identify the homotopy parameter in 
       * the ParameterVector object.
       */
      const std::string conParamLabel;

      /*!
       * \brief Tracks whether the LOCA::Homotopy::Group method 
       * augmentJacobianForHomotopy is implemented.  If not, the augmentation 
       * is applied during the applyJacobian assuming a matrix-free 
       * implementation.
       */
      bool augmentJacForHomotopyNotImplemented;

    }; // class Group

  } // namespace Homotopy

} // namespace LOCA

#endif
