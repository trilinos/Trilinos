 // $Id$ 
 // $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
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

#ifndef NOX_ABSTRACT_GROUP_H
#define NOX_ABSTRACT_GROUP_H

#include "NOX_Abstract_Vector.H" // for NOX::CopyType
#include "NOX_Common.H" // for std::string 
#include "Teuchos_RCP.hpp" 
#include "Teuchos_Assert.hpp"

namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Abstract {
    class MultiVector;
  }
}

namespace NOX { 
namespace Abstract { 

 /*!  
   \brief %NOX pure abstract interface to a "group"; i.e., a
   solution vector and the corresponding F-vector, Jacobian matrix,
   gradient vector, and Newton vector.

   This class is a member of the namespace NOX::Abstract. 

   The user should implement their own concrete implementation of this
   class or use one of the implementations provided by us. Typically
   the implementation is also tied to a particular
   NOX::Abstract::Vector implementation.

   \note The group may be implemented so that multiple groups can
   share underlying memory space. This is particularly important when
   it comes to the Jacobian, which is often to big to be replicated
   for every group. Thus, we have included instructions on how
   <em>shared data</em> should be treated for the operator=() and
   clone() functions.
 */

 class Group {

 public:

   /*! 
     \brief The computation of, say, the Newton direction in
     computeNewton() may fail in many different ways, so we have
     included a variety of return codes to describe the failures. Of
     course, we also have a code for success.

     \note These return types may be expanded in future releases.
    */
   enum ReturnType {

     //! Computation completed successfully
     Ok,

     //! This function is not implemented
     NotDefined,

     //! Data dependencies not satisfied
     BadDependency,

     //! Unable to satisfy convergence criteria
     NotConverged,

     //! Any other type of failure
     Failed

   };

   //! Constructor.
   /*!  
     \note Constructors for any derived object should always define a default
     x-value so that getX() is always defined.
   */
   Group() {};

   //! Destructor.
   virtual ~Group() {};

   /*! 
     \brief Copies the source group into this group.

     \note Any <em>shared data</em> owned by the source should have
     its ownership transfered to this group. This may result in
     a secret modification to the source object.
   */
   virtual NOX::Abstract::Group& operator=(const NOX::Abstract::Group& source) = 0;

   //@{ \name "Compute" functions.

   //! Set the solution vector x to y.
   /*!
     \note
     This should invalidate the function value, Jacobian,
     gradient, and Newton direction.

     \note
     Throw an error if the copy fails.

     \return
     Reference to this object 
   */
   virtual void setX(const NOX::Abstract::Vector& y) = 0;

   //! Compute x = grp.x + step * d.
   /*! 
     Let \f$x\f$ denote this group's solution vector.
     Let \f$\hat x\f$ denote the result of grp.getX().
     Then set
     \f[ 
     x = \hat x + \mbox{step} \; d. 
     \f]

     \note
     This should invalidate the function value, Jacobian,
     gradient, and Newton direction.

     \note
     Throw an error if the copy fails.

     \return
     Reference to this object 

   */
   virtual void computeX(const NOX::Abstract::Group& grp, 
			 const NOX::Abstract::Vector& d, double step) = 0;

   //! Compute and store F(x). 
   /*! 
     \note 
     It's generally useful to also compute and store the 2-norm of F(x)
     at this point for later access by the getNormF() function.

     \return
     <ul>
     <li> NOX::Abstract::Group::Failed - If the computation fails in any way
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeF() = 0;

   //! Compute and store Jacobian.
   /*! 
     Recall that 
     \f[
     F(x) = \left[ \begin{array}{c} F_1(x) \\ F_2(x) \\ \vdots \\ F_n(x) \\ \end{array} \right].
     \f]
     
     The \b Jacobian
     is denoted by \f$J\f$ and defined by
     \f[
     J_{ij} = \frac{\partial F_i}{\partial x_j} (x).
     \f]
     
     \note 
     If this is a <em>shared object</em>, this group should taken
     ownership of the Jacobian before it computes it.

     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::Failed - If the computation fails in any other way
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeJacobian();

   //! Compute and store gradient.
   /*! 
     We can pose the nonlinear equation problem \f$F(x) = 0\f$ as an optimization problem as follows:
     \f[
     \min f(x) \equiv \frac{1}{2} \|F(x)\|_2^2.
     \f]

     In that case, the \b gradient (of \f$f\f$) is defined as 
     \f[
     g \equiv J^T F.
     \f]

     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If either \f$F\f$ or \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::Failed - If the computation fails in any other way
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeGradient();

   //! Compute the Newton direction, using parameters for the linear solve.
   /*! 
     The <b>Newton direction</b> is the solution, s, of
     \f[
     J s = -F.
     \f]
     
     The parameters are from the "Linear %Solver" sublist of the
     "Direction" sublist that is passed to solver during construction.
     
     The "Tolerance" parameter may be added/modified in the sublist of
     "Linear Solver" parameters that is passed into this function. The
     solution should be such that
     \f[
     \frac{\| J s - (-F) \|_2}{\max \{ 1, \|F\|_2\} } < \mbox{Tolerance}
     \f]
     
     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If either \f$F\f$ or \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance" 
	  specified in \c params
     <li> NOX::Abstract::Group::Failed - If the computation fails in any other way
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeNewton(Teuchos::ParameterList& params);

   //@}

   /** @name Jacobian operations.

     Operations using the Jacobian matrix. 
   */

   //@{

   //! Applies Jacobian to the given input vector and puts the answer in the result.
   /*! 
     Computes 
     \f[ v = J u, \f]
     where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result vector.

     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If the Jacobian \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType 
   applyJacobian(const NOX::Abstract::Vector& input, 
		 NOX::Abstract::Vector& result) const;

   //! Applies Jacobian-Transpose to the given input vector and puts the answer in the result.

   /*!  
     Computes 
     \f[ v = J^T u, \f]
      where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result vector.

     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType 
   applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			  NOX::Abstract::Vector& result) const;

   /*!
     \brief Applies the inverse of the Jacobian matrix to the given
     input vector and puts the answer in result.

     Computes 
     \f[ v = J^{-1} u, \f]
     where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result vector.

     The "Tolerance" parameter specifies that the
     solution should be such that
     \f[
     \frac{\| J v - u \|_2}{\max \{ 1, \|u\|_2\} } < \mbox{Tolerance}
     \f]

     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance" 
	  specified in \c params
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>

     The parameter "Tolerance" may be added/modified in the list of
     parameters - this is the ideal solution tolerance for an iterative
     linear solve.
   */
   virtual NOX::Abstract::Group::ReturnType 
   applyJacobianInverse(Teuchos::ParameterList& params, 
			const NOX::Abstract::Vector& input, 
			NOX::Abstract::Vector& result) const;

   //! Apply right preconditiong to the given input vector
   /*!
     Let \f$M\f$ be a right preconditioner for the Jacobian \f$J\f$; in
     other words, \f$M\f$ is a matrix such that
     \f[ JM \approx I. \f]

     Compute
     \f[ u = M^{-1} v, \f]
     where \f$u\f$ is the input vector and \f$v\f$ is the result vector.

     If <em>useTranspose</em> is true, then the transpose of the 
     preconditioner is applied:  
     \f[ u = {M^{-1}}^T v, \f]
     The transpose preconditioner is currently only required for 
     Tensor methods. 

     The "Tolerance" parameter specifies that the
     solution should be such that
     \f[
     \frac{\| M v - u \|_2}{\max \{ 1, \|u\|_2\} } < \mbox{Tolerance}
     \f]

     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance" 
	  specified in \c params
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>

     The parameters are from the "Linear %Solver" sublist of the
     "Direction" sublist that is passed to solver during construction.
   */
   virtual NOX::Abstract::Group::ReturnType 
   applyRightPreconditioning(bool useTranspose,
			     Teuchos::ParameterList& params,
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;
  //@}

   /** @name Block Jacobian operations.

   Operations using the Jacobian matrix. 
   */

   //@{

   //! applyJacobian for multiple right-hand sides 
   /*! 
    * The default implementation here calls applyJacobian() for
    * each right hand side serially but should be overloaded if a 
    * block method is available.
    */
   virtual NOX::Abstract::Group::ReturnType
   applyJacobianMultiVector(const NOX::Abstract::MultiVector& input, 
			    NOX::Abstract::MultiVector& result) const;

   //! applyJacobianTranspose for multiple right-hand sides 
   /*! 
    * The default implementation here calls applyJacobianTranspose() for
    * each right hand side serially but should be overloaded if a 
    * block method is available.
    */
   virtual NOX::Abstract::Group::ReturnType
   applyJacobianTransposeMultiVector(const NOX::Abstract::MultiVector& input, 
				     NOX::Abstract::MultiVector& result) const;

   //! applyJacobianInverse for multiple right-hand sides
   /*! 
    * The default implementation here calls applyJacobianInverse() for
    * each right hand side serially but should be overloaded if a 
    * block solver is available.
    */
   virtual NOX::Abstract::Group::ReturnType 
   applyJacobianInverseMultiVector(Teuchos::ParameterList& params, 
				   const NOX::Abstract::MultiVector& input, 
				   NOX::Abstract::MultiVector& result) const;
   
   //! applyRightPreconditioning for multiple right-hand sides
   /*! 
    * The default implementation here calls applyRightPreconditioning() for
    * each right hand side serially but should be overloaded if a 
    * block method is available.
    */
   virtual NOX::Abstract::Group::ReturnType 
   applyRightPreconditioningMultiVector(
				    bool useTranspose,
				    Teuchos::ParameterList& params,
				    const NOX::Abstract::MultiVector& input, 
				    NOX::Abstract::MultiVector& result) const;

   //@}

  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns true
    if the corresponding "compute" function has been called since the
    last change to the solution vector.
  */

  //@{
  
  //! Return true if F is valid.
  virtual bool isF() const = 0;

  //! Return true if the Jacobian is valid.
  /*! \note Default implementation in NOX::Abstract::Group returns false. */
  virtual bool isJacobian() const;

  //! Return true if the gradient is valid.
  /*! \note Default implementation in NOX::Abstract::Group returns false. */
  virtual bool isGradient() const;

  //! Return true if the Newton direction is valid.
  /*! \note Default implementation in NOX::Abstract::Group returns false. */
  virtual bool isNewton() const;
  //@}

  /** @name "Get" functions.
   
    Note that these function do not check whether or not the vectors
    are valid. Must use the "Is" functions for that purpose.
  */
  //@{ 

  //! Return solution vector.  
  virtual const NOX::Abstract::Vector& getX() const = 0;

  //! Return F(x)
  virtual const NOX::Abstract::Vector& getF() const = 0;

  //! Return 2-norm of F(x). 
  /*! In other words, \f[ \sqrt{\sum_{i=1}^n F_i^2} \f] */
  virtual double getNormF() const = 0;

  //! Return gradient.
  virtual const NOX::Abstract::Vector& getGradient() const = 0;

  //! Return Newton direction.
  virtual const NOX::Abstract::Vector& getNewton() const = 0;

  //@{ \name PyTrilinos Interface functions (require RCPs

  //! Return RCP to solution vector.  
   virtual Teuchos::RCP< const NOX::Abstract::Vector > getXPtr() const
   {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			"NOX::Abstract::Group::getXPtr() not implemented!");
     return Teuchos::null;
   }

  //! Return RCP to F(x)
   virtual Teuchos::RCP< const NOX::Abstract::Vector > getFPtr() const
   {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			"NOX::Abstract::Group::getFPtr() not implemented!");
     return Teuchos::null;
   }

  //! Return RCP to gradient.
   virtual Teuchos::RCP< const NOX::Abstract::Vector > getGradientPtr() const
   {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			"NOX::Abstract::Group::getGradientPtr() not implemented!");
     return Teuchos::null;
   }

  //! Return RCP to Newton direction.
   virtual Teuchos::RCP< const NOX::Abstract::Vector > getNewtonPtr() const
   {
     TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
			"NOX::Abstract::Group::getNewtonPtr() not implemented!");
     return Teuchos::null;
   }

  //@}

  /*!
    \brief 
    Return the norm of the last linear solve residual as the result of
    either a call to computeNewton() or applyJacobianInverse().
  
     \return
     <ul>
     <li> NOX::Abstract::Group::NotDefined - Returned by default implementation in NOX::Abstract::Group
     <li> NOX::Abstract::Group::BadDependency - If no linear solve has been calculated
     <li> NOX::Abstract::Group::Failed - Any other type of failure
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
  */
   virtual NOX::Abstract::Group::ReturnType 
   getNormLastLinearSolveResidual(double& residual) const;

  //@}


  //@{ \name Creating new Groups.

  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a ref count pointer to the new group.  

    If type is NOX::DeepCopy, then we need to create an exact replica
    of "this". Otherwise, if type is NOX::ShapeCopy, we need only
    replicate the shape of "this" (only the memory is allocated, the
    values are not copied into the vectors and Jacobian). Returns NULL
    if clone is not supported.

    \note Any <em>shared data</em> should have its ownership transfered to this
    group from the source for a NOX::DeepCopy.
  */
   virtual Teuchos::RCP<NOX::Abstract::Group> 
   clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

  //@}

private:

};
} // namespace Abstract
} // namespace NOX

#endif
