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
 // This program is free software; you can redistribute it and/or modify
 // it under the terms of the GNU General Public License as published by
 // the Free Software Foundation; either version 2, or (at your option)
 // any later version.
 //   
 // This program is distributed in the hope that it will be useful, but
 // WITHOUT ANY WARRANTY; without even the implied warranty of
 // MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 // General Public License for more details.
 //   
 // You should have received a copy of the GNU General Public License
 // along with this program; if not, write to the Free Software
 // Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 // 
 // Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
 // (rppawlo@sandia.gov).
 // 
 // ************************************************************************
 //@HEADER

 #ifndef NOX_ABSTRACT_GROUP_H
 #define NOX_ABSTRACT_GROUP_H

 #include "NOX_Abstract_Vector.H" // for NOX::CopyType
 #include "NOX_Common.H" // for string 

 namespace NOX {
 namespace Parameter {
 class List;
 }
 }

 namespace NOX { 
 namespace Abstract { 

 /*!  \brief %NOX pure abstract interface to a "group"; i.e., a
   solution vector and the corresponding F-vector, Jacobian matrix,
   gradient vector, and Newton vector.

   This class groups together the objects relating to a particular
   solution vector. Each time the solution vector is updated via the
   "computeX" objects, all the previously computed values become
   invalid.

   Recall that our goal is to solve
   \f[
   F(x) = \left[ \begin{array}{c} F_1(x) \\ F_2(x) \\ \vdots \\ F_n(x) \\ \end{array} \right] = 0.
   \f]
   We may also consider the equivalent optimization problem
   \f[
   \min f(x) \equiv \frac{1}{2} \|F(x)\|_2^2.
   \f]

   The \b Jacobian
   is denoted by \f$J(x)\f$ and defined by
   \f[
   J_{ij}(x) = \frac{\partial F_i}{\partial x_j} (x).
   \f]
   The \b Newton \b direction is the solution, \f$s\f$ of
   \f[
   J s = -F.
   \f]
   The \b gradient (of \f$f\f$) is defined as 
   \f[
   g \equiv J^T F.
   \f]
 */

 class Group {

 public:

   //! Return types for group object computations
   enum ReturnType {

     //! Computation completed successfully
     Ok,

     //! Function not defined in derived object
     NotDefined,

     //! Data dependencies not satisfied
     BadDependency,

     //! Unable to satisfy convergence criteria
     NotConverged,

     //! Any other type of failure
     Failed

   };

   //! Constructor.
   /*!  Constructors for derived object should always define a default
     \f$x\f$ value so that getX() is always defined.
   */
   Group() {};

   //! Destructor.
   virtual ~Group() {};

   /*! 
     \brief Copies the values of all vectors and any other data in
     the source group to this group.  

     \note Any <em>shared data</em> should have its ownership transfered to this
     group from the source.
   */
   virtual Group& operator=(const NOX::Abstract::Group& source) = 0;

   //@{ \name "Compute" functions.

   //! Set the solution vector, x, to y.
   /*!
     This should invalidate the function value, Jacobian,
     gradient, and Newton direction.

     Throw an error if the copy fails.
   */
   virtual void setX(const NOX::Abstract::Vector& y) = 0;

   //! Compute and return solution vector, x, where this.x = grp.x + step * d.
   /*! Let \f$x\f$ denote this group's solution and \f$\hat x\f$ denote
     the solution associated with group \c grp. Then we compute 
     \f[ 
     x = \hat x + \mbox{step} \; d. 
     \f]

     This should invalidate the function value, Jacobian, gradient,
     and Newton direction.

     Throw an error if the computation fails.
   */
   virtual void computeX(const NOX::Abstract::Group& grp, 
			 const NOX::Abstract::Vector& d, double step) = 0;

   //! Compute and return F(x). 
   /*! 

     In other words, return \f$F \equiv F(x)\f$. It's generally useful
     to also compute and stores norm of F(x) at this point for later
     access by the getNormF() function.

     \retval
     <ul>
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeF() = 0;

   //! Compute Jacobian.
   /*! 
     In other words, compute \f$  J_{ij} = \frac{\partial F_i}{\partial x_j} (x)\f$.

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeJacobian();

   //! Compute and return gradient.
   //! Returns false if either F or the Jacobian have not been computed.
   /*! 
     In other words, compute \f$g = J^T F\f$. 

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::BadDependency - If either \f$F\f$ or \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeGradient();

   //! Compute and return Newton direction, using parameters for the linear solve.
   /*! Solve \f$Js = -F\f$ for \f$s\f$.  

   The parameters are from the "Linear %Solver" sublist of the
   "%Direction" sublist that is passed to solver during construction.

   The parameter "Tolerance" may be added/modified in the list of
   parameters - this is the ideal solution tolerance for an iterative
   linear solve.

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::BadDependency - If either \f$F\f$ or \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance" 
	  specified in \c params
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType computeNewton(NOX::Parameter::List& params);

   //@}

   /** @name Jacobian operations.

     Operations using the Jacobian matrix. 
   */

   //@{

   //! Applies Jacobian to the given input vector and puts the answer in the result.
   /*! 
     Computes 
     \f[ \mbox{result} = J \; \mbox{input}, \f]
     and returns true.

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed
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
     \f[ \mbox{result} = J^T \; \mbox{input}, \f]
     and returns true.

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::BadDependency - If \f$J\f$ has not been computed
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
   */
   virtual NOX::Abstract::Group::ReturnType 
   applyJacobianTranspose(const NOX::Abstract::Vector& input, 
			  NOX::Abstract::Vector& result) const;

   /*!
     \brief Applies the inverse of the Jacobian matrix to the given input vector.
   */
   /*! 
     Computes 
     \f[ \mbox{result} = J^{-1} \; \mbox{input}, \f]
     and returns true.

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
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
   applyJacobianInverse(NOX::Parameter::List& params, 
			const NOX::Abstract::Vector& input, 
			NOX::Abstract::Vector& result) const;

   //! Apply right preconditiong to the given input vector
   /*!
     Let \f$M\f$ be a right preconditioner for the Jacobian \f$J\f$; in
     other words, \f$M\f$ is a matrix such that
     \f[ JM \approx I. \f]
     Computes 
     \f[ \mbox{result} = M^{-1} \; \mbox{input}, \f]
     and returns true.

     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::NotConverged - If the linear solve fails to satisfy the "Tolerance" 
	  specified in \c params
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>

     The parameters are from the "Linear %Solver" sublist of the
     "%Direction" sublist that is passed to solver during construction.
   */
   virtual NOX::Abstract::Group::ReturnType 
   applyRightPreconditioning(NOX::Parameter::List& params,
			     const NOX::Abstract::Vector& input, 
			     NOX::Abstract::Vector& result) const;
  //@}

  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns true
    if the corresponding "compute" function has been called since the
    last update to the solution vector (via instantiation or
    computeX).
  */

  //@{
  
  //! Return true if F is valid.
  virtual bool isF() const = 0;

  //! Return true if the Jacobian is valid.
  /*! \note If not defined in derived object, returns false. */
  virtual bool isJacobian() const;

  //! Return true if the gradient is valid.
  /*! \note If not defined in derived object, returns false. */
  virtual bool isGradient() const;

  //! Return true if the Newton direction is valid.
  /*! \note If not defined in derived object, returns false. */
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

  /*!
    \brief 
    Return the norm of the last linear solve residual as the result of
    either a call to computeNewton() or applyJacobianInverse().
  
     \retval
     <ul>
     <li> NOX::Abstract::Group::NotDefined - If not defined in derived object
     <li> NOX::Abstract::Group::BadDependency - If no linear solve has not been calculated
     <li> NOX::Abstract::Group::Failed - If the computation fails
     <li> NOX::Abstract::Group::Ok - Otherwise
     </ul>
  */
   virtual NOX::Abstract::Group::ReturnType 
   getNormLastLinearSolveResidual(double& residual) const;

  //@}


  //@{ \name Creating new Groups.

  /*! 
    \brief Create a new %Group of the same derived type as this one by
    cloning this one, and return a pointer to the new group.  

    If type is NOX::DeepCopy, then we need to create an exact replica of
    "this". Otherwise, if type is NOX::ShapeCopy, we need only replicate
    the shape of "this". Returns NULL if clone is not supported.

    \note Any shared data should have its ownership transfered to this
    group from the source.
  */
  virtual NOX::Abstract::Group* clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

  //@}

private:

};
} // namespace Abstract
} // namespace NOX

#endif
