 //@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
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

#ifndef NOX_MULTIPHYSICS_GROUP_H
#define NOX_MULTIPHYSICS_GROUP_H

#include "NOX_Abstract_Group.H"         // Base class
#include "NOX_MultiVector.H" 
#include "NOX_Common.H"                 // for string 
#include "NOX.H"                        // for NOX::Solver::Manager
#include "Teuchos_RefCountPtr.hpp" 

namespace NOX {
  namespace Parameter {
    class List;
  }
  namespace Abstract {
    class MultiVector;
  }
}

namespace NOX { 
namespace Multiphysics { 

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

 class Group : public NOX::Abstract::Group {

 public:

   //! Constructor.
   /*!  
     \note Constructors for any derived object should always define a default
     x-value so that getX() is always defined.
   */
   Group( const Teuchos::RefCountPtr< vector<NOX::Solver::Manager*> >& solvers, 
          const Teuchos::RefCountPtr<NOX::StatusTest::Generic>& t, 
          const Teuchos::RefCountPtr<Teuchos::ParameterList>& p);

   //! Constructor.
   Group( const Group & grp, NOX::CopyType typ );

   //! Destructor.
   virtual ~Group();

   /*! 
     \brief Copies the source group into this group.

     \note Any <em>shared data</em> owned by the source should have
     its ownership transfered to this group. This may result in
     a secret modification to the source object.
   */
   virtual NOX::Abstract::Group & operator=(const NOX::Abstract::Group& source);

   virtual NOX::Abstract::Group & operator=(const NOX::Multiphysics::Group& source);

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
   virtual void setX(const NOX::Abstract::Vector& y);

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
			 const NOX::Abstract::Vector& d, double step);

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
   virtual NOX::Abstract::Group::ReturnType computeF();


  /** @name "Is" functions.
   
    Checks to see if various objects have been computed. Returns true
    if the corresponding "compute" function has been called since the
    last change to the solution vector.
  */

  //@{
  
  //! Return true if F is valid.
  virtual bool isF() const;

  /** @name "Get" functions.
   
    Note that these function do not check whether or not the vectors
    are valid. Must use the "Is" functions for that purpose.
  */
  //@{ 

  //! Return solution vector.  
  virtual const NOX::Abstract::Vector& getX() const;

  //! Return F(x)
  virtual const NOX::Abstract::Vector& getF() const;

  //! Return 2-norm of F(x). 
  /*! In other words, \f[ \sqrt{\sum_{i=1}^n F_i^2} \f] */
  virtual double getNormF() const;

  //! Return gradient.
  virtual const NOX::Abstract::Vector& getGradient() const;

  //! Return Newton direction.
  virtual const NOX::Abstract::Vector& getNewton() const;

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
   virtual Teuchos::RefCountPtr<NOX::Abstract::Group> 
   clone(NOX::CopyType type = NOX::DeepCopy) const;

  //@}

protected:

  //! resets the isValid flags to false
  virtual void resetIsValid();

protected:

  //! Pointer to the container of solvers for each problem to be coupled
  Teuchos::RefCountPtr<vector<NOX::Solver::Manager*> > solversVecPtr;

  /** @name IsValid flags 
   *  
   * True if the current solution is up-to-date with respect to the
   * currect xVector. */
  //@{
  bool isValidRHS;

  //! 2-Norm of RHS
  double normRHS;

};
} // namespace Multiphysics
} // namespace NOX

#endif
