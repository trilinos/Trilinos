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
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov), Sandia National Laboratories.
// 
// ************************************************************************
//@HEADER

#ifndef LOCA_CONTINUATION_ABSTRACTGROUP_H
#define LOCA_CONTINUATION_ABSTRACTGROUP_H

#include "NOX_Abstract_Group.H"    // base class

namespace NOX {
  namespace Parameter {
    class List;
  }
}

namespace LOCA {
  class ParameterVector;
}

/*! 
 * \brief Library of continuation algorithms namespace.
 */
/*! 
 * All classes, functions, and global variables written for %LOCA reside
 * in the LOCA namespace.
*/
namespace LOCA { 

  //! Groups and vectors for continuation
  namespace Continuation { 

    /*! 
     * \brief %LOCA abstract interface for continuation, derived from the 
     * NOX::Abstract::Group.  This abstract class provides 
     * the interface necessary to perform continuation, i.e., compute 
     * families of solutions to \f$ F(x,p) = 0 \f$.
    */
    /*! 
     * Concrete implemenations of this interface must provide implementations 
     * of all of the methods in the NOX::Abstract::Group interface as well as 
     * the additional interface defined here.
    */

    class AbstractGroup : public virtual NOX::Abstract::Group {

    public:
  
      //! Default constructor.
      AbstractGroup() {}
 
      //! Destructor.
      virtual ~AbstractGroup() {}

      /*! 
       * @name Pure virtual methods 
       * These methods must be defined by any concrete implementation
       */
      //@{

      //! Assignment operator
      virtual AbstractGroup& operator=(const AbstractGroup& source) = 0;

      //! Set the parameter vector in the group to p (pVector = p).  
      virtual void setParams(const ParameterVector& p) = 0;

       //! Set parameter indexed by (integer) paramID
      virtual void setParam(int paramID, double val) = 0;

      //! Set parameter indexed by (string) paramID
      virtual void setParam(string paramID, double val) = 0;

      //! Return a const reference to the ParameterVector owned by the group. 
      virtual const ParameterVector& getParams() const = 0;

      //! Return copy of parameter indexed by (integer) paramID
      virtual double getParam(int paramID) const = 0;

      //! Return copy of parameter indexed by (string) paramID
      virtual double getParam(string paramID) const = 0;

      /*! 
       * Compute \f$\partial F/\partial p\f$ where \f$p\f$ is the 
       * parameter indexed by paramID.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDp(int paramID, NOX::Abstract::Vector& result) = 0;

      //@}

      
      /*! 
       * @name Virtual methods with default implementations 
       * These methods should be overloaded in a concrete implementation
       * if more appropriate/efficient approaches are available.
       */
      //@{

      //! applyJacobianInverse for multiple right-hand sides
      /*! 
       * The default implementation here calls applyJacobianInverse() for
       * each right hand side serially but should be overloaded if a 
       * block solver is available.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyJacobianInverseMulti(NOX::Parameter::List& params,
			const NOX::Abstract::Vector* const* inputs,
			NOX::Abstract::Vector** outputs, int nVecs) const;

      //! Compute a scaled dot product
      /*! 
       * The default implementation here just computes a.dot(b) but should
       * be overloaded for any problem that his difficult scaling.
       */
      virtual double
      computeScaledDotProduct(const NOX::Abstract::Vector& a,
			      const NOX::Abstract::Vector& b) const;

      //@}

      /*! 
       * @name Virtual methods with empty or trivial implementations 
       * These methods should be overloaded in a concrete implementation
       * but their implementation is not critical to the rest of %LOCA and
       * therefore have empty or trivial implementations.
       */
      //@{

      //! Compute eigenvalues/eigenvectors
      /*! 
       * Default implementation prints an error message and throws an 
       * error string.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeEigenvalues(NOX::Parameter::List& params);

      //! Function to print out solution and parameter after successful step
      /*! 
       * Empty default definition. 
       */
      virtual void printSolution(const double conParam) const {};

      //! Function to print out a vector and parameter after successful step
      /*! 
       * Empty default definition. 
       */
      virtual void printSolution(const NOX::Abstract::Vector& x_,
                                 const double conParam) const  {}

      //! Solve arclength continuation equations using Householder projection
      /*!
       * This method should compute a solution to the following system
       * of equations
       * \f[
       *    \begin{bmatrix}
       *      J & \frac{df}{dp}
       *   \end{bmatrix} Q y = f 
       * \f]
       * where \f$Q\in\Re^{n+1\times n}\f$ is an orthogonal matrix such that
       * for any \f$z\in\Re^n\f$,
       * \f[
       *    Q z = P \begin{bmatrix} z \\ 0 \end{bmatrix}.
       * \f]
       * Here \f$P = I - \beta u u^T\in\Re^{n+1\times n+1}\f$ is a Householder
       * transformation such that \f$P v = \|v\|_2 e_{n+1}\f$ and 
       * \f$u_{n+1} = 1\f$ where \f$v\in\Re^{n+1}\f$ is the predictor 
       * vector and \f$e_{n+1}\f$ is the \f$n+1\f$ column of the identity.
       *
       * The default implementation here is to throw an error.  This method
       * is only called by the Householder continuation group, and should
       * be implemented to use this continuation method.
       */
      virtual NOX::Abstract::Group::ReturnType
      applyHouseholderJacobianInverse(NOX::Parameter::List& params,
				      const NOX::Abstract::Vector& f,
				      const NOX::Abstract::Vector& dfdp,
				      const NOX::Abstract::Vector& ux,
				      double up, double beta,
				      NOX::Abstract::Vector& result_x,
				      double& result_p) const;

      //! Scales a vector using scaling vector
      /*!
       * The default definition here is to do nothing, i.e., no scaling
       */
      virtual void
      scaleVector(NOX::Abstract::Vector& x) const;

      //@}

    };
  } // namespace Continuation
} // namespace LOCA

#endif
