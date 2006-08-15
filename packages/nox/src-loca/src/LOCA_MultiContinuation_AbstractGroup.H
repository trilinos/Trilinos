// $Id$ 
// $Source$ 

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

#ifndef LOCA_MULTICONTINUATION_ABSTRACTGROUP_H
#define LOCA_MULTICONTINUATION_ABSTRACTGROUP_H

#include "NOX_Abstract_Group.H" // base class
#include "NOX_Abstract_MultiVector.H"
#include "LOCA_Abstract_Iterator.H"

namespace LOCA {
  class ParameterVector;
}

namespace LOCA { 

  //! Groups and vectors for multi-parameter continuation
  namespace MultiContinuation { 

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

      //! Copy the group (replaces operator = )
      virtual void copy(const NOX::Abstract::Group& source) = 0;

      //! Set parameters indexed by (integer) paramIDs
      virtual void setParamsMulti(
		     const vector<int>& paramIDs, 
		     const NOX::Abstract::MultiVector::DenseMatrix& vals) = 0;

      //! Set the parameter vector in the group to p (pVector = p).  
      virtual void setParams(const LOCA::ParameterVector& p) = 0;

       //! Set parameter indexed by (integer) paramID
      virtual void setParam(int paramID, double val) = 0;

      //! Set parameter indexed by (string) paramID
      virtual void setParam(string paramID, double val) = 0;

      //! Return a const reference to the ParameterVector owned by the group. 
      virtual const LOCA::ParameterVector& getParams() const = 0;

      //! Return copy of parameter indexed by (integer) paramID
      virtual double getParam(int paramID) const = 0;

      //! Return copy of parameter indexed by (string) paramID
      virtual double getParam(string paramID) const = 0;

      /*! 
       * Compute \f$\partial F/\partial p\f$ for each parameter \f$p\f$ 
       * indexed by paramIDs.  The first column of \em dfdp holds F,
       * which is valid if \em isValidF is true.  Otherwise F must be
       * computed.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDfDpMulti(const vector<int>& paramIDs, 
		       NOX::Abstract::MultiVector& dfdp, 
		       bool isValidF) = 0;

      //@}

      
      /*! 
       * @name Virtual methods with default implementations 
       * These methods should be overloaded in a concrete implementation
       * if more appropriate/efficient approaches are available.
       */
      //@{

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.  The default implementation to empty.
       */
      virtual void 
      preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus);

      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.  The default implementation to empty.
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
       * The default implementation is the max norm of the vector.
       */
      virtual void projectToDraw(const NOX::Abstract::Vector& x,
				 double *px) const;

      //! Returns the dimension of the project to draw array
      /*!
       * The default implementation is to return 1 since the default
       * projection is the max norm of the vector (a scalar).
       */
      virtual int projectToDrawDimension() const;

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

      //! Scales a vector using scaling vector
      /*!
       * The default definition here is to do nothing, i.e., no scaling
       */
      virtual void
      scaleVector(NOX::Abstract::Vector& x) const;

      //@}

      //! Bring NOX::Abstract::Group::operator=() into scope
      using NOX::Abstract::Group::operator=;

    };
  } // namespace MultiContinuation
} // namespace LOCA

#endif
