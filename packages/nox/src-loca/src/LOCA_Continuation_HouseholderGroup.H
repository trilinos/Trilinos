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

#ifndef LOCA_CONTINUATION_HOUSEHOLDERGROUP_H
#define LOCA_CONTINUATION_HOUSEHOLDERGROUP_H

#include "LOCA_Continuation_ArcLengthGroup.H" // base class

namespace LOCA { 

  namespace Continuation { 

    /*! 
     * \brief Arc-length continuation group using Householder projection
     */
    /*!
     * This group provides an alternative method for computing Newton steps
     * for the pseudo-arclength equations
     * \f[
     *    \begin{bmatrix}
     *       J & \frac{df}{dp} \\
     *       v_x^T & v_p
     *    \end{bmatrix}
     *    \begin{bmatrix}
     *       \Delta x \\
     *       \Delta p
     *    \end{bmatrix} = 
     *    \begin{bmatrix}
     *       -f \\
     *        0
     *    \end{bmatrix}
     * \f]
     * via Householder projection.  In particular, these equations are
     * replaced by
     * \f[
     *    \begin{bmatrix}
     *      J & \frac{df}{dp}
     *   \end{bmatrix} Q y = f 
     * \f]
     * where \f$Q\in\Re^{n+1\times n}$\f$ is an orthogonal matrix such that
     * for any \f$z\in\Re^n\f$,
     * \f[
     *    Q z = P \begin{bmatrix} z \\ 0 \end{bmatrix}.
     * \f]
     * Here \f$P = I - \beta u u^T\in\Re^{n+1\times n+1\f$ is a Householder
     * transformation such that \f$P v = \|v\|_2 e_{n+1}\f$ and 
     * \f$u_{n+1} = 1\f$ where \f$v\in\Re^{n+1}\f$ is the predictor 
     * vector and \f$e_{n+1}\f$ is the \f$n+1\f$ column of the identity.
     * This group takes care of setting up the linear system of equations
     * and computing the Householder vector \f$u\f$, but it is assumed
     * the underlying group provides some means of solving the linear system.
     *
     * See Homer F. Walker, "An Adaptation of Krylov Subspace Methods to 
     * Path Following Problems", SIAM J. Sci. Comput. 21(3), 1999, 
     * pp. 1191-1198, for a description of this method.
     */
    class HouseholderGroup : public LOCA::Continuation::ArcLengthGroup {

    public:
  
      //! Constructor with integer continuation parameter id
      HouseholderGroup(LOCA::Continuation::AbstractGroup& g, int paramID,
		       NOX::Parameter::List& params);

      //! Constructor with string continuation parameter id
      HouseholderGroup(LOCA::Continuation::AbstractGroup& g, string paramID,
		       NOX::Parameter::List& params);

      //! Copy constructor.
      HouseholderGroup(const HouseholderGroup& source, 
		       NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~HouseholderGroup();

      //! Assignment operator
      virtual LOCA::Continuation::ArcLengthGroup& 
      operator=(const LOCA::Continuation::ArcLengthGroup& source);
  
      //! Assignment operator
      virtual LOCA::Continuation::ExtendedGroup& 
      operator=(const LOCA::Continuation::ExtendedGroup& source);

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Assignment operator
      virtual HouseholderGroup& 
      operator=(const HouseholderGroup& source);
     
      //! Cloning function
      virtual NOX::Abstract::Group* 
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      //@{ \name "Compute" functions.

      //! Overload of computeF for Householder projection
      /*!
       * Since the Householder projection method does not use the
       * arc-length equation in the Newton solve, we set the residual
       * for the arc-length equation to zero.
       */
      virtual NOX::Abstract::Group::ReturnType computeF();

      //! Overload of computeNewton to use Householder projection.
      virtual NOX::Abstract::Group::ReturnType 
      computeNewton(NOX::Parameter::List& params);
	
      //@}

    protected:

      //! Computes the Householder vector from the predictor vector
      /*!
       * This method computes a vector \f$u\in\Re^{n+1}\f$ and scalar 
       * \f$\beta\f$ such that \f$P=I-\beta u u^T\f$ is orthogonal,
       * \f$u_{n+1} = 1\f$, and \fPv = \|v\|_2 e_{n+1}\f$ where \f$v\f$
       * is the predictor vector and \f$e_{n+1}\f$ is the \f$n+1\f$
       * column of the identity.
       *
       * The algorithm employed here is taken from Algorithm 5.1.1 from
       * Gene Golub and Charles Van Loan, "Matrix Computations", 3rd Edition,
       * 1996, John Hopkins University Press, Baltimore.
       */
      void computeHouseholderVector();

    protected:
	
      //! Stores the Householder vector
      LOCA::Continuation::ExtendedVector houseVec;

      //! Stores scalar beta in Householder projection
      double beta;

    };
  } // namespace Continuation
} // namespace LOCA

#endif
