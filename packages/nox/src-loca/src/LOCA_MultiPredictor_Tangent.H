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

#ifndef LOCA_MULTIPREDICTOR_TANGENT_H
#define LOCA_MULTIPREDICTOR_TANGENT_H

#include "LOCA_MultiPredictor_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA{
  class GlobalData;
}

namespace LOCA {

  namespace MultiPredictor {

    //! %Tangent predictor strategy
    /*! 
     * This class implements a predictor strategy based on computing 
     * the tangent to the continuation curve.  If \f$p\f$ is the vector
     * of continuation parameters, then the solution component of the 
     * tangent vectors \f$v_x\f$ are computed by solving
     * \f[
     *      J v_x = - \frac{\partial f}{\partial p}.
     * \f]
     * The parameter component \f$v_p\f$ is set to the identity matrix.
     */
    class Tangent : public LOCA::MultiPredictor::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object
       * \param predParams [in] Predictor parameters. Currently no
       * parameters are used by the Tangent predictor.
       * \param solverParams [in] Linear solver parameters used in linear
       * solve to compute tangent vectors \f$v\f$.
       */
      Tangent(const Teuchos::RCP<LOCA::GlobalData>& global_data,
	      const Teuchos::RCP<Teuchos::ParameterList>& predParams,
	      const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

      //! Destructor
      virtual ~Tangent();

      //! Copy constructor
      Tangent(const Tangent& source, NOX::CopyType type = NOX::DeepCopy);
  
      //! Assignment operator
      virtual LOCA::MultiPredictor::AbstractStrategy& 
      operator=(const LOCA::MultiPredictor::AbstractStrategy& source);

      //! Clone function
      virtual Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      /*! 
       * \brief Compute the predictor given the current and
       * previous solution vectors.  Set \em baseOnSecant to \em false
       * if the predictor orientation should not be based on the
       * secant vector (first or last steps of a continuation run).
       */
      /*!
       * This method actually implements the predictor solve described above
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, const vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec);

      //! Evaluate predictor with step size \c stepSize.
      /*!
       * This method computes
       * \c result[i] = \c xVec[i] + \c stepSize[i] * v[i]
       * for each i, where \c v[i] is the ith predictor direction.
       */
      virtual NOX::Abstract::Group::ReturnType 
      evaluate(const vector<double>& stepSize,
	       const LOCA::MultiContinuation::ExtendedVector& xVec,
	       LOCA::MultiContinuation::ExtendedMultiVector& result) const;

      //! Compute tangent to predictor and store in \c tangent.
      virtual NOX::Abstract::Group::ReturnType 
      computeTangent(LOCA::MultiContinuation::ExtendedMultiVector& tangent);

      //! Is the tangent vector for this predictor scalable
      /*!
       * For the tangent predictor, this always returns true.
       */
      virtual bool isTangentScalable() const;

    protected:

      //! LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Stores linear solver parameters
      Teuchos::RCP<Teuchos::ParameterList> linSolverParams;

      //! Stores derivative of residual w.r.t. parameters
      Teuchos::RCP<NOX::Abstract::MultiVector> fdfdp;

      //! Stores tangent vector
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedMultiVector> tangent;

      //! Stores secant vector for setting orientation
      Teuchos::RCP<LOCA::MultiContinuation::ExtendedVector> secant;

      //! Flag indicating whether vectors have been initialized
      bool initialized;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
