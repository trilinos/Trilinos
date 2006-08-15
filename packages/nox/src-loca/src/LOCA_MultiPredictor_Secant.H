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

#ifndef LOCA_MULTIPREDICTOR_SECANT_H
#define LOCA_MULTIPREDICTOR_SECANT_H

#include "LOCA_MultiPredictor_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace MultiPredictor {

    //! %Secant predictor strategy
    /*! 
     * This class implements a predictor strategy based on computing 
     * the secant vector \f$v\f$ to the continuation curve given by
     * \f[
     *      v = x - x_o
     * \f]
     * where \f$x\f$ is the current solution vector and \f$x_o\f$ is
     * the previous solution vector.  Note that for multi-parameter 
     * continuation, the solution component for each secant direction is
     * given as above, with the parameter components given by the
     * identity matrix.
     *
     * For the first step of a continuation run, \f$x_o\f$ is not defined,
     * and so a different predictor is used for this step.  This predictor
     * is specified via the "First Step Predictor" sublist of the "Predictor"
     * sublist.  This predictor is instantiated using the LOCA::Factory 
     * as usual.
     *
     * The parameters used by this class supplied in the constructor are:
     * <ul>
     * <li> "First Step Predictor" Predictor sublist for first step
     * </ul>
     */
    class Secant : public LOCA::MultiPredictor::AbstractStrategy {

    public:

      //! Constructor. 
      /*!
       * \param global_data [in] Global data object used for LOCA::Factory
       * \param topParams [in] Parsed top-level parameter list used when
       * creating first step predictor
       * \param predParams [in] Predictor parameters used to obtain 
       * "First Step Predictor" as described above.
       */
      Secant(
	 const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	 const Teuchos::RefCountPtr<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RefCountPtr<Teuchos::ParameterList>& predParams);

      //! Destructor
      virtual ~Secant();

      //! Copy constructor
      Secant(const Secant& source, NOX::CopyType type = NOX::DeepCopy);
  
      //! Assignment operator
      virtual LOCA::MultiPredictor::AbstractStrategy& 
      operator=(const LOCA::MultiPredictor::AbstractStrategy& source);

      //! Clone function
      virtual Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy>
      clone(NOX::CopyType type = NOX::DeepCopy) const;

      /*! 
       * \brief Compute the predictor given the current and
       * previous solution vectors.  Set \em baseOnSecant to \em false
       * if the predictor orientation should not be based on the
       * secant vector (first or last steps of a continuation run).
       */
      /*!
       * This method actually implements the secant calculation described above
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
       * For the secant predictor, this always returns true.
       */
      virtual bool isTangentScalable() const;

    protected:

      //! LOCA global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! First step predictor
      Teuchos::RefCountPtr<LOCA::MultiPredictor::AbstractStrategy> firstStepPredictor;

      //! Flag indicating whether this is the first step
      bool isFirstStep;

      //! Flag indicating whether first step predictor has been computed
      bool isFirstStepComputed;

      //! Stores predictor vector
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedMultiVector> predictor;

      //! Stores secant vector for setting orientation
      Teuchos::RefCountPtr<LOCA::MultiContinuation::ExtendedVector> secant;

      //! Flag indicating whether vectors have been initialized
      bool initialized;

    };
  } // namespace Predictor
} // namespace LOCA

#endif
