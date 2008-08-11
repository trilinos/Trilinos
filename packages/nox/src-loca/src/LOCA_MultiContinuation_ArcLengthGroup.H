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

#ifndef LOCA_MULTICONTINUATION_ARCLENGTHGROUP_H
#define LOCA_MULTICONTINUATION_ARCLENGTHGROUP_H

#include "LOCA_MultiContinuation_ExtendedGroup.H"       // base class

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Specialization of LOCA::MultiContinuation::ExtendedGroup
     * to pseudo-arclength continuation.
     */
    /*!
     * Pseudo arc-length continuation corresponds to a continuation equation
     * \f$g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=0\f$ with \f$g\f$ given by
     * \f[
     *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
     *         (x-x^\ast)^Tv_x + (p-p^\ast) v_p - \Delta s
     * \f]
     * where \f$v_x\f$, \f$v_p\f$ are the solution and  parameter components 
     * of the predictor direction \f$v\f$ respectively.
     * This corresponds geometrically to constraining the nonlinear solver 
     * steps used in calculating \f$F(x,p)=0\f$ to be orthogonal to the 
     * predictor direction \f$v\f$.  The arclength constraint \f$g\f$ is
     * represented by a LOCA::MultiContinuation::ArcLengthConstraint object.
     *
     * This class also reimplements the scaleTangent() and 
     * computeScaledDotProduct() methods to implement a scaling method 
     * that tries to ensure the solution and parameter contributions to the
     * arc-length equation are of the same order.  Specifically, the
     * arc-length equation is replaced by
     * \f[
     *      (x-x^\ast)^Tv_x + \theta^2(p-p^\ast) v_p - \Delta s = 0
     * \f]
     * where \f$\theta\f$ is chosen so that \f$\theta^2 v_p\f$
     * is equal to a target value, 0.5 by default.  Parameters for this
     * scaling method are passed through the \c continuationParams argument
     * to the constructor and are:
     * <ul>
     * <li> "Enable Arc Length Scaling" -- [bool] (default: true)
     * <li> "Initial Scale Factor" -- [double] (default: 1.0)
     * <li> "Goal Arc Length Parameter Contribution" -- [double] (default: 0.5)
     * <li> "Max Arc Length Parameter Contribution" -- [double] (default: 0.8)
     * <li> "Min Scale Factor" -- [double] (default: 1.0e-3)
     * </ul>
     * Whether this scaling method is used is determined by the
     * "Enable Arc Length Scaling", and the initial value for \f$\theta\f$
     * is given by "Initial Scale Factor".  A new value of \f$\theta\f$ is
     * chosen only if \f$\theta^2 v_p\f$ is larger than the value given by
     * "Max Arc Length Parameter Contribution" and "Min Scale Factor" 
     * provides a minimum value for \f$\theta\f$.
     */
    class ArcLengthGroup : 
      public virtual LOCA::MultiContinuation::ExtendedGroup {

    public:
  
      //! Constructor
      /*!
       * \param global_data [in] Global data object
       * \param topParams [in] Parsed top-level parameter list.
       * \param continuationParams [in] Continuation parameters as described 
       * above.
       * \param grp [in] Group representing \f$F\f$.
       * \param pred [in] Predictor strategy.
       * \param paramIDs [in] Parameter IDs of continuation parameters.
       */
      ArcLengthGroup(
      const Teuchos::RCP<LOCA::GlobalData>& global_data,
      const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
      const Teuchos::RCP<Teuchos::ParameterList>& continuationParams,
      const Teuchos::RCP<LOCA::MultiContinuation::AbstractGroup>& grp,
      const Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>& pred,
      const vector<int>& paramIDs);

      //! Copy constructor
      ArcLengthGroup(const ArcLengthGroup& source, 
		     NOX::CopyType type = NOX::DeepCopy);

      //! Destructor.
      virtual ~ArcLengthGroup();

      /*! 
       * @name Implementation of NOX::Abstract::Group virtual methods 
       */
      //@{

      //! Assignment operator
      virtual NOX::Abstract::Group& 
      operator=(const NOX::Abstract::Group& source);

      //! Clone function
      virtual Teuchos::RCP<NOX::Abstract::Group>
      clone(NOX::CopyType type = NOX::DeepCopy) const;
      
      //@}

      /*! 
       * @name Implementation of LOCA::MultiContinuation::AbstractStrategy virtual methods 
       */
      //@{

      //! Copy
      virtual void copy(const NOX::Abstract::Group& source);

      //! Scales predictor
      virtual void scaleTangent();

      //! Computes a scaled dot product between two continuation vectors
      virtual double computeScaledDotProduct(
			       const NOX::Abstract::Vector& x, 
			       const NOX::Abstract::Vector& y) const;

      //! Calculates scale factors
      virtual void recalculateScaleFactor(double dpds, double thetaOld,
					  double& thetaNew);

      //@}

    private:

      //! Prohibit generation and use of operator=()
      ArcLengthGroup& operator=(const ArcLengthGroup& source);

    protected:

      //! Stores scaling factor for each arclength equation
      vector<double> theta;

      //! Flag indicating whether to do arc-length scaling
      bool doArcLengthScaling;

      //! Goal value of dp/ds squared
      double gGoal;

      //! Minimum value for dp/ds for which rescaling is applied
      double gMax;

      //! Maximum value of scale factor
      double thetaMin;

      //! Flag indicating whether this is the first rescaling of predictor
      bool isFirstRescale;
    };

  } // namespace MultiContinuation
} // namespace LOCA

#endif
