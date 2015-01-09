//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
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
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_GAASP_ERROR_ESTIMATOR_H
#define Rythmos_GAASP_ERROR_ESTIMATOR_H

#include "Rythmos_GAASPInterface.hpp"
#include "Rythmos_GAASPGModel_ThyraModelEvaluator.hpp"
#include "Rythmos_GAASPErrorEstimate.hpp"
#include "Rythmos_ErrorEstimatorBase.hpp"
#include "Rythmos_ErrorEstimateBase.hpp"
#include "GModelBase.h"
#include "InitGaaspOO.h"
#include "GForwardSolve.h"
#include "Teuchos_RCPBoostSharedPtrConversions.hpp"

namespace Rythmos {

class GAASPErrorEstimator : virtual public ErrorEstimatorBase<double> {
public:

  // Destructor
  ~GAASPErrorEstimator() {};

  // Constructor
  GAASPErrorEstimator();

  // Redefined from ErrorEstimatorBase
  // setModel
  void setModel( const Teuchos::RCP<Thyra::ModelEvaluator<double> > model );
  
  // setQuantityOfInterest
  void setQuantityOfInterest( ERROR_QUANTITY_OF_INTEREST qtyOfInterest );

  // getErrorEstimate
  Teuchos::RCP<const ErrorEstimateBase<double> > getErrorEstimate();

  // New functions:

  // Global Error Control
  Teuchos::RCP<const ErrorEstimateBase<double> > controlGlobalError(double uTOL);

  // Get valid parameter list
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

  // Redefined from Teuchos::Describable
  /** \brief . */
  std::string description() const;

  /** \brief . */
  void describe(
    Teuchos::FancyOStream       &out
    ,const Teuchos::EVerbosityLevel      verbLevel
    ) const;

  /// Redefined from Teuchos::ParameterListAcceptor
  /** \brief . */
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const& paramList);

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  
private:
  // initialize 
  void initialize_();

private:
  Teuchos::RCP<Thyra::ModelEvaluator<double> > model_;
  ERROR_QUANTITY_OF_INTEREST qtyOfInterest_;
  Teuchos::RCP<GAASPInterface> gaaspInterfacePtr_;
  bool isInitialized_;
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  static const std::string GAASPInterface_name_;
};


} // namespace Rythmos


#endif // Rythmos_GAASP_ERROR_ESTIMATOR_H
