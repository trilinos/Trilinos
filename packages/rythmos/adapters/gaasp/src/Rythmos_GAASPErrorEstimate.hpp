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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef Rythmos_GAASP_ERROR_ESTIMATE_H
#define Rythmos_GAASP_ERROR_ESTIMATE_H

#include "Rythmos_ErrorEstimateBase.hpp"

namespace Rythmos {

class GAASPErrorEstimate : public virtual ErrorEstimateBase<double> {
  public:
    
    // Constructor
    GAASPErrorEstimate();

    // Redefined from Rythmos::ErrorEstimateBase
    double getTotalError() const;
    
    // GAASP::ErrorEstimate specific functions:
    void setErrorEstimate(double errorEstimate);
    void setIntervalErrorContributions(double **intError);

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

    // Get valid parameter list
    Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  
  private:
    Teuchos::RCP<Teuchos::ParameterList> paramList_;
    double totalError_;
    double **intervalErrorContributions_;
};

} // namespace Rythmos


#endif // Rythmos_GAASP_ERROR_ESTIMATE_H
