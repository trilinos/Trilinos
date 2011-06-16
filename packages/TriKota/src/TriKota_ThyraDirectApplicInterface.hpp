// @HEADER
// ************************************************************************
// 
//        TriKota: A Trilinos Wrapper for the Dakota Framework
//                  Copyright (2009) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov), Sandia
// National Laboratories.
// 
// ************************************************************************
// @HEADER

#ifndef TRIKOTA_DIRECTAPPLICINTERFACE
#define TRIKOTA_DIRECTAPPLICINTERFACE

// Have to do this first to pull in all of Dakota's #define's
#include "TriKota_ConfigDefs.hpp"

#include "DirectApplicInterface.H"
#include "CommandLineHandler.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"

#include "Thyra_ModelEvaluatorDefaultBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_TestForException.hpp"

//!  TriKota namespace
namespace TriKota {

/*! \brief Adapter class that transates from a Trilinos interface toa Dakota interface. 
  An object of this class IS a Dakota::DirectApplicInterface
  and wraps an EpetraExt::ModelEvaluator. It can then be passed in
  as the argument to the TriKota::Driver::run method.
*/
class ThyraDirectApplicInterface : public Dakota::DirectApplicInterface
{
public:

  //! Constructor that takes the Model Evaluator to wrap

   ThyraDirectApplicInterface(Dakota::ProblemDescDB& problem_db_,
                         const Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > App_);
                         //const Teuchos::RCP<Thyra::ModelEvaluator<double> > App_);

  ~ThyraDirectApplicInterface() {};

protected:

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_ac(const Dakota::String& ac_name);

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_of(const Dakota::String& of_name);

  //int derived_map_if(const Dakota::String& if_name);

private:

  // Data
    Teuchos::RCP<Thyra::ModelEvaluatorDefaultBase<double> > App;
    Teuchos::RCP<Thyra::VectorBase<double> > model_p;
    Teuchos::RCP<Thyra::VectorBase<double> > model_g;
    Teuchos::RCP<Thyra::MultiVectorBase<double> > model_dgdp;
    Thyra::ModelEvaluatorBase::EDerivativeMultiVectorOrientation orientation;
    unsigned int numParameters;
    unsigned int numResponses;
    bool supportsSensitivities;
};

} // namespace TriKota

#endif //TRIKOTA_DIRECTAPPLICINTERFACE
