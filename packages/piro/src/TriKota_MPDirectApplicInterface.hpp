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

#ifndef TRIKOTA_MPDIRECTAPPLICINTERFACE
#define TRIKOTA_MPDIRECTAPPLICINTERFACE

// Have to do this first to pull in all of Dakota's #define's
#include "TriKota_ConfigDefs.hpp"

#include "DirectApplicInterface.H"
#include "CommandLineHandler.H"
#include "DakotaStrategy.H"
#include "DakotaModel.H"
#include "ParallelLibrary.H"
#include "ProblemDescDB.H"

#include "Piro_Epetra_StokhosMPSolver.hpp"
#include "Stokhos_ProductEpetraVector.hpp"
#include "Stokhos_ProductEpetraMultiVector.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

//!  TriKota namespace
namespace TriKota {

/*! \brief Adapter class that transates from a Trilinos interface to a 
 * Dakota interface. 
 *
 * This differs from the main TriKota interface in that it wraps multipe
 * evaluations into a single multi-point evaluation.
 */
class MPDirectApplicInterface : public Dakota::DirectApplicInterface {
public:

  //! Constructor for multi-point evaluation
  MPDirectApplicInterface(
    Dakota::ProblemDescDB& problem_db_,
    const Teuchos::RCP<Piro::Epetra::StokhosMPSolver>& model_,
    int p_index = 0, int g_index = 0);

  ~MPDirectApplicInterface() {};

protected:

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  int derived_map_ac(const Dakota::String& ac_name);

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  void derived_map_asynch(const Dakota::ParamResponsePair& pair);

  //! Virtual function redefinition from Dakota::DirectApplicInterface
  //int derived_map_of(const Dakota::String& of_name);

  //int derived_map_if(const Dakota::String& if_name);

  //! evaluate the batch of jobs contained in prp_queue
  void derived_synch(Dakota::PRPQueue& prp_queue);

private:

  // Data
  Teuchos::RCP<Piro::Epetra::StokhosMPSolver> model;
  int p_index;
  int g_index;
  Teuchos::RCP<Teuchos::FancyOStream> out;
  Teuchos::RCP<Stokhos::ProductEpetraVector> model_p;
  Teuchos::RCP<Stokhos::ProductEpetraVector> model_g;
  Teuchos::RCP<Stokhos::ProductEpetraMultiVector> model_dgdp;
  unsigned int numParameters;
  unsigned int numResponses;
  bool supportsSensitivities;
  EpetraExt::ModelEvaluator::EDerivativeMultiVectorOrientation orientation;
  unsigned int block_size;

};

} // namespace TriKota

#endif //TRIKOTA_MPDIRECTAPPLICINTERFACE
