// @HEADER
// ************************************************************************
// 
//        Piro: Strategy package for embedded analysis capabilitites
//                  Copyright (2010) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include "DirectApplicInterface.hpp"
#include "CommandLineHandler.hpp"
#include "DakotaStrategy.hpp"
#include "DakotaModel.hpp"
#include "ParallelLibrary.hpp"
#include "ProblemDescDB.hpp"

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

  void derived_synch_nowait(Dakota::PRPQueue& prp_queue) { 
    derived_synch(prp_queue); 
  }

  void check_configuration(int max_iterator_concurrency) {}
  void set_communicators_checks(int max_iterator_concurrency) {}

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
