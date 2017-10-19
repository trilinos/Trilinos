// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_ASSEMBLY_ENGINE_INARGS_HPP
#define PANZER_ASSEMBLY_ENGINE_INARGS_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "Panzer_GlobalEvaluationDataContainer.hpp"

class Epetra_Vector;
class Epetra_CrsMatrix;
class Epetra_Map;

namespace panzer {

  class LinearObjContainer;

  class AssemblyEngineInArgs {
  public:

    AssemblyEngineInArgs(const Teuchos::RCP<panzer::LinearObjContainer> & ghostedContainer,
                         const Teuchos::RCP<panzer::LinearObjContainer> & container)
       : ghostedContainer_(ghostedContainer), container_(container) 
       , alpha(Teuchos::ScalarTraits<double>::nan())   // also setup some painful and
       , beta(Teuchos::ScalarTraits<double>::nan())    // hopefully loud initial values
       , time(Teuchos::ScalarTraits<double>::nan())
       , step_size(Teuchos::ScalarTraits<double>::nan())
       , stage_number(Teuchos::ScalarTraits<double>::one())
       , evaluate_transient_terms(false)
       , first_sensitivities_name("")
       , second_sensitivities_name("")
       , apply_dirichlet_beta(false)
       , dirichlet_beta(0.0)
    { }

    AssemblyEngineInArgs()
       : ghostedContainer_(Teuchos::null), container_(Teuchos::null) 
       , alpha(Teuchos::ScalarTraits<double>::nan())   // also setup some painful and
       , beta(Teuchos::ScalarTraits<double>::nan())    // hopefully loud initial values
       , time(Teuchos::ScalarTraits<double>::nan())
       , step_size(Teuchos::ScalarTraits<double>::nan())
       , stage_number(Teuchos::ScalarTraits<double>::one())
       , evaluate_transient_terms(false)
       , first_sensitivities_name("")
       , second_sensitivities_name("")
       , apply_dirichlet_beta(false)
       , dirichlet_beta(0.0)
    { }

    Teuchos::RCP<panzer::LinearObjContainer> ghostedContainer_;
    Teuchos::RCP<panzer::LinearObjContainer> container_;

    double alpha;
    double beta;
    double time;
    double step_size;
    double stage_number;
    std::vector<double> gather_seeds; // generic gather seeds
    bool evaluate_transient_terms;
    std::string first_sensitivities_name;
    std::string second_sensitivities_name;

    bool apply_dirichlet_beta;
    double dirichlet_beta;

    /** Add a global evaluation data object to be used in all FieldManager
      * evaluate calls.
      *
      * \param[in] key Key to pair with global evaluation data object
      * \param[in] ged Pointer to global evaluation data object
      */
    void addGlobalEvaluationData(const std::string & key,const Teuchos::RCP<GlobalEvaluationData> & ged)
    {
       TEUCHOS_TEST_FOR_EXCEPTION(ged_map.find(key)!=ged_map.end(),std::logic_error,
                                  "AssemblyEngine::addGlobalEvaluationData: Method cannot over write existing "
                                  "data object with key \"" + key + "\"");
  
       ged_map[key] = ged;
    }

    void addGlobalEvaluationData(const GlobalEvaluationDataContainer & gedc)
    {
       ged_map.insert(gedc.begin(),gedc.end()); 
    }

    //! Using internal map fill the global evaluation data container object
    void fillGlobalEvaluationDataContainer(GlobalEvaluationDataContainer & gedc) const
    {
       std::map<std::string,Teuchos::RCP<GlobalEvaluationData> >::const_iterator itr;
       for(itr=ged_map.begin();itr!=ged_map.end();++itr)
          gedc.addDataObject(itr->first,itr->second);
    }

    const std::map<std::string,Teuchos::RCP<GlobalEvaluationData> > &
    getGlobalEvaluationDataMap() const
    { return ged_map; }

  private:
    std::map<std::string,Teuchos::RCP<GlobalEvaluationData> > ged_map;
  };


  /** Streaming output for the AssemblyEngineInArgs. This is helpful
    * in debugging exactly what is in the input arguments data structure.
    */
  inline std::ostream & operator<<(std::ostream & os,const AssemblyEngineInArgs & in)
  {
    os << "AE Inargs:\n"
       << "  alpha         = " << in.alpha << "\n"
       << "  beta          = "  << in.beta << "\n"
       << "  time          = "  << in.time << "\n"
       << "  step_size     = "  << in.step_size << "\n"
       << "  stage_number  = "  << in.stage_number << "\n"
       << "  eval_tran     = " << in.evaluate_transient_terms << "\n"
       << "  1st sens_name = " << in.first_sensitivities_name << "\n"
       << "  2nd sens_name = " << in.second_sensitivities_name << "\n"
       << "  apply_db      = " << in.apply_dirichlet_beta << "\n"
       << "  db            = " << in.dirichlet_beta << "\n"
       << "  seeds         = ";
    for(std::size_t i=0;i<in.gather_seeds.size();i++)
      os << in.gather_seeds[i] << " ";
    os << "\n";

    const std::map<std::string,Teuchos::RCP<GlobalEvaluationData> > & ged_map
        = in.getGlobalEvaluationDataMap();
    os << "  ged_map   = ";
    for(std::map<std::string,Teuchos::RCP<GlobalEvaluationData> >::const_iterator itr=ged_map.begin();
        itr!=ged_map.end();++itr) {
      os << "              \"" << itr->first << "\": ";
      itr->second->print(os);
      os << "\n";
    }
    os << std::endl;

    return os;
  }

}

#endif
