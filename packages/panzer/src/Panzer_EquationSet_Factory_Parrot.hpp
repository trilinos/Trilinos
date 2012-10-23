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

#ifndef __Panzer_EquationSet_Factory_Parrot_hpp__
#define __Panzer_EquationSet_Factory_Parrot_hpp__ 

#include "Panzer_config.hpp"

#include "Teuchos_RCP.hpp"

#include "Panzer_Traits.hpp"
#include "Panzer_InputEquationSet.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_EquationSet_TemplateManager.hpp"
#include "Panzer_GlobalData.hpp"
#include "Panzer_EquationSet_Parrot.hpp"

#include <boost/unordered_map.hpp>

namespace panzer {

  /** \brief Allocates and initializes an equation set template manager
             based on a privous set of equation set template managers.

     \param ies [in] Input equation set
     \param cell_data [in] The cell data
     \param global_data [in] Global data
     \param build_transient_support [in] If true, the transient evaluators will be built, registered, and required in the Phalanx evaluation graph.

     Returns an RCP to a newly allocated EquationSet_TemplateManager.  
  */
  class EquationSet_FactoryParrot : public EquationSetFactory {
  public:
    EquationSet_FactoryParrot(const std::vector<panzer::InputEquationSet> & inputEqSets,
                              const std::vector< Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> > > & eqSets)
    { 
      TEUCHOS_ASSERT(inputEqSets.size()==eqSets.size()); 
      for(std::size_t i=0;i<inputEqSets.size();i++)
        eqSets_[inputEqSets[i]] = eqSets[i];        
    }

    virtual ~EquationSet_FactoryParrot() {}

    virtual Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >
    buildEquationSet(const panzer::InputEquationSet& ies,
		     const panzer::CellData& cell_data,
		     const Teuchos::RCP<panzer::GlobalData>& global_data,
		     bool build_transient_support) const
    { 
      // lookup input equation set to see if it was built into parrot factory
      HashTable::const_iterator itr = eqSets_.find(ies);
      TEUCHOS_ASSERT(itr!=eqSets_.end());
      Teuchos::RCP<EquationSet_TemplateManager<panzer::Traits> > inEqSet_tm = itr->second;

      // setup builder for template manager using equation set TM that is to be parroted
      EqSetParrot_TemplateBuilder builder(*inEqSet_tm,ies,cell_data,global_data,build_transient_support);

      // build new template manager that contains all parrotted equation set factories
      Teuchos::RCP<EquationSet_TemplateManager<panzer::Traits> > outEqSet_tm 
        = Teuchos::rcp(new EquationSet_TemplateManager<panzer::Traits>);   
      outEqSet_tm->buildObjects(builder);

      return outEqSet_tm;
    }

  private:

    EquationSet_FactoryParrot();    
    EquationSet_FactoryParrot(const EquationSet_FactoryParrot &);    

    // A parrot template builder, that uses a passed in TM to build parrot equations
    class EqSetParrot_TemplateBuilder {                                   
      const panzer::EquationSet_TemplateManager<panzer::Traits> & eqSet_tm_;
      const panzer::InputEquationSet& ies_;
      const panzer::CellData& cell_data_;
      Teuchos::RCP<panzer::GlobalData> global_data_;
      bool build_transient_support_;

    public:
      EqSetParrot_TemplateBuilder(const panzer::EquationSet_TemplateManager<panzer::Traits> & eqSet_tm,
                                  const panzer::InputEquationSet& ies,
                                  const panzer::CellData& cell_data,
                                  const Teuchos::RCP<panzer::GlobalData>& global_data,
                                  bool build_transient_support)
          : eqSet_tm_(eqSet_tm), ies_(ies), cell_data_(cell_data), global_data_(global_data), build_transient_support_(build_transient_support) {}
      
      template<typename EvalT>                                     
      Teuchos::RCP<panzer::EquationSetBase> build() const 
      { return Teuchos::rcp(new EquationSet_Parrot<EvalT>(*eqSet_tm_.getAsBase<EvalT>(),ies_,cell_data_,global_data_,build_transient_support_)); }
    };


    struct IESHash {
       boost::hash<std::string> hash;
       std::size_t operator()(const InputEquationSet & ies) const
       { return hash(ies.name+"_"+ies.prefix); }
    };

    struct IESEquality {
       bool operator()(const InputEquationSet & ies1,const InputEquationSet & ies2) const
       { return (ies1.name==ies2.name) && (ies1.prefix==ies2.prefix); }
    };

    typedef boost::unordered_map<panzer::InputEquationSet,Teuchos::RCP<panzer::EquationSet_TemplateManager<panzer::Traits> >,IESHash,IESEquality> HashTable;

    // hash table lookup for finding correct eq set TM to parrot
    HashTable eqSets_; 
  };

}

#endif
