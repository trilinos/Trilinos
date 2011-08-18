#ifndef __Panzer_ResponseLibrary_hpp__
#define __Panzer_ResponseLibrary_hpp__

#include "Panzer_ResponseLibraryDef.hpp"

#include "Phalanx_FieldTag_Tag.hpp"
#include "Panzer_ScatterResponses.hpp"

namespace panzer {

template <typename Traits,typename EvalT>
void ResponseLibrary::
registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                  const std::string & blockId,PHX::FieldManager<Traits> & fm) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   std::vector<std::string> names;
   this->getCheckedOutVolumeResponses(blockId,names);

   // add scatter residual
   {
      Teuchos::ParameterList p;
      p.set("Comm", comm);
      p.set("Names", Teuchos::rcpFromRef(names));
      p.set("Workset Size", worksetSize);
 
      RCP<panzer::ScatterResponses<EvalT,panzer::Traits> > e 
         = rcp(new panzer::ScatterResponses<EvalT,panzer::Traits>(p));

      fm.template requireField<EvalT>(e->getRequiredFieldTag());
      fm.template registerEvaluator<EvalT>(e);
   }
}

}

#endif
