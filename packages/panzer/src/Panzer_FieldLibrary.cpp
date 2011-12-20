#include "Panzer_FieldLibrary.hpp"

namespace panzer {

void FieldLayoutLibrary::addFieldAndLayout(const std::string & fieldName,
                                      const Teuchos::RCP<panzer::Basis> & layout)
{
   fieldToLayout_[fieldName] = layout; 
}

Teuchos::RCP<panzer::Basis> FieldLayoutLibrary::lookup(const std::string & fieldName) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::Basis> > Map;
   Map::const_iterator itr = fieldToLayout_.find(fieldName);
   if(itr!=fieldToLayout_.end())
      return itr->second;
 
   return Teuchos::null;
}

Teuchos::RCP<panzer::PureBasis> FieldLibrary::lookup(const std::string & fieldName) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;
   Map::const_iterator itr = fieldToBasis_.find(fieldName);
   if(itr!=fieldToBasis_.end())
      return itr->second;
 
   return Teuchos::null;
}

void FieldLibrary::addFieldAndBasis(const std::string & fieldName,
                                    const Teuchos::RCP<panzer::PureBasis> & basis)
{
   fieldToBasis_[fieldName] = basis;
}

Teuchos::RCP<const FieldLayoutLibrary> FieldLibrary::buildFieldLayoutLibrary(panzer::IntegrationRule & ir) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;

   Teuchos::RCP<FieldLayoutLibrary> layoutLibrary = Teuchos::rcp(new FieldLayoutLibrary);

   // loop over each member of the map, create a new FieldLayout and addit to layout library
   for(Map::const_iterator itr=fieldToBasis_.begin();itr!=fieldToBasis_.end();++itr) {
      Teuchos::RCP<Basis> layout = Teuchos::rcp(new Basis(*itr->second,ir));
      layoutLibrary->addFieldAndLayout(itr->first,layout);
   }
 
   return layoutLibrary;
}

}
