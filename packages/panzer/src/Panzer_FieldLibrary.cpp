#include "Panzer_FieldLibrary.hpp"

namespace panzer {

Teuchos::RCP<const panzer::PureBasis> FieldLayoutLibrary::lookupBasis(const std::string & fieldName) const
{
   Teuchos::RCP<panzer::BasisIRLayout> layout = lookupLayout(fieldName);

   return layout->getBasis();
}

void FieldLayoutLibrary::uniqueBases(std::list<Teuchos::RCP<const panzer::PureBasis> > & bases) const
{
   bases.clear();
   
   // simply loop over map of basis name to pointers and add them to the list
   std::map<std::string,Teuchos::RCP<const panzer::PureBasis> >::const_iterator itr;
   for(itr=basisNameToPointer_.begin();itr!=basisNameToPointer_.end();++itr) 
      bases.push_back(itr->second);
}

void FieldLayoutLibrary::addFieldAndLayout(const std::string & fieldName,
                                      const Teuchos::RCP<panzer::BasisIRLayout> & layout)
{
   fieldToLayout_[fieldName] = layout; 
   basisNameToPointer_[layout->getBasis()->name()] = layout->getBasis();
}

Teuchos::RCP<panzer::BasisIRLayout> FieldLayoutLibrary::lookupLayout(const std::string & fieldName) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::BasisIRLayout> > Map;
   Map::const_iterator itr = fieldToLayout_.find(fieldName);
   if(itr!=fieldToLayout_.end())
      return itr->second;
 
   return Teuchos::null;
}

void FieldLayoutLibrary::print(std::ostream & os) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::BasisIRLayout> > Map;
   
   for(Map::const_iterator itr=fieldToLayout_.begin();itr!=fieldToLayout_.end();++itr) {
      std::string fieldName = itr->first; 
      Teuchos::RCP<BasisIRLayout> basis = itr->second;

      os << "\"" << fieldName << "\"" << " {" << basis->name() 
         << "(dim=" << basis->getDimension() 
         << ",cells=" << basis->getNumCells() 
         << ",irdeg=" << basis->integrationRuleDegree() << ")} ";
   }
}

Teuchos::RCP<const panzer::PureBasis> FieldLibrary::lookupBasis(const std::string & fieldName) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;
   Map::const_iterator itr = fieldToBasis_.find(fieldName);
   if(itr!=fieldToBasis_.end())
      return itr->second;
 
   return Teuchos::null;
}

void FieldLibrary::uniqueBases(std::list<Teuchos::RCP<const panzer::PureBasis> > & bases) const
{
   bases.clear();
   
   // simply loop over map of basis name to pointers and add them to the list
   std::map<std::string,Teuchos::RCP<const panzer::PureBasis> >::const_iterator itr;
   for(itr=basisNameToPointer_.begin();itr!=basisNameToPointer_.end();++itr) 
      bases.push_back(itr->second);
}

void FieldLibrary::addFieldAndBasis(const std::string & fieldName,
                                    const Teuchos::RCP<panzer::PureBasis> & basis)
{
   fieldToBasis_[fieldName] = basis;
   basisNameToPointer_[basis->name()] = basis;
}

Teuchos::RCP<const FieldLayoutLibrary> FieldLibrary::buildFieldLayoutLibrary(panzer::IntegrationRule & ir) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;

   Teuchos::RCP<FieldLayoutLibrary> layoutLibrary = Teuchos::rcp(new FieldLayoutLibrary);

   // loop over each member of the map, create a new FieldLayout and addit to layout library
   for(Map::const_iterator itr=fieldToBasis_.begin();itr!=fieldToBasis_.end();++itr) {
      Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(itr->second,ir));
      layoutLibrary->addFieldAndLayout(itr->first,layout);
   }
 
   return layoutLibrary;
}

void FieldLibrary::print(std::ostream & os) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;
   
   for(Map::const_iterator itr=fieldToBasis_.begin();itr!=fieldToBasis_.end();++itr) {
      std::string fieldName = itr->first; 
      Teuchos::RCP<PureBasis> basis = itr->second;

      os << "\"" << fieldName << "\"" << " {" << basis->name() 
         << "(dim=" << basis->getDimension() 
         << ",cells=" << basis->getNumCells() << ") ";
   }
}

}
