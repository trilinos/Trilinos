#include "Panzer_ResponseLibraryDef.hpp"

#include "Teuchos_ParameterEntry.hpp"

namespace panzer {

ResponseLibrary::ResponseLibrary() 
   : previousCheckout_(false)
{} // do nothing constructor

////////////////////////////////////
// Volume and BC methods
void ResponseLibrary::
checkoutResponses(const Teuchos::ParameterList & p)
{
   using Teuchos::ParameterList;
   using Teuchos::ParameterEntry;

   TEST_FOR_EXCEPTION(previousCheckout_,std::logic_error,
                      "Cannot call ResponseLibrary checkoutResponses more than once!");

   ParameterList::ConstIterator itr;
   for(itr=p.begin();itr!=p.end();++itr) {
      const std::string & name = itr->first;
      // const ParameterEntry & entry = itr->second; // this is meaningless now!

      // first search volume elements
      VolumeMap::const_iterator nameBlockPair = volumeResponses_.find(name);

      if(nameBlockPair!=volumeResponses_.end()) {
         const std::set<std::string> & elementBlocks = nameBlockPair->second;

         // insert field tag name into checked out responses map
         for(std::set<std::string>::const_iterator eItr=elementBlocks.begin();
             eItr!=elementBlocks.end();++eItr) 
            coVolumeResponses_[(*eItr)].insert(name);
      }
      
      // must also check boundary condition responses
      TEST_FOR_EXCEPTION(nameBlockPair==volumeResponses_.end(),std::logic_error,
                         "Invalid response requested, cannot find \""+name+"\" in response library."); 
   }

   previousCheckout_ = true;
}

void ResponseLibrary::
printAvailableResponses(std::ostream & out) const
{
   // first print volume responses
   out << "Volume Responses: \n";

   std::vector<std::pair<std::string,std::set<std::string> > > volResponsePairs;

   this->getAvailableVolumeResponses(volResponsePairs);
   for(std::size_t i=0;i<volResponsePairs.size();i++) {
      std::string name = volResponsePairs[i].first;
      const std::set<std::string> & eBlocks = volResponsePairs[i].second;
      
      out << "   Response: \"" << name << "\": Active Blocks = ";
     
      // print active element blocks (where Response is registered)
      std::set<std::string>::const_iterator itr=eBlocks.begin();

      while(itr!=eBlocks.end()) {
         out << "\"" << *(itr++) << "\""; 
 
         // only add comma if this is not the last element block
         if(itr!=eBlocks.end())
            out << ", ";
      }
      out << "\n";
   }

   out.flush(); // force output
}

void ResponseLibrary::
printCheckedOutResponses(std::ostream & out) const
{
   // first print volume responses
   out << "Volume Responses: \n";

   for(VolumeMap::const_iterator itr=coVolumeResponses_.begin();
       itr!=coVolumeResponses_.end();++itr) {
      const std::string & eBlock = itr->first;
      const std::set<std::string> & names = itr->second;
      
      out << "   Element Block: \"" << eBlock << "\": Responses = ";
     
      // print active element blocks (where Response is registered)
      std::set<std::string>::const_iterator nmItr=names.begin();
      while(nmItr!=names.end()) {
         out << "\"" << *(nmItr++) << "\""; 
 
         // only add comma if this is not the last element block
         if(nmItr!=names.end())
            out << ", ";
      }
      out << "\n";
   }

   out.flush(); // force output
}

void ResponseLibrary::
print(std::ostream & out) const
{
   Teuchos::FancyOStream fout(Teuchos::rcpFromRef(out));

   fout << "Available Responses\n";
   fout.pushTab(3);
   this->printAvailableResponses(fout);
   fout.popTab();

   fout << "\nChecked out Responses \n";
   fout.pushTab(3);
   this->printCheckedOutResponses(fout);
   fout.popTab();

   out.flush(); // force output
}

////////////////////////////////////
// Volume methods
void ResponseLibrary::
addResponse(const PHX::FieldTag & ft, const std::string & blockId)
{
   volumeResponses_[ft.name()].insert(blockId);
}

void ResponseLibrary::
getAvailableVolumeResponses(std::vector<std::string> & name) const
{
   name.clear();

   // a simple loop transcribing responses
   for(VolumeMap::const_iterator itr=volumeResponses_.begin();
       itr!=volumeResponses_.end();++itr)
      name.push_back(itr->first);
}

void ResponseLibrary::
getAvailableVolumeResponses(std::vector<std::pair<std::string,std::set<std::string> > > & responsePairs) const
{
   responsePairs.clear();

   // a simple loop transcribing responses
   for(VolumeMap::const_iterator itr=volumeResponses_.begin();
       itr!=volumeResponses_.end();++itr)
      responsePairs.push_back(std::make_pair(itr->first,itr->second));
}

void ResponseLibrary::
getCheckedOutVolumeResponses(std::vector<std::string> & names) const
{
   names.clear();

   std::set<std::string> finalSet;

   for(VolumeMap::const_iterator itr=coVolumeResponses_.begin();
       itr!=coVolumeResponses_.end();++itr) 
      finalSet.insert(itr->second.begin(),itr->second.end());

   // copy set to names
   names.insert(names.end(),finalSet.begin(),finalSet.end());
}

void ResponseLibrary::
getCheckedOutVolumeResponses(const std::string & eBlock,std::vector<std::string> & names) const
{
   names.clear();

   VolumeMap::const_iterator itr = coVolumeResponses_.find(eBlock); 
   if(itr!=coVolumeResponses_.end())
      names.insert(names.end(),itr->second.begin(),itr->second.end());
}

}
