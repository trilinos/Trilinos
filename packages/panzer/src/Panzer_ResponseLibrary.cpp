#include "Panzer_ResponseLibrary.hpp"

#include "Teuchos_ParameterEntry.hpp"

namespace panzer {

ResponseLibrary::ResponseLibrary() 
   : isInitialized_(false)
{
} // do nothing constructor

////////////////////////////////////
// Volume and BC methods

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
printReservedResponses(std::ostream & out) const
{
   using Teuchos::RCP;

   // first print volume responses
   out << "Volume Responses: \n";

   for(std::map<std::string,RCP<RespContManager> >::const_iterator itr=rsvdVolResp_.begin();
       itr!=rsvdVolResp_.end();++itr) {
      std::string eBlock = itr->first;
      Teuchos::RCP<const RespContManager> respContainerManager = itr->second;

      out << "   Element Block \"" << eBlock << "\"" << std::endl;
 
      // loop over response types
      for(RespContManager::const_iterator rspItr=respContainerManager->begin();
          rspItr!=respContainerManager->end();++rspItr) {

         // get names of fields and response types
         std::string rspType = rspItr->getResponseType();
         std::vector<std::string> names;
         rspItr->getReserved(names);        

         out << "      " << rspType << ": ";
         // loop over names printing information
         for(std::size_t n=0;n<names.size();n++) 
            out << names[n] << " ";
         out << std::endl;
      }
   }

   out << "ending" << std::endl;
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

   fout << "\nReserved Responses \n";
   fout.pushTab(3);
   this->printReservedResponses(fout);
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
getAvailableVolumeResponses(std::vector<std::pair<std::string,std::set<std::string> > > & responsePairs) const
{
   responsePairs.clear();

   // a simple loop transcribing responses
   for(VolumeMap::const_iterator itr=volumeResponses_.begin();
       itr!=volumeResponses_.end();++itr)
      responsePairs.push_back(std::make_pair(itr->first,itr->second));
}

void ResponseLibrary::
initializeReservedVolumes() 
{
   if(isInitialized_)
      return;

   rsvdVolResp_.clear();

   // extract unique element block names
   std::set<std::string> eBlocks;
   for(VolumeMap::const_iterator itr=volumeResponses_.begin();itr!=volumeResponses_.end();itr++)
      eBlocks.insert(itr->second.begin(),itr->second.end());

   // build response container objects for this element block
   for(std::set<std::string>::const_iterator itr=eBlocks.begin();
       itr!=eBlocks.end();++itr) {
      rsvdVolResp_[*itr] = Teuchos::rcp(new RespContManager);
      rsvdVolResp_[*itr]->buildObjects();
   }

   isInitialized_ = true;
}

}
