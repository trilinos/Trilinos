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

#ifdef HAVE_TEKO

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk_classic {

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::ParameterListCallbackBlocked(
                      const Teuchos::RCP<const panzer_stk_classic::STKConnManager<GlobalOrdinalT> > & connManager, 
                      const Teuchos::RCP<const panzer::BlockedDOFManager<int,GlobalOrdinalT> > & blocked_ugi)
   : connManager_(connManager), blocked_ugi_(blocked_ugi)
{
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<Teuchos::ParameterList> 
ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::request(const Teko::RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm)); // design by contract

   // loop over parameter list and set the field by a particular key
   Teuchos::RCP<Teuchos::ParameterList> outputPL = rcp(new Teuchos::ParameterList);
   Teuchos::RCP<const Teuchos::ParameterList> inputPL = rm.getParameterList();
   Teuchos::ParameterList::ConstIterator itr;
   for(itr=inputPL->begin();itr!=inputPL->end();++itr) {
      std::string * str_ptr = 0; // just used as a template specifier
      std::string field = inputPL->entry(itr).getValue(str_ptr);
      setFieldByKey(itr->first,field,*outputPL);
   }

   return outputPL;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
bool ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::handlesRequest(const Teko::RequestMesg & rm)
{
   // check if is a parameter list message, and that the parameter
   // list contains the right fields
   if(rm.getName()=="Parameter List") {
     bool isHandled = true;
     Teuchos::RCP<const Teuchos::ParameterList> pl = rm.getParameterList();
     std::string field;
     if(pl->isType<std::string>("x-coordinates")) {
       field = pl->get<std::string>("x-coordinates");
       if(!isField(field)) {
         return false;
       }
     }
     if(pl->isType<std::string>("y-coordinates")) {
       // we assume that the fields must be the same
       if(field != pl->get<std::string>("y-coordinates")) {
         return false;
       }
     }
     if(pl->isType<std::string>("z-coordinates")) {
       // we assume that the fields must be the same
       if(field != pl->get<std::string>("z-coordinates")) {
         return false;
       }
     }

     return isHandled;
   }
   else return false;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::preRequest(const Teko::RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm)); // design by contract

   const std::string & field = getHandledField(*rm.getParameterList());
   int block = blocked_ugi_->getFieldBlock(blocked_ugi_->getFieldNum(field));

   // empty...nothing to do
   buildArrayToVector(block,field);
   buildCoordinates(field);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::setFieldByKey(const std::string & key,const std::string & field,Teuchos::ParameterList & pl) const
{
   if(key=="x-coordinates") {
      double * x = const_cast<double *>(&getCoordinateByField(0,field)[0]);
      pl.set<double*>(key,x);
   }
   else if(key=="y-coordinates") {
      double * y = const_cast<double *>(&getCoordinateByField(1,field)[0]);
      pl.set<double*>(key,y);
   }
   else if(key=="z-coordinates") {
      double * z = const_cast<double *>(&getCoordinateByField(2,field)[0]);
      pl.set<double*>(key,z);
   }
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                         "ParameterListCallback cannot handle key=\"" << key << "\"");
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::buildArrayToVector(int block,const std::string & field)
{
   if(arrayToVector_[field]==Teuchos::null) {
      Teuchos::RCP<const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > ugi = blocked_ugi_->getFieldDOFManagers()[block];
      arrayToVector_[field] = Teuchos::rcp(new panzer::ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>(ugi));
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::buildCoordinates(const std::string & field)
{
   std::map<std::string,Intrepid::FieldContainer<double> > data;

   Teuchos::RCP<const panzer::IntrepidFieldPattern> fieldPattern = getFieldPattern(field);

   std::vector<std::string> elementBlocks;
   blocked_ugi_->getElementBlockIds(elementBlocks);
   for(std::size_t i=0;i<elementBlocks.size();++i) {
      std::string blockId = elementBlocks[i];
      std::vector<std::size_t> localCellIds;

      // allocate block of data to store coordinates
      Intrepid::FieldContainer<double> & fieldData = data[blockId];
      fieldData.resize(connManager_->getElementBlock(blockId).size(),fieldPattern->numberIds());

      if(fieldPattern->supportsInterpolatoryCoordinates()) {
         // get degree of freedom coordiantes
         connManager_->getDofCoords(blockId,*fieldPattern,localCellIds,fieldData);
      }
      else {
         Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
         out.setOutputToRootOnly(-1);
         out << "WARNING: In ParameterListCallback::buildCoordinates(), the Intrepid::FieldPattern in "
             << "block \"" << blockId << "\" does not support interpolatory coordinates. "
             << "This may be fine if coordinates are not actually needed. However if they are then bad things "
             << "will happen. Enjoy!" << std::endl;

         return;
      }
   }

   Teuchos::RCP<Tpetra::MultiVector<double,int,GlobalOrdinalT,Node> > resultVec 
      = arrayToVector_[field]->template getDataVector<double>(field,data);

   switch(resultVec->getNumVectors()) {
   case 3:
      zcoords_[field].resize(resultVec->getLocalLength()); 
      resultVec->getVector(2)->get1dCopy(Teuchos::arrayViewFromVector(zcoords_[field]));
   case 2:
      ycoords_[field].resize(resultVec->getLocalLength()); 
      resultVec->getVector(1)->get1dCopy(Teuchos::arrayViewFromVector(ycoords_[field]));
   case 1:
      xcoords_[field].resize(resultVec->getLocalLength()); 
      resultVec->getVector(0)->get1dCopy(Teuchos::arrayViewFromVector(xcoords_[field]));
      break;
   default:
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
             "ParameterListCallback::buildCoordinates: Constructed multivector has nonphysical dimensions.");
      break;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
std::string ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::
getHandledField(const Teuchos::ParameterList & pl) const
{
  // because this method assumes handlesRequest is true, this call will succeed
  return pl.get<std::string>("x-coordinates");
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
const std::vector<double> & ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::
getCoordinateByField(int dim,const std::string & field) const
{
  TEUCHOS_ASSERT(dim>=0);
  TEUCHOS_ASSERT(dim<=2);

  // get the coordinate vector you want
  const std::map<std::string,std::vector<double> > * coord;
  if(dim==0) coord = &xcoords_;
  else if(dim==1) coord = &ycoords_;
  else if(dim==2) coord = &zcoords_;

  std::map<std::string,std::vector<double> >::const_iterator itr;
  itr = coord->find(field);

  TEUCHOS_TEST_FOR_EXCEPTION(itr==coord->end(),std::runtime_error,
                      "ParameterListCallbackBlocked::getCoordinateByField: Coordinates for field \"" + field +
                      "\" dimension " << dim << " have not been built!");

  return itr->second;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const panzer::IntrepidFieldPattern> ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>
::getFieldPattern(const std::string & fieldName) const
{
  std::vector<std::string> elementBlocks;
  blocked_ugi_->getElementBlockIds(elementBlocks);

  for(std::size_t e=0;e<elementBlocks.size();e++) {
    std::string blockId = elementBlocks[e];

    if(blocked_ugi_->fieldInBlock(fieldName,blockId))
      return Teuchos::rcp_dynamic_cast<const panzer::IntrepidFieldPattern>(blocked_ugi_->getFieldPattern(blockId,fieldName),true);
  }

  return Teuchos::null;
}


} 

#endif
