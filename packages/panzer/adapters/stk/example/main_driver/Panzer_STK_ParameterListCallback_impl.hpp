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

namespace panzer_stk {

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::ParameterListCallback(
                                             const std::string & coordFieldName,
                                             const std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fps,
                                             const Teuchos::RCP<const panzer_stk::STKConnManager> & connManager, 
                                             const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi)
   : coordFieldName_(coordFieldName), fieldPatterns_(fps), connManager_(connManager), ugi_(ugi), coordinatesBuilt_(false)
{ }

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<Teuchos::ParameterList> ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::request(const Teko::RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm)); // design by contract

   // loop over parameter list and set the field by a particular key
   Teuchos::RCP<Teuchos::ParameterList> outputPL = rcp(new Teuchos::ParameterList);
   Teuchos::RCP<const Teuchos::ParameterList> inputPL = rm.getParameterList();
   Teuchos::ParameterList::ConstIterator itr;
   for(itr=inputPL->begin();itr!=inputPL->end();++itr)
      setFieldByKey(itr->first,*outputPL);

   return outputPL;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
bool ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::handlesRequest(const Teko::RequestMesg & rm)
{
   // check if is a parameter list message, and that the parameter
   // list contains the right fields
   if(rm.getName()=="Parameter List") return true;
   else return false;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::preRequest(const Teko::RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm)); // design by contract

   // empty...nothing to do
   buildArrayToVector();
   buildCoordinates();
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::setFieldByKey(const std::string & key,Teuchos::ParameterList & pl) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(!coordinatesBuilt_,std::runtime_error,
                      "ParameterListCallback::setFieldByKey: Coordinates have not been built!");

   double * x = const_cast<double *>(&xcoords_[0]);
   double * y = const_cast<double *>(&ycoords_[0]);
   double * z = const_cast<double *>(&zcoords_[0]);

   if(key=="x-coordinates") 
      pl.set<double*>(key,x);
   else if(key=="y-coordinates") 
      pl.set<double*>(key,y);
   else if(key=="z-coordinates") 
      pl.set<double*>(key,z);
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                         "ParameterListCallback cannot handle key=\"" << key << "\"");
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::buildArrayToVector()
{
   if(arrayToVector_==Teuchos::null)
      arrayToVector_ = Teuchos::rcp(new panzer::ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>(ugi_));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::buildCoordinates()
{
   std::map<std::string,Intrepid::FieldContainer<double> > data;

   std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> >::const_iterator itr; 
   for(itr=fieldPatterns_.begin();itr!=fieldPatterns_.end();++itr) {
      std::string blockId = itr->first;
      Teuchos::RCP<const panzer::IntrepidFieldPattern> fieldPattern = itr->second;
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

         coordinatesBuilt_ = true;
         return;
      }
   }

   Teuchos::RCP<Tpetra::MultiVector<double,int,GlobalOrdinalT,Node> > resultVec 
      = arrayToVector_->template getDataVector<double>(coordFieldName_,data);

   switch(resultVec->getNumVectors()) {
   case 3:
      zcoords_.resize(resultVec->getLocalLength()); 
      resultVec->getVector(2)->get1dCopy(Teuchos::arrayViewFromVector(zcoords_));
   case 2:
      ycoords_.resize(resultVec->getLocalLength()); 
      resultVec->getVector(1)->get1dCopy(Teuchos::arrayViewFromVector(ycoords_));
   case 1:
      xcoords_.resize(resultVec->getLocalLength()); 
      resultVec->getVector(0)->get1dCopy(Teuchos::arrayViewFromVector(xcoords_));
      break;
   default:
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
             "ParameterListCallback::buildCoordinates: Constructed multivector has nonphysical dimensions.");
      break;
   }

   coordinatesBuilt_ = true;
}

} 

#endif
