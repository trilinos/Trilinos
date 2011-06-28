#ifdef HAVE_TEKO

using Teuchos::RCP;
using Teuchos::rcp;

namespace panzer_stk {

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
ParameterListCallback<LocalOrdinalT,GlobalOrdinalT,Node>::ParameterListCallback(
                                             const std::string & coordFieldName,
                                             const map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fps,
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
   TEST_FOR_EXCEPTION(!coordinatesBuilt_,std::runtime_error,
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
      TEST_FOR_EXCEPTION(true,std::runtime_error,
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

      // get degree of freedom coordiantes
      connManager_->getDofCoords(blockId,*fieldPattern,localCellIds,fieldData);
   }

   Teuchos::RCP<Tpetra::MultiVector<double,std::size_t,GlobalOrdinalT,Node> > resultVec 
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
      TEST_FOR_EXCEPTION(true,std::logic_error,
             "ParameterListCallback::buildCoordinates: Constructed multivector has nonphysical dimensions.");
      break;
   }

   coordinatesBuilt_ = true;
}

} 

#endif
