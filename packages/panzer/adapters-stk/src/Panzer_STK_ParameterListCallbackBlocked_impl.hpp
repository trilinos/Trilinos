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

#ifdef PANZER_HAVE_TEKO

namespace panzer_stk {

using Teuchos::RCP;
using Teuchos::rcp;

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::ParameterListCallbackBlocked(
                      const Teuchos::RCP<const panzer_stk::STKConnManager<GlobalOrdinalT> > & connManager, 
                      const Teuchos::RCP<const panzer::BlockedDOFManager<int,GlobalOrdinalT> > & blocked_ugi,
                      const Teuchos::RCP<const panzer::BlockedDOFManager<int,GlobalOrdinalT> > & aux_blocked_ugi)
   : connManager_(connManager), blocked_ugi_(blocked_ugi), aux_blocked_ugi_(aux_blocked_ugi)
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
     if(pl->isType<std::string>("Coordinates")){
       field = pl->get<std::string>("Coordinates");
     }
     if(pl->isType<std::string>("Coordinates-Epetra")){
       field = pl->get<std::string>("Coordinates-Epetra");
     }

     return isHandled;
   }
   else return false;
}

template <typename LocalOrdinalT, typename GlobalOrdinalT, typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT, GlobalOrdinalT, Node>::
preRequest(const Teko::RequestMesg & rm)
{
  TEUCHOS_ASSERT(handlesRequest(rm)); // design by contract

  const std::string& field(getHandledField(*rm.getParameterList()));

  // Check if the field is in the main UGI.  If it's not, assume it's in the
  // auxiliary UGI.
  bool useAux(true);
  std::vector<Teuchos::RCP<panzer::UniqueGlobalIndexer<int, GlobalOrdinalT>>>
    fieldDOFMngrs = blocked_ugi_->getFieldDOFManagers();
  for (int b(0); b < static_cast<int>(fieldDOFMngrs.size()); ++b)
  {
    for (int f(0); f < fieldDOFMngrs[b]->getNumFields(); ++f)
    {
      if (fieldDOFMngrs[b]->getFieldString(f) == field)
        useAux = false;
    }
  }

  int block(-1);
  if (useAux)
    block =
      aux_blocked_ugi_->getFieldBlock(aux_blocked_ugi_->getFieldNum(field));
  else
    block = blocked_ugi_->getFieldBlock(blocked_ugi_->getFieldNum(field));

  // Empty...  Nothing to do.
  if (rm.getParameterList()->isType<std::string>("Coordinates-Epetra")) {
    buildArrayToVectorEpetra(block, field, useAux);
    buildCoordinatesEpetra(field, useAux);
  } else {
    buildArrayToVectorTpetra(block, field, useAux);
    buildCoordinatesTpetra(field, useAux);
  }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::setFieldByKey(const std::string & key,const std::string & field,Teuchos::ParameterList & pl) const
{
   // x-, y-, z-coordinates needed for ML, Coordinates needed for MueLu
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
   } else if(key == "Coordinates") {
     pl.set<Teuchos::RCP<Tpetra::MultiVector<double,LocalOrdinalT,GlobalOrdinalT,Node> > >(key,coordsVecTp_);
   } else if(key == "Coordinates-Epetra") {
      pl.set<Teuchos::RCP<Epetra_MultiVector> >("Coordinates",coordsVecEp_);
      // pl.remove("Coordinates-Epetra");
   }
   else
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::runtime_error,
                         "ParameterListCallback cannot handle key=\"" << key << "\"");
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::buildArrayToVectorTpetra(int block,const std::string & field, const bool useAux)
{
   if(arrayToVectorTpetra_[field]==Teuchos::null) {
      Teuchos::RCP<const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > ugi;
      if(useAux)
        ugi = aux_blocked_ugi_->getFieldDOFManagers()[block];
      else
        ugi = blocked_ugi_->getFieldDOFManagers()[block];
      arrayToVectorTpetra_[field] = Teuchos::rcp(new panzer::ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node>(ugi));
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::buildArrayToVectorEpetra(int block,const std::string & field, const bool useAux)
{
   if(arrayToVectorEpetra_[field]==Teuchos::null) {
      Teuchos::RCP<const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > ugi;
      if(useAux)
        ugi = aux_blocked_ugi_->getFieldDOFManagers()[block];
      else
        ugi = blocked_ugi_->getFieldDOFManagers()[block];
      arrayToVectorEpetra_[field] = Teuchos::rcp(new panzer::ArrayToFieldVectorEpetra<LocalOrdinalT,GlobalOrdinalT,Node>(ugi));
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::buildCoordinatesTpetra(const std::string & field, const bool useAux)
{
   std::map<std::string,Kokkos::DynRankView<double,PHX::Device> > data;

   Teuchos::RCP<const panzer::Intrepid2FieldPattern> fieldPattern = getFieldPattern(field,useAux);

   std::vector<std::string> elementBlocks;
   if(useAux)
     aux_blocked_ugi_->getElementBlockIds(elementBlocks);
   else
     blocked_ugi_->getElementBlockIds(elementBlocks);
   for(std::size_t i=0;i<elementBlocks.size();++i) {
      std::string blockId = elementBlocks[i];
      std::vector<std::size_t> localCellIds;

      // allocate block of data to store coordinates
      Kokkos::DynRankView<double,PHX::Device> & fieldData = data[blockId];
      fieldData = Kokkos::DynRankView<double,PHX::Device>("fieldData",connManager_->getElementBlock(blockId).size(),fieldPattern->numberIds());

      if(fieldPattern->supportsInterpolatoryCoordinates()) {
         // get degree of freedom coordiantes
         connManager_->getDofCoords(blockId,*fieldPattern,localCellIds,fieldData);
      }
      else {
         Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
         out.setOutputToRootOnly(-1);
         out << "WARNING: In ParameterListCallback::buildCoordinates(), the Intrepid2::FieldPattern in "
             << "block \"" << blockId << "\" does not support interpolatory coordinates. "
             << "This may be fine if coordinates are not actually needed. However if they are then bad things "
             << "will happen. Enjoy!" << std::endl;

         return;
      }
   }

   coordsVecTp_ = arrayToVectorTpetra_[field]->template getDataVector<double>(field,data);

   switch(coordsVecTp_->getNumVectors()) {
   case 3:
      zcoords_[field].resize(coordsVecTp_->getLocalLength());
      coordsVecTp_->getVector(2)->get1dCopy(Teuchos::arrayViewFromVector(zcoords_[field]));
      // Intentional fall-through.
   case 2:
      ycoords_[field].resize(coordsVecTp_->getLocalLength());
      coordsVecTp_->getVector(1)->get1dCopy(Teuchos::arrayViewFromVector(ycoords_[field]));
      // Intentional fall-through.
   case 1:
      xcoords_[field].resize(coordsVecTp_->getLocalLength());
      coordsVecTp_->getVector(0)->get1dCopy(Teuchos::arrayViewFromVector(xcoords_[field]));
      break;
   default:
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
             "ParameterListCallback::buildCoordinates: Constructed multivector has nonphysical dimensions.");
      break;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::buildCoordinatesEpetra(const std::string & field, const bool useAux)
{
   std::map<std::string,Kokkos::DynRankView<double,PHX::Device> > data;

   Teuchos::RCP<const panzer::Intrepid2FieldPattern> fieldPattern = getFieldPattern(field,useAux);

   std::vector<std::string> elementBlocks;
   if(useAux)
     aux_blocked_ugi_->getElementBlockIds(elementBlocks);
   else
     blocked_ugi_->getElementBlockIds(elementBlocks);
   for(std::size_t i=0;i<elementBlocks.size();++i) {
      std::string blockId = elementBlocks[i];
      std::vector<std::size_t> localCellIds;

      // allocate block of data to store coordinates
      Kokkos::DynRankView<double,PHX::Device> & fieldData = data[blockId];
      fieldData = Kokkos::DynRankView<double,PHX::Device>("fieldData",connManager_->getElementBlock(blockId).size(),fieldPattern->numberIds());

      if(fieldPattern->supportsInterpolatoryCoordinates()) {
         // get degree of freedom coordiantes
         connManager_->getDofCoords(blockId,*fieldPattern,localCellIds,fieldData);
      }
      else {
         Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
         out.setOutputToRootOnly(-1);
         out << "WARNING: In ParameterListCallback::buildCoordinates(), the Intrepid2::FieldPattern in "
             << "block \"" << blockId << "\" does not support interpolatory coordinates. "
             << "This may be fine if coordinates are not actually needed. However if they are then bad things "
             << "will happen. Enjoy!" << std::endl;

         return;
      }
   }

   coordsVecEp_ = arrayToVectorEpetra_[field]->template getDataVector<double>(field,data);
}


template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
std::string ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>::
getHandledField(const Teuchos::ParameterList & pl) const
{
  // because this method assumes handlesRequest is true, this call will succeed
  if(pl.isType<std::string>("x-coordinates"))
    return pl.get<std::string>("x-coordinates");
  else if(pl.isType<std::string>("Coordinates"))
    return pl.get<std::string>("Coordinates");
  else if(pl.isType<std::string>("Coordinates-Epetra"))
    return pl.get<std::string>("Coordinates-Epetra");
  else
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,"Neither x-coordinates nor Coordinates or Coordinates-Epetra field provided.");
    
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
Teuchos::RCP<const panzer::Intrepid2FieldPattern> ParameterListCallbackBlocked<LocalOrdinalT,GlobalOrdinalT,Node>
::getFieldPattern(const std::string & fieldName, const bool useAux) const
{
  std::vector<std::string> elementBlocks;
  if(useAux)
    aux_blocked_ugi_->getElementBlockIds(elementBlocks);
  else
    blocked_ugi_->getElementBlockIds(elementBlocks);

  for(std::size_t e=0;e<elementBlocks.size();e++) {
    std::string blockId = elementBlocks[e];

    if(blocked_ugi_->fieldInBlock(fieldName,blockId))
      return Teuchos::rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(blocked_ugi_->getFieldPattern(blockId,fieldName),true);

    if(aux_blocked_ugi_->fieldInBlock(fieldName,blockId))
      return Teuchos::rcp_dynamic_cast<const panzer::Intrepid2FieldPattern>(aux_blocked_ugi_->getFieldPattern(blockId,fieldName),true);
  }

  return Teuchos::null;
}


} 

#endif
