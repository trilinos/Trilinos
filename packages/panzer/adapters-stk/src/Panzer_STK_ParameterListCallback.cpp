// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "PanzerAdaptersSTK_config.hpp"
#ifdef PANZER_HAVE_TEKO

#include "Panzer_STK_ParameterListCallback.hpp"

namespace panzer_stk {

using Teuchos::RCP;
using Teuchos::rcp;

ParameterListCallback::ParameterListCallback(const std::string & coordFieldName,
                                             const std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > & fps,
                                             const Teuchos::RCP<const panzer_stk::STKConnManager> & connManager,
                                             const Teuchos::RCP<const panzer::GlobalIndexer> & ugi)
   : coordFieldName_(coordFieldName), fieldPatterns_(fps), connManager_(connManager), ugi_(ugi), coordinatesBuilt_(false)
{ }

Teuchos::RCP<Teuchos::ParameterList> ParameterListCallback::request(const Teko::RequestMesg & rm)
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

bool ParameterListCallback::handlesRequest(const Teko::RequestMesg & rm)
{
   // check if is a parameter list message, and that the parameter
   // list contains the right fields
   if(rm.getName()=="Parameter List") return true;
   else return false;
}

void ParameterListCallback::preRequest(const Teko::RequestMesg & rm)
{
   TEUCHOS_ASSERT(handlesRequest(rm)); // design by contract

   // empty...nothing to do
   buildArrayToVector();
   buildCoordinates();
}

void ParameterListCallback::setFieldByKey(const std::string & key,Teuchos::ParameterList & pl) const
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

void ParameterListCallback::buildArrayToVector()
{
   if(arrayToVector_==Teuchos::null)
      arrayToVector_ = Teuchos::rcp(new panzer::ArrayToFieldVector(ugi_));
}

void ParameterListCallback::buildCoordinates()
{
   TEUCHOS_ASSERT(fieldPatterns_.size()>0); // must be at least one field pattern

   std::map<std::string,Kokkos::DynRankView<double,PHX::Device> > data;

   std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> >::const_iterator itr;
   for(itr=fieldPatterns_.begin();itr!=fieldPatterns_.end();++itr) {
      std::string blockId = itr->first;
      Teuchos::RCP<const panzer::Intrepid2FieldPattern> fieldPattern = itr->second;
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

         coordinatesBuilt_ = true;
         return;
      }
   }

   Teuchos::RCP<Tpetra::MultiVector<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > resultVec
      = arrayToVector_->template getDataVector<double>(coordFieldName_,data);

   switch(resultVec->getNumVectors()) {
   case 3:
      zcoords_.resize(resultVec->getLocalLength());
      resultVec->getVector(2)->get1dCopy(Teuchos::arrayViewFromVector(zcoords_));
      // Intentional fall-through.
   case 2:
      ycoords_.resize(resultVec->getLocalLength());
      resultVec->getVector(1)->get1dCopy(Teuchos::arrayViewFromVector(ycoords_));
      // Intentional fall-through.
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
