// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_ParameterListCallback_hpp__
#define __Panzer_STK_ParameterListCallback_hpp__

#include "PanzerAdaptersSTK_config.hpp"
#ifdef PANZER_HAVE_TEKO

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Teko_RequestCallback.hpp"

#include "Panzer_STKConnManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"

#include <vector>
#include <map>

namespace panzer_stk {

class STKConnManager;

/** Implements an interface used by the Teko request handler mechanism.
  * This particular class is usesd most frequently with an ML preconditioner that
  * requres the nodal coordinates for repartitioning.
  */
class ParameterListCallback : public Teko::RequestCallback<Teuchos::RCP<Teuchos::ParameterList> > {
public:
  ParameterListCallback(const std::string & coordFieldName,
                        const std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > & fp,
                        const Teuchos::RCP<const panzer_stk::STKConnManager> & connManager,
                        const Teuchos::RCP<const panzer::GlobalIndexer> & ugi);

   Teuchos::RCP<Teuchos::ParameterList> request(const Teko::RequestMesg & rm);

   bool handlesRequest(const Teko::RequestMesg & rm);

   void preRequest(const Teko::RequestMesg & rm);

   const std::vector<double> & getCoordsVector(unsigned dim) const
   { switch(dim) {
     case 0:
       return getXCoordsVector();
     case 1:
       return getYCoordsVector();
     case 2:
       return getZCoordsVector();
     default:
       TEUCHOS_ASSERT(false);
     }
     TEUCHOS_ASSERT(false);
     return xcoords_; // should never get here!
   }
   const std::vector<double> & getXCoordsVector() const { return xcoords_; }
   const std::vector<double> & getYCoordsVector() const { return ycoords_; }
   const std::vector<double> & getZCoordsVector() const { return zcoords_; }

   //! Return the internally stored arraytofieldvector object. May return null if none constructed
   Teuchos::RCP<const panzer::ArrayToFieldVector> getArrayToFieldVector() const
   { return arrayToVector_; }

   void buildCoordinates();
   void buildArrayToVector();

   //! Store a vector solely for the purpose of making it persist with this object
   void storeExtraVector(const Teuchos::RCP<const std::vector<double> > & extra)
   { extraVecs_.push_back(extra); }

private:

   void setFieldByKey(const std::string & key,Teuchos::ParameterList & pl) const;

   std::string coordFieldName_;
   std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > fieldPatterns_;
   Teuchos::RCP<const panzer_stk::STKConnManager> connManager_;
   Teuchos::RCP<const panzer::GlobalIndexer> ugi_;
   bool coordinatesBuilt_;

   std::vector<double> xcoords_;
   std::vector<double> ycoords_;
   std::vector<double> zcoords_;

   mutable Teuchos::RCP<const panzer::ArrayToFieldVector> arrayToVector_;
   std::vector<Teuchos::RCP<const std::vector<double> > > extraVecs_;
};

}

#endif // PANZER_HAVE_TEKO

#endif
