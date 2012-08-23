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

#ifndef __Panzer_STK_ParameterListCallback_hpp__
#define __Panzer_STK_ParameterListCallback_hpp__

#ifdef HAVE_TEKO

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Teko_RequestCallback.hpp"

#include "Panzer_STKConnManager.hpp"
#include "Panzer_UniqueGlobalIndexer_Utilities.hpp"

#include <vector>

namespace panzer_stk {

class STKConnManager;

/** Implements an interface used by the Teko request handler mechanism.
  * This particular class is usesd most frequently with an ML preconditioner that
  * requres the nodal coordinates for repartitioning.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node=Kokkos::DefaultNode::DefaultNodeType>
class ParameterListCallback : public Teko::RequestCallback<Teuchos::RCP<Teuchos::ParameterList> > {
public:
  ParameterListCallback(const std::string & coordFieldName,
                        const map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > & fp,
                        const Teuchos::RCP<const panzer_stk::STKConnManager> & connManager, 
                        const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi);

   Teuchos::RCP<Teuchos::ParameterList> request(const Teko::RequestMesg & rm);

   bool handlesRequest(const Teko::RequestMesg & rm);

   void preRequest(const Teko::RequestMesg & rm);

   const std::vector<double> & getXCoordsVector() const { return xcoords_; }
   const std::vector<double> & getYCoordsVector() const { return ycoords_; }
   const std::vector<double> & getZCoordsVector() const { return zcoords_; }

   //! Return the internally stored arraytofieldvector object. May return null if none constructed
   Teuchos::RCP<const panzer::ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node> > getArrayToFieldVector() const
   { return arrayToVector_; }

   void buildCoordinates();
   void buildArrayToVector();

private:

   void setFieldByKey(const std::string & key,Teuchos::ParameterList & pl) const;

   std::string coordFieldName_;
   std::map<std::string,Teuchos::RCP<const panzer::IntrepidFieldPattern> > fieldPatterns_;
   Teuchos::RCP<const panzer_stk::STKConnManager> connManager_;
   Teuchos::RCP<const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > ugi_;
   bool coordinatesBuilt_;
 
   std::vector<double> xcoords_;
   std::vector<double> ycoords_;
   std::vector<double> zcoords_;

   mutable Teuchos::RCP<const panzer::ArrayToFieldVector<LocalOrdinalT,GlobalOrdinalT,Node> > arrayToVector_;
};

}

#include "Panzer_STK_ParameterListCallback_impl.hpp"

#endif // HAVE_TEKO

#endif
