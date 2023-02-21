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

#ifndef __Panzer_STK_ParameterListCallbackBlocked_hpp__
#define __Panzer_STK_ParameterListCallbackBlocked_hpp__

#include "PanzerAdaptersSTK_config.hpp"
#ifdef PANZER_HAVE_TEKO

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Teko_RequestCallback.hpp"

#include "Panzer_STKConnManager.hpp"
#include "Panzer_GlobalIndexer_Utilities.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_GlobalIndexer_EpetraUtilities.hpp"
#endif
#include "Panzer_BlockedDOFManager.hpp"

#include <vector>
#include <map>

namespace panzer_stk {

class STKConnManager;

/** Implements an interface used by the Teko request handler mechanism.
  * This particular class is usesd most frequently with an ML preconditioner that
  * requres the nodal coordinates for repartitioning.
  */
class ParameterListCallbackBlocked : public Teko::RequestCallback<Teuchos::RCP<Teuchos::ParameterList> > {
public:
  ParameterListCallbackBlocked(const Teuchos::RCP<const panzer_stk::STKConnManager> & connManager,
                               const Teuchos::RCP<const panzer::BlockedDOFManager> & blkDofs,
                               const Teuchos::RCP<const panzer::BlockedDOFManager> & auxBlkDofs=Teuchos::null);

  Teuchos::RCP<Teuchos::ParameterList> request(const Teko::RequestMesg & rm);

  bool handlesRequest(const Teko::RequestMesg & rm);

  void preRequest(const Teko::RequestMesg & rm);

private:

  bool isField(const std::string& field) const
  {
    // Check both the main and auxiliary UGIs.
    bool useAux(true);
    std::vector<Teuchos::RCP<panzer::GlobalIndexer>>
      fieldDOFMngrs = blocked_ugi_->getFieldDOFManagers();
    for (int b(0); b < static_cast<int>(fieldDOFMngrs.size()); ++b)
    {
      for (int f(0); f < fieldDOFMngrs[b]->getNumFields(); ++f)
      {
        if (fieldDOFMngrs[b]->getFieldString(f) == field)
          useAux = false;
      }
    }
    if (useAux)
      return (aux_blocked_ugi_->getFieldNum(field) != -1);
    else
      return (blocked_ugi_->getFieldNum(field) != -1);
  }

  void buildArrayToVectorTpetra(int block,const std::string & field, const bool useAux = false);
  void buildCoordinatesTpetra(const std::string & field, const bool useAux = false);

#ifdef PANZER_HAVE_EPETRA_STACK
  void buildArrayToVectorEpetra(int block,const std::string & field, const bool useAux = false);
  void buildCoordinatesEpetra(const std::string & field, const bool useAux = false);
#endif

  // this method assumes handlesRequest(rm)==true
  std::string getHandledField(const Teuchos::ParameterList & pl) const;

  void setFieldByKey(const std::string & key,const std::string & field,Teuchos::ParameterList & pl) const;

  // Get the coordinate vector by field, note that it will check to make sure "buildCoordinates" has
  // been called.
  const std::vector<double> & getCoordinateByField(int dim,const std::string & field) const;

  // Access the field pattern associated with this field (this will not work if the field pattern
  // and element blocks are different
  Teuchos::RCP<const panzer::Intrepid2FieldPattern> getFieldPattern(const std::string & fieldName, const bool useAux = false) const;

  // Generally used members
  Teuchos::RCP<const panzer_stk::STKConnManager> connManager_;
  Teuchos::RCP<const panzer::BlockedDOFManager> blocked_ugi_;
  Teuchos::RCP<const panzer::BlockedDOFManager> aux_blocked_ugi_;

  std::map<std::string,Teuchos::RCP<const panzer::Intrepid2FieldPattern> > fieldPatterns_;

  // look up by field name (field name to coordinates

  std::map<std::string,std::vector<double> > xcoords_;
  std::map<std::string,std::vector<double> > ycoords_;
  std::map<std::string,std::vector<double> > zcoords_;

  mutable std::map<std::string,Teuchos::RCP<const panzer::ArrayToFieldVector> > arrayToVectorTpetra_;
#ifdef PANZER_HAVE_EPETRA_STACK
  mutable std::map<std::string,Teuchos::RCP<const panzer::ArrayToFieldVectorEpetra> > arrayToVectorEpetra_;
#endif

  Teuchos::RCP<Tpetra::MultiVector<double,int,panzer::GlobalOrdinal,panzer::TpetraNodeType> > coordsVecTp_;
#ifdef PANZER_HAVE_EPETRA_STACK
  Teuchos::RCP<Epetra_MultiVector> coordsVecEp_;
#endif

  bool returnTpetraObjects_;
};

}

#endif // PANZER_HAVE_TEKO

#endif
