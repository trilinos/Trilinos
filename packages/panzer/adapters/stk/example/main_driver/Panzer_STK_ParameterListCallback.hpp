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

private:
   void buildCoordinates();
   void buildArrayToVector();

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

#include "Panzer_STK_ParameterListCallbackT.hpp"

#endif // HAVE_TEKO

#endif
