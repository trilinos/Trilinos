#ifndef __Panzer_DOFLibrary_hpp__
#define __Panzer_DOFLibrary_hpp__

#include "Panzer_config.hpp"

#include "Teuchos_RCP.hpp"

#include "Panzer_Basis.hpp"
#include "Panzer_PureBasis.hpp"

namespace panzer {

/** There is one of these objects per equation set.
  */
class FieldLayoutLibrary {
public:
   /** Add a field associated with a basis to the library.
     */
   void addFieldAndLayout(const std::string & fieldName,
                         const Teuchos::RCP<panzer::Basis> & basis);   

   //! Get the basis associated with a particular field.
   Teuchos::RCP<panzer::Basis> lookup(const std::string & fieldName) const;

private:

   //! Basic mapped storage.
   std::map<std::string,Teuchos::RCP<panzer::Basis> > fieldToLayout_;
};

/** Build a container that holds, and provides
  * easy lookups for each fields basis. This provides
  * further functionality by providing a class that
  * oversees the marriage of the basis and integration
  * rule objects. There is one of these objects per
  * physics block.
  */
class FieldLibrary {
public:

   //! Get the basis associated with a particular field.
   Teuchos::RCP<panzer::PureBasis> lookup(const std::string & fieldName) const;

   /** Add a field associated witha basis to the library.
     */
   void addFieldAndBasis(const std::string & fieldName,
                 const Teuchos::RCP<panzer::PureBasis> & basis);   

   /** Given an integration rule build a BasisIRLibrary which
     * oversees the marriage of the integration rule and the basis
     * into a BasisIRLayout.
     */
   Teuchos::RCP<const FieldLayoutLibrary> buildFieldLayoutLibrary(panzer::IntegrationRule & ir) const;

private:

   //! Basic mapped storage.
   std::map<std::string,Teuchos::RCP<panzer::PureBasis> > fieldToBasis_;
};

}

#endif
