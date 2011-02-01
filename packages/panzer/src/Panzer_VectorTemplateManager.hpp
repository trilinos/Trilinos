#ifndef __Panzer_VectorTemplateManager_hpp__
#define __Panzer_VectorTemplateManager_hpp__

#include <vector>

#include "Teuchos_RCP.hpp"

#include "Phalanx_TemplateManager.hpp"

#include "boost/mpl/apply.hpp"

namespace panzer {

/** A vector template manager that handles the construction
  * of objects templated on a particular set of types.
  */
template <typename SeqTypes,typename BaseT,typename ObjectT>
class VectorTemplateManager {
public:
   VectorTemplateManager() {}
   ~VectorTemplateManager() {}

   /** Add the contents of a builder to this vector object
     */
   template <typename BuilderT>
   void buildAndPushBackObjects(const BuilderT & builder);

   /** Get a vector to the base object pointers
     */
   template <typename EvalT>
   void getAsBase(std::vector<Teuchos::RCP<BaseT> > & baseObjects);

   /** Get a vector to the base object pointers
     */
   template <typename EvalT>
   void getAsBase(std::vector<Teuchos::RCP<const BaseT> > & baseObjects) const;

   /** Get a vector to the object pointers
     */
   template <typename EvalT>
   void getAsObject(std::vector<Teuchos::RCP<typename boost::mpl::apply<ObjectT,EvalT>::type > > & objects);

   /** Get a vector to the object pointers
     */
   template <typename EvalT>
   void getAsObject(std::vector<Teuchos::RCP<const typename boost::mpl::apply<ObjectT,EvalT>::type > > & objects) const;

private:
   typedef PHX::TemplateManager<SeqTypes,BaseT,ObjectT> ManagerType;

   std::vector<Teuchos::RCP<ManagerType> > managers_;
};

#include "Panzer_VectorTemplateManagerT.hpp"

}

#endif
