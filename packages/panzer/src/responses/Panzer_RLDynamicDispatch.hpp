#ifndef __Panzer_RLDynamicDispatch_hpp__
#define __Panzer_RLDynamicDispatch_hpp__

#include <string>

#include "Teuchos_Ptr.hpp"
#include "Sacado_mpl_for_each.hpp"
#include "Phalanx_TemplateManager.hpp"

namespace panzer {

// forward declaration
template <typename Traits> class ResponseLibrary;

/** \file Panzer_RLDynamicDispatch.hpp 
  * \brief Of use only by the internals of the response library.
  *
  * This provides a mechnism where a templated member function
  * can be called a run time using a dynamic call. This is specialized
  * to the ResponseLibrary, is meant to handle those functions that
  * use templating on the evaluation type.
  */

/** Abstract base class which acts as a surrogate for the
  * dynamic/static call.
  */
template <typename Traits>
class RLDynamicDispatchBase {
public:
   virtual ~RLDynamicDispatchBase() {}

   virtual void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType) = 0;
};

/** A concrete instantiation of the response library
  * dynamic dispatch mechanism. This is the step that
  * calls the compile time polymorphic type from the dynamic string.
  * (CT stands for Compile Time)
  */
template <typename EvalT,typename Traits>
class RLDynamicDispatchCT : public RLDynamicDispatchBase<Traits> {
public:
   RLDynamicDispatchCT(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl)
      : responseLibrary_(rl) {}

   virtual ~RLDynamicDispatchCT() {}

   virtual void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType)
   { 
      // sanity check
      TEST_FOR_EXCEPTION(PHX::TypeString<EvalT>::value!=evalType,std::logic_error,
                         "panzer::RLDynamicDispatch: Somethihng is seriously wrong, dynamic evaluation type \""+evalType+"\" "
                         "is not equal to static type \""+PHX::TypeString<EvalT>::value+"\"!"); 
      responseLibrary_->template reserveVolumeResponse<EvalT>(rid,eBlock); 
   }

private:
   RLDynamicDispatchCT();
   RLDynamicDispatchCT(const RLDynamicDispatchCT<EvalT,Traits> &);

   Teuchos::Ptr<ResponseLibrary<Traits> > responseLibrary_;
};

template <typename Traits>
class RLDynamicDispatch : public RLDynamicDispatchBase<Traits> {
public:
   RLDynamicDispatch();

   void buildObjects(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl);

   virtual void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock,const std::string & evalType)
   { dynToStaticMap_[evalType]->reserveVolumeResponse(rid,eBlock,evalType); }

   const std::vector<std::string> & getEvaluationNames() const
   { return evaluationNames_; }

private:
   /** Class to turn an mpl::vector into a string vector
     * using the <code>PHX::TypeString</code> class.
     */
   struct NameFill {
      //! mutable so that we can modify this and use the for_each in Sacado
      mutable std::vector<std::string> & vec;

      NameFill(std::vector<std::string> & v) : vec(v) {}
      NameFill(const NameFill & n) : vec(n.vec) {}

      template <typename T> void operator()(T x) const
      { vec.push_back(PHX::TypeString<T>::value); }
      
   private:
      NameFill();
   };

   //! For use with a template manager for building dynamic dispatch CT objects
   struct RLDynDispCTBuilder {
      Teuchos::Ptr<ResponseLibrary<Traits> > responseLibrary_;

      RLDynDispCTBuilder(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl)
         : responseLibrary_(rl) {}
 
      template <typename EvalT>
      Teuchos::RCP<RLDynamicDispatchBase<Traits> > build() const
      { return Teuchos::rcp(new RLDynamicDispatchCT<EvalT,Traits>(responseLibrary_)); }
   };

   std::vector<std::string> evaluationNames_;

   std::map<std::string,Teuchos::RCP<RLDynamicDispatchBase<Traits> > > dynToStaticMap_;
};

template <typename Traits>
RLDynamicDispatch<Traits>::RLDynamicDispatch()
{
   typedef typename Traits::EvalTypes EvalTypes;
   Sacado::mpl::for_each<EvalTypes>(NameFill(evaluationNames_));
}

template <typename Traits>
void RLDynamicDispatch<Traits>::buildObjects(const Teuchos::Ptr<ResponseLibrary<Traits> > & rl)
{
   typedef PHX::TemplateManager<typename Traits::EvalTypes,RLDynamicDispatchBase<Traits>,RLDynamicDispatchCT<_,Traits> > 
           DynDispatchTM;

   DynDispatchTM dynDispatchTM;
   dynDispatchTM.buildObjects(RLDynDispCTBuilder(rl));

   // build map from string to dynamic dispatch object
   std::vector<std::string>::const_iterator nameItr = evaluationNames_.begin();
   typename DynDispatchTM::iterator dynDspItr(dynDispatchTM.begin());
   for(; nameItr!=evaluationNames_.end() && dynDspItr!=dynDispatchTM.end();
       ++nameItr,++dynDspItr) {
      dynToStaticMap_[*nameItr] = dynDspItr.rcp();
   }
}

}

#endif
