#ifndef __Panzer_ResponseAggregator_Manager_hpp__
#define __Panzer_ResponseAggregator_Manager_hpp__

#include "Panzer_ResponseAggregatorBase.hpp"

#include "Phalanx_TemplateManager.hpp"

namespace panzer {

template <typename TraitsT>
class ResponseAggregator_Manager {
public:
   typedef PHX::TemplateManager<typename TraitsT::EvalTypes,ResponseAggregatorBase<TraitsT>,ResponseAggregator<_,TraitsT> >
           AggregatorManager;

   /** Statically access an aggregator of a particular type
     *
     * \note This method throws no aggregator exists with the 
     *       type string and evaluation type template parameter.
     */
   template <typename EvalT>
   const ResponseAggregatorBase<TraitsT> & getAggregator(const std::string & type) const;

   /** Add a user defined aggregator to this container with a particular type.
     * The user must provide a builder as defined in the <code>PHX::TemplateManager</code>
     * documentation. This method will excercise the builder and construct a list of
     * aggregators for each evaluation type. In this case it is safe, and expected,
     * that the builder return Teuchos::null if no aggregator is appropriate for that evalution
     * type. This is not the case for the <code>PHX::TemplateManager</code>
     *
     * \note This method throws if this type has already been defined
     */
   template <typename Builder>
   void defineAggregatorTypeFromBuilder(const std::string & type,const Builder & builder);

   /** \defgroup isAggregator
     * Look ups for if an aggregator is contained in this manager.
     * @{
     */

   //! Is there an aggretator of a particular type.
   bool isAggregator(const std::string & type) const;

   //! Is there an aggregator that uses a particular evaluation type
   template <typename EvalT>
   bool isAggregator() const;

   //! Is there an aggregator that is a particular type and uses a specified evaluation type
   template <typename EvalT>
   bool isAggregator(const std::string & type) const;

   /** @} */

   /** Get aggregator manager associated with a type string
     *
     * \note This method throws if no aggregator manager is found.
     */
   const AggregatorManager & getAggregatorManager(const std::string & type) const;


   /** Add a set of panzer defined aggregators to the default aggregator manager.
     * In paricular the aggregator "Functional" is added as a Residual evaluator.
     */
   static void defineDefaultAggregators(ResponseAggregator_Manager<TraitsT> & aggMngr);
private:

   std::map<std::string,Teuchos::RCP<AggregatorManager> > aggregators_;
};

} // end namespace panzer

#include "Panzer_ResponseAggregator_ManagerT.hpp"

#endif
