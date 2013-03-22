#ifndef __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_impl_hpp__
#define __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_impl_hpp__

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ScatterFields.hpp"
#include "Panzer_STK_ScatterVectorFields.hpp"
#include "Panzer_PointValues_Evaluator.hpp"
#include "Panzer_BasisValues_Evaluator.hpp"
#include "Panzer_DOF.hpp"

#include <boost/unordered_set.hpp>

namespace panzer_stk {

namespace {
   //! A dummy response for local use, is only used by the response library
   class Response_STKDummy : public panzer::ResponseBase {
   public:
     Response_STKDummy(const std::string & rn)
       : ResponseBase(rn) {}
     virtual void scatterResponse() {}
     virtual void initializeResponse() {}
   private:
     Response_STKDummy();
     Response_STKDummy(const Response_STKDummy &);
   };
}

template <typename EvalT> 
Teuchos::RCP<panzer::ResponseBase> ResponseEvaluatorFactory_SolutionWriter<EvalT>:: 
buildResponseObject(const std::string & responseName) const
{
   return Teuchos::rcp(new Response_STKDummy(responseName));
}
   
template <typename EvalT> 
void ResponseEvaluatorFactory_SolutionWriter<EvalT>:: 
buildAndRegisterEvaluators(const std::string & responseName,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList & user_data) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  // this will help so we can print out any unused scaled fields as a warning
  boost::unordered_set<std::string> scaledFieldsHash = scaledFieldsHash_;

  const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases = physicsBlock.getBases();
  std::map<std::string,std::vector<std::string> > basisBucket;
  bucketByBasisType(physicsBlock.getProvidedDOFs(),basisBucket);

  // add this for HCURL and HDIV basis, only want to add them once: evaluate vector fields at centroid
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  RCP<panzer::PointRule> centroidRule;
  for(std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator itr=bases.begin();
      itr!=bases.end();++itr) {

    if(itr->second->isVectorBasis()) {
      centroidRule = rcp(new panzer::PointRule("Centroid",1,physicsBlock.cellData()));

      // compute centroid
      Intrepid::FieldContainer<double> centroid;
      computeReferenceCentroid(bases,physicsBlock.cellData().baseCellDimension(),centroid);

      // build pointe values evaluator
      RCP<PHX::Evaluator<panzer::Traits> > evaluator  = 
         rcp(new panzer::PointValues_Evaluator<EvalT,panzer::Traits>(centroidRule,centroid));
      fm.template registerEvaluator<EvalT>(evaluator);

      break; // get out of the loop, only need one evaluator
    }
  }

  // add evaluators for each field
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(std::map<std::string,std::vector<std::string> >::const_iterator itr=basisBucket.begin();
      itr!=basisBucket.end();++itr) {

    std::string basisName = itr->first;
    const std::vector<std::string> & fields = itr->second;

    std::map<std::string,Teuchos::RCP<panzer::PureBasis> >::const_iterator found = bases.find(basisName);
    TEUCHOS_ASSERT(found!=bases.end());
    Teuchos::RCP<const panzer::PureBasis> basis = found->second;
    
    // write out nodal fields
    if(basis->getElementSpace()==panzer::PureBasis::HGRAD) {
      
      // determine if user has modified field scalar for each field to be written to STK
      std::vector<double> scalars(fields.size(),1.0); // fill with 1.0 
      for(std::size_t f=0;f<fields.size();f++) { 
        boost::unordered_map<std::string,double>::const_iterator f2s_itr = fieldToScalar_.find(fields[f]);

        // if scalar is found, include it in the vector and remove the field from the
        // hash table so it won't be included in the warning message.
        if(f2s_itr!=fieldToScalar_.end()) { 
          scalars[f] = f2s_itr->second;
          scaledFieldsHash.erase(fields[f]);
        }
      }

      Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval = 
        Teuchos::rcp(new ScatterFields<EvalT,panzer::Traits>("STK HGRAD Scatter Basis " +basis->name(),
                                                      mesh_, basis, fields,scalars));

      // register and require evaluator fields
      fm.template registerEvaluator<EvalT>(eval);
      fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
    }
    else if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
      TEUCHOS_ASSERT(centroidRule!=Teuchos::null);

      // register basis values evaluator
      {
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer::BasisValues_Evaluator<EvalT,panzer::Traits>(centroidRule,basis));
        fm.template registerEvaluator<EvalT>(evaluator);
      }

      // add a DOF_PointValues for each field
      std::vector<std::string> pointFields;
      for(std::size_t f=0;f<fields.size();f++) {
        Teuchos::ParameterList p;
        p.set("Name",fields[f]);
        p.set("Basis",basis);
        p.set("Point Rule",centroidRule.getConst());
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer::DOF_PointValues<EvalT,panzer::Traits>(p));

        fm.template registerEvaluator<EvalT>(evaluator);

        pointFields.push_back(fields[f]+"_"+centroidRule->getName());
      }

      // add the scatter field evaluator for this basis
      {
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer_stk::ScatterVectorFields<EvalT,panzer::Traits>("STK HCURL Scatter Basis " +basis->name(),
                                                                              mesh_,centroidRule,fields));

        fm.template registerEvaluator<EvalT>(evaluator);
        fm.template requireField<EvalT>(*evaluator->evaluatedFields()[0]); // require the dummy evaluator
      }
    }
  }

  // print warning message for any unused scaled fields
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);

  for(boost::unordered_set<std::string>::const_iterator itr=scaledFieldsHash.begin();
      itr!=scaledFieldsHash.end();itr++) { 
    out << "WARNING: STK Solution Writer did not scale the field \"" << *itr << "\" "
        << "because it was not written." << std::endl;
  }
}

template <typename EvalT>
void ResponseEvaluatorFactory_SolutionWriter<EvalT>::
bucketByBasisType(const std::vector<panzer::StrPureBasisPair> & providedDofs,
                  std::map<std::string,std::vector<std::string> > & basisBucket)
{
   // this should be self explanatory
   for(std::size_t i=0;i<providedDofs.size();i++) {
      std::string fieldName = providedDofs[i].first;
      Teuchos::RCP<const panzer::PureBasis> basis = providedDofs[i].second;

      basisBucket[basis->name()].push_back(fieldName);
   }
}

template <typename EvalT>
void ResponseEvaluatorFactory_SolutionWriter<EvalT>::
computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases,
                         int baseDimension,
                         Intrepid::FieldContainer<double> & centroid) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   centroid.resize(1,baseDimension);

   // loop over each possible basis
   for(std::map<std::string,RCP<panzer::PureBasis> >::const_iterator itr=bases.begin();
       itr!=bases.end();++itr) {

      // see if this basis has coordinates
      RCP<Intrepid::Basis<double,Intrepid::FieldContainer<double> > > intrepidBasis = itr->second->getIntrepidBasis();
      RCP<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<double> > > basisCoords 
         = rcp_dynamic_cast<Intrepid::DofCoordsInterface<Intrepid::FieldContainer<double> > >(intrepidBasis);

      if(basisCoords==Teuchos::null) // no coordinates...move on
         continue;

      // we've got coordinates, lets commpute the "centroid"
      Intrepid::FieldContainer<double> coords(intrepidBasis->getCardinality(),
                                              intrepidBasis->getBaseCellTopology().getDimension());
      basisCoords->getDofCoords(coords);
      TEUCHOS_ASSERT(coords.rank()==2);
      TEUCHOS_ASSERT(coords.dimension(1)==baseDimension);

      for(int i=0;i<coords.dimension(0);i++)
         for(int d=0;d<coords.dimension(1);d++)
            centroid(0,d) += coords(i,d);

      // take the average
      for(int d=0;d<coords.dimension(1);d++)
         centroid(0,d) /= coords.dimension(0);

      return;
   }

   // no centroid was found...die
   TEUCHOS_ASSERT(false);
}

template <typename EvalT>
void ResponseEvaluatorFactory_SolutionWriter<EvalT>::
scaleField(const std::string & fieldName,double fieldScalar)
{
  fieldToScalar_[fieldName] = fieldScalar;
}

}

#endif
