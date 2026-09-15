// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_impl_hpp__
#define __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_impl_hpp__

#include "Panzer_STK_Interface.hpp"
#include "Panzer_STK_ScatterFields.hpp"
#include "Panzer_STK_ScatterVectorFields.hpp"
#include "Panzer_PointValues_Evaluator.hpp"
#include "Panzer_BasisValues_Evaluator.hpp"
#include "Panzer_DOF.hpp"
#include "Panzer_DOF_PointValues.hpp"

#include <unordered_set>

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
buildAndRegisterEvaluators(const std::string& /* responseName */,
                           PHX::FieldManager<panzer::Traits> & fm,
                           const panzer::PhysicsBlock & physicsBlock,
                           const Teuchos::ParameterList& /* user_data */) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  typedef std::pair<std::string,RCP<const panzer::PureBasis> > StrConstPureBasisPair;

  // this will help so we can print out any unused scaled fields as a warning
  std::unordered_set<std::string> scaledFieldsHash = scaledFieldsHash_;

  std::map<std::string,RCP<const panzer::PureBasis> > bases;
  std::map<std::string,std::vector<std::string> > basisBucket;
  { 
    const std::map<std::string,RCP<panzer::PureBasis> > & nc_bases = physicsBlock.getBases();
    bases.insert(nc_bases.begin(),nc_bases.end());
  }

  std::vector<StrConstPureBasisPair> allFields;

  // only add in solution fields if required

  if(!addCoordinateFields_ && addSolutionFields_) {
    // inject all the fields, including the coordinates (we will remove them shortly)
    allFields.insert(allFields.end(),physicsBlock.getProvidedDOFs().begin(),physicsBlock.getProvidedDOFs().end());
   

    // get a list of strings with fields to remove
    std::vector<std::string> removedFields;
    const std::vector<std::vector<std::string> > & coord_fields = physicsBlock.getCoordinateDOFs();
    for(std::size_t c=0;c<coord_fields.size();c++)
      for(std::size_t d=0;d<coord_fields[c].size();d++)
        removedFields.push_back(coord_fields[c][d]);

    // remove all coordinate fields
    deleteRemovedFields(removedFields,allFields); 
  }
  else if(addCoordinateFields_ && !addSolutionFields_) {
    Teuchos::RCP<const panzer::FieldLibraryBase> fieldLib = physicsBlock.getFieldLibraryBase();
    const std::vector<std::vector<std::string> > & coord_fields = physicsBlock.getCoordinateDOFs();
    
    // get the basis and field for each coordiante
    for(std::size_t c=0;c<coord_fields.size();c++) {
      for(std::size_t d=0;d<coord_fields[c].size();d++) {
        Teuchos::RCP<panzer::PureBasis> basis = // const_cast==yuck!
            Teuchos::rcp_const_cast<panzer::PureBasis>(fieldLib->lookupBasis(coord_fields[c][d]));

        // make sure they are inserted in the allFields list
        allFields.push_back(std::make_pair(coord_fields[c][d],basis));
      }
    }
  }
  else if(addSolutionFields_)
    allFields.insert(allFields.end(),physicsBlock.getProvidedDOFs().begin(),physicsBlock.getProvidedDOFs().end());

  // Add in tangent fields
  if(addSolutionFields_)
    allFields.insert(allFields.end(),physicsBlock.getTangentFields().begin(),physicsBlock.getTangentFields().end());

  // add in bases for any addtional fields
  for(std::size_t i=0;i<additionalFields_.size();i++)
    bases[additionalFields_[i].second->name()] = additionalFields_[i].second;

  allFields.insert(allFields.end(),additionalFields_.begin(),additionalFields_.end());

  deleteRemovedFields(removedFields_,allFields);

  bucketByBasisType(allFields,basisBucket);

  // add this for HCURL and HDIV basis, only want to add them once: evaluate vector fields at centroid
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////
  RCP<panzer::PointRule> centroidRule;
  for(std::map<std::string,Teuchos::RCP<const panzer::PureBasis> >::const_iterator itr=bases.begin();
      itr!=bases.end();++itr) {

    if(itr->second->isVectorBasis()) {
      centroidRule = rcp(new panzer::PointRule("Centroid",1,physicsBlock.cellData()));

      // compute centroid
      Kokkos::DynRankView<double,PHX::Device> centroid;
      computeReferenceCentroid(bases,physicsBlock.cellData().baseCellDimension(),centroid);

      // build pointe values evaluator
      RCP<PHX::Evaluator<panzer::Traits> > evaluator  = 
         rcp(new panzer::PointValues_Evaluator<EvalT,panzer::Traits>(centroidRule,centroid));
      this->template registerEvaluator<EvalT>(fm, evaluator);

      break; // get out of the loop, only need one evaluator
    }
  }

  // add evaluators for each field
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////

  for(std::map<std::string,std::vector<std::string> >::const_iterator itr=basisBucket.begin();
      itr!=basisBucket.end();++itr) {

    std::string basisName = itr->first;
    const std::vector<std::string> & fields = itr->second;

    std::map<std::string,Teuchos::RCP<const panzer::PureBasis> >::const_iterator found = bases.find(basisName);
    TEUCHOS_TEST_FOR_EXCEPTION(found==bases.end(),std::logic_error,
                               "Could not find basis \""+basisName+"\"!");
    Teuchos::RCP<const panzer::PureBasis> basis = found->second;

    // determine if user has modified field scalar for each field to be written to STK
    std::vector<double> scalars(fields.size(),1.0); // fill with 1.0
    for(std::size_t f=0;f<fields.size();f++) {
      std::unordered_map<std::string,double>::const_iterator f2s_itr = fieldToScalar_.find(fields[f]);

      // if scalar is found, include it in the vector and remove the field from the
      // hash table so it won't be included in the warning message.
      if(f2s_itr!=fieldToScalar_.end()) {
        scalars[f] = f2s_itr->second;
        scaledFieldsHash.erase(fields[f]);
      }
    }

    // write out nodal fields
    if(basis->getElementSpace()==panzer::PureBasis::HGRAD ||
       basis->getElementSpace()==panzer::PureBasis::CONST) {

      std::string fields_concat = "";
      for(std::size_t f=0;f<fields.size();f++) {
        fields_concat += fields[f];
      }

      Teuchos::RCP<PHX::Evaluator<panzer::Traits> > eval = 
        Teuchos::rcp(new ScatterFields<EvalT,panzer::Traits>("STK HGRAD Scatter Basis " +basis->name()+": "+fields_concat,
                                                      mesh_, basis, fields,scalars));

      // register and require evaluator fields
      this->template registerEvaluator<EvalT>(fm, eval);
      fm.template requireField<EvalT>(*eval->evaluatedFields()[0]);
    }
    else if(basis->getElementSpace()==panzer::PureBasis::HCURL) {
      TEUCHOS_ASSERT(centroidRule!=Teuchos::null);

      // register basis values evaluator
      {
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer::BasisValues_Evaluator<EvalT,panzer::Traits>(centroidRule,basis));
        this->template registerEvaluator<EvalT>(fm, evaluator);
      }

      // add a DOF_PointValues for each field
      std::string fields_concat = "";
      for(std::size_t f=0;f<fields.size();f++) {
        Teuchos::ParameterList p;
        p.set("Name",fields[f]);
        p.set("Basis",basis);
        p.set("Point Rule",centroidRule.getConst());
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer::DOF_PointValues<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, evaluator);

        fields_concat += fields[f];
      }

      // add the scatter field evaluator for this basis
      {
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new ScatterVectorFields<EvalT,panzer::Traits>("STK HCURL Scatter Basis " +basis->name()+": "+fields_concat,
                                                                        mesh_,centroidRule,fields,scalars));

        this->template registerEvaluator<EvalT>(fm, evaluator);
        fm.template requireField<EvalT>(*evaluator->evaluatedFields()[0]); // require the dummy evaluator
      }
    }
    else if(basis->getElementSpace()==panzer::PureBasis::HDIV) {
      TEUCHOS_ASSERT(centroidRule!=Teuchos::null);

      // register basis values evaluator
      {
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer::BasisValues_Evaluator<EvalT,panzer::Traits>(centroidRule,basis));
        this->template registerEvaluator<EvalT>(fm, evaluator);
      }

      // add a DOF_PointValues for each field
      std::string fields_concat = "";
      for(std::size_t f=0;f<fields.size();f++) {
        Teuchos::ParameterList p;
        p.set("Name",fields[f]);
        p.set("Basis",basis);
        p.set("Point Rule",centroidRule.getConst());
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new panzer::DOF_PointValues<EvalT,panzer::Traits>(p));

        this->template registerEvaluator<EvalT>(fm, evaluator);

        fields_concat += fields[f];
      }

      // add the scatter field evaluator for this basis
      {
        Teuchos::RCP<PHX::Evaluator<panzer::Traits> > evaluator  
           = Teuchos::rcp(new ScatterVectorFields<EvalT,panzer::Traits>("STK HDIV Scatter Basis " +basis->name()+": "+fields_concat,
                                                                        mesh_,centroidRule,fields,scalars));

        this->template registerEvaluator<EvalT>(fm, evaluator);
        fm.template requireField<EvalT>(*evaluator->evaluatedFields()[0]); // require the dummy evaluator
      }
    }
  }

  // print warning message for any unused scaled fields
  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);

  for(std::unordered_set<std::string>::const_iterator itr=scaledFieldsHash.begin();
      itr!=scaledFieldsHash.end();itr++) { 
    out << "WARNING: STK Solution Writer did not scale the field \"" << *itr << "\" "
        << "because it was not written." << std::endl;
  }
}

template <typename EvalT>
void ResponseEvaluatorFactory_SolutionWriter<EvalT>::
bucketByBasisType(const std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > & providedDofs,
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
computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<const panzer::PureBasis> > & bases,
                         int baseDimension,
                         Kokkos::DynRankView<double,PHX::Device> & centroid) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   centroid = Kokkos::DynRankView<double,PHX::Device>("centroid",1,baseDimension);
   auto l_centroid = centroid;

   // loop over each possible basis
   for(std::map<std::string,RCP<const panzer::PureBasis> >::const_iterator itr=bases.begin();
       itr!=bases.end();++itr) {

      RCP<Intrepid2::Basis<PHX::exec_space,double,double>> intrepidBasis = itr->second->getIntrepid2Basis();

      // we've got coordinates, lets commpute the "centroid"
      Kokkos::DynRankView<double,PHX::Device> coords("coords",intrepidBasis->getCardinality(),
						     intrepidBasis->getBaseCellTopology().getDimension());
      intrepidBasis->getDofCoords(coords);
      TEUCHOS_ASSERT(coords.rank()==2);
      TEUCHOS_ASSERT(coords.extent_int(1)==baseDimension);
      
      Kokkos::parallel_for(coords.extent_int(0), KOKKOS_LAMBDA (int i)  {
	  for(int d=0;d<coords.extent_int(1);d++) 
            l_centroid(0,d) += coords(i,d)/coords.extent(0);
	});
      Kokkos::fence();

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
  scaledFieldsHash_.insert(fieldName);
}

template <typename EvalT>
bool ResponseEvaluatorFactory_SolutionWriter<EvalT>::
typeSupported() const
{
  if(PHX::print<EvalT>()==PHX::print<panzer::Traits::Residual>())
    return true;

  return false;
}

template <typename EvalT>
void ResponseEvaluatorFactory_SolutionWriter<EvalT>::
addAdditionalField(const std::string & fieldName,const Teuchos::RCP<const panzer::PureBasis> & basis)
{
  additionalFields_.push_back(std::make_pair(fieldName,basis));
}

template <typename EvalT>
void ResponseEvaluatorFactory_SolutionWriter<EvalT>::
deleteRemovedFields(const std::vector<std::string> & removedFields,
                    std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > & fields) const
{
  RemovedFieldsSearchUnaryFunctor functor;
  functor.removedFields_ = removedFields;

  // This is the Erase-Remove Idiom: see http://en.wikipedia.org/wiki/Erase-remove_idiom
  fields.erase(std::remove_if(fields.begin(),fields.end(),functor),fields.end());
}

}

#endif
