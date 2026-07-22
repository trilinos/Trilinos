// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_hpp__
#define __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_hpp__

#include <string>

#include "PanzerAdaptersSTK_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"

#include "Panzer_STK_Interface.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <unordered_map>
#include <unordered_set>

namespace panzer_stk {

/** This class defines a response based solution writer.
  */
template <typename EvalT> 
class ResponseEvaluatorFactory_SolutionWriter : public panzer::ResponseEvaluatorFactory<EvalT> {
public:

   ResponseEvaluatorFactory_SolutionWriter(const Teuchos::RCP<STK_Interface> & mesh)
     : mesh_(mesh), addSolutionFields_(true), addCoordinateFields_(true) {}

   virtual ~ResponseEvaluatorFactory_SolutionWriter() {}
 
   /** Build the response object used by this factory. This object
     * assumes the role of the scatter target and will be accessible
     * by all the evaluators in the field managers. 
     *
     * \param[in] responseName Name of response to be built. This
     *                         name will be used for looking up
     *                         the response in the <code>GlobalEvaluationDataContainer</code>
     *                         object.
     */
   virtual Teuchos::RCP<panzer::ResponseBase> buildResponseObject(const std::string & responseName) const;

   virtual Teuchos::RCP<panzer::ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<panzer::WorksetDescriptor>& /* wkstDesc */) const 
   { return buildResponseObject(responseName); }
   
   /** Build and register evaluators for a response on a particular physics
     * block. 
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be fuild with the evaluators.
     * \param[in] physicsBlock What physics block is being used for constructing
     *                         the evaluators
     * \param[in] user_data The user data parameter list, this stores things
     *                      that the user may find useful.
     */
   virtual void buildAndRegisterEvaluators(const std::string & responseName,
                                           PHX::FieldManager<panzer::Traits> & fm,
                                           const panzer::PhysicsBlock & physicsBlock,
                                           const Teuchos::ParameterList & user_data) const;

   /** Is this evaluation type supported by the factory. This is used to determine cases
     * where a response may support a particular evaluation type, however at runtime the user
     * decides not to enable the (say) Jacobian evaluation of this response.
     *
     * Note that use of this mechanism is complementary to having the builder return 
     * <code>Teuchos::null</code> for a particular evaluation type.
     */
   virtual bool typeSupported() const;

   /** Take a vector of (std::string (field name), RCP<PureBasis>) pairs and bucket them
     * by basis name. What is returned is a map pairing the basis to a vector of field names.
     */
   static void bucketByBasisType(const std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > & providedDofs,
                                 std::map<std::string,std::vector<std::string> > & basisBucket);

   /** Scale a field before writing out to STK by a prescribed value (the implicit default is 1.0). 
     * A warning will be printed if no field is found of that name when <code>buildAndRegisterEvaluators</code>
     * is called.
     *
     * \note Currently only works for HGRAD fields.
     *
     * \param[in] fieldName HGRAD field to scale.
     * \param[in] fieldScalar Value to scale the field by.
     */
   void scaleField(const std::string & fieldName,double fieldScalar);

   /** Add an additional (solution) field to write out that is not in the physics blocks.
     */
   void addAdditionalField(const std::string & fieldName,const Teuchos::RCP<const panzer::PureBasis> & basis);

   /** Enable/disable addition of solution fields. Note that this "true" by default.
     */
   void setAddSolutionFields(bool asf) 
   { addSolutionFields_ = asf; }

   /** Enable/disable addition of coordinate fields. Note that this "true" by default.
     */
   void setAddCoordinateFields(bool acf) 
   { addCoordinateFields_ = acf; }

   /** Remove a field (even a solution field) from the response. Note that even if a field has not
     * been added, it will be removed by any previous or following call to removeField. This implies
     * that removeField takes precedence.
     */
   void removeField(const std::string & fieldName)
   { removedFields_.push_back(fieldName); }

  // should be private but needs a lambda
   void computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<const panzer::PureBasis> > & bases,
                                 int baseDimension,
                                 Kokkos::DynRankView<double,PHX::Device> & centroid) const;

private:
   //! Delete from the argument all the fields that are in the removedFields array
   void deleteRemovedFields(const std::vector<std::string> & removedFields,
                            std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > & fields) const;

   struct RemovedFieldsSearchUnaryFunctor {
     std::vector<std::string> removedFields_;
     bool operator() (const std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > & field) 
     { return std::find(removedFields_.begin(),removedFields_.end(),field.first)!=removedFields_.end(); }
   };

   Teuchos::RCP<STK_Interface> mesh_;

   std::unordered_map<std::string,double> fieldToScalar_;
   std::unordered_set<std::string> scaledFieldsHash_; // used to print the warning about unused scaling

   std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > additionalFields_;
   std::vector<std::string> removedFields_;
   bool addSolutionFields_;
   bool addCoordinateFields_;
};

/** A simple builder for this the SolutionWriter response factory, simply set the mesh 
  * and this will build the response factories for you. (Pass into ResponseLibrary::addResponse)
  */
struct RespFactorySolnWriter_Builder {
  RespFactorySolnWriter_Builder() : addSolutionFields_(true), addCoordinateFields_(true) {}

  Teuchos::RCP<panzer_stk::STK_Interface> mesh;

  void scaleField(const std::string & fieldName,double fieldScalar)
  { fieldToScalar_[fieldName] = fieldScalar; }

  void addAdditionalField(const std::string & fieldName,const Teuchos::RCP<const panzer::PureBasis> & basis)
  { additionalFields_.push_back(std::make_pair(fieldName,basis)); }

  /** Remove a field (even a solution field) from the response. Note that even if a field has not
    * been added, it will be removed by any previous or following call to removeField. This implies
    * that removeField takes precedence.
    */
  void removeField(const std::string & fieldName)
  { removedFields_.push_back(fieldName); }

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { 
    Teuchos::RCP<ResponseEvaluatorFactory_SolutionWriter<T> > ref = 
        Teuchos::rcp(new panzer_stk::ResponseEvaluatorFactory_SolutionWriter<T>(mesh)); 
 
    // disable/enable the solution fields
    ref->setAddSolutionFields(addSolutionFields_);

    // disable/enable the coordinate fields
    ref->setAddCoordinateFields(addCoordinateFields_);

    // add all additional fields
    for(std::size_t i=0;i<additionalFields_.size();i++)
      ref->addAdditionalField(additionalFields_[i].first,additionalFields_[i].second);

    for(std::size_t i=0;i<removedFields_.size();i++)
      ref->removeField(removedFields_[i]);

    // set all scaled field values
    for(std::unordered_map<std::string,double>::const_iterator itr=fieldToScalar_.begin();
        itr!=fieldToScalar_.end();++itr) 
      ref->scaleField(itr->first,itr->second);

    return ref;
  }

   /** Enable/disable addition of solution fields. Note that this "true" by default.
     */
   void setAddSolutionFields(bool asf) 
   { addSolutionFields_ = asf; }

   /** Enable/disable addition of coordinate fields. Note that this "true" by default.
     */
   void setAddCoordinateFields(bool acf) 
   { addCoordinateFields_ = acf; }

private:
  std::unordered_map<std::string,double> fieldToScalar_;
  std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > additionalFields_;
  std::vector<std::string> removedFields_;
  bool addSolutionFields_;
  bool addCoordinateFields_;
};

}

#endif
