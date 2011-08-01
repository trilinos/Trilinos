// FEI includes
#include "fei_Factory_Trilinos.hpp"

#include <map>

#include "Panzer_GeometricAggFieldPattern.hpp"

#include "Teuchos_DefaultMpiComm.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class DOFManager
// ************************************************************

template <typename LocalOrdinalT,typename GlobalOrdinalT>
DOFManager<LocalOrdinalT,GlobalOrdinalT>::DOFManager()
   : numFields_(0), fieldsRegistered_(false)
{ }

template <typename LocalOrdinalT,typename GlobalOrdinalT>
DOFManager<LocalOrdinalT,GlobalOrdinalT>::DOFManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm)
   : numFields_(0), fieldsRegistered_(false)
{
   setConnManager(connMngr,mpiComm);
}

/** \brief Set the connection manager and MPI_Comm objects.
  *
  * Set the connection manager and MPI_Comm objects. If this method
  * is called more than once, the behavior is to reset the indices in
  * the DOF manager.  However, the fields will be the same (this assumes
  * that the element blocks are consistent with the fields). The indices
  * will need to be rebuilt by calling <code>buildGlobalUnknowns</code>.
  *
  * \param[in] connMngr Connection manager to use.
  * \param[in] mpiComm  Communicator to use.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::setConnManager(const Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm)
{
   // this kills any old connection manager as well as the old FEI objects
   resetIndices();

   connMngr_ = connMngr;

   // build fei components
   feiFactory_ = Teuchos::rcp(new Factory_Trilinos(mpiComm));

   // build fei components
   vectorSpace_ = feiFactory_->createVectorSpace(mpiComm,"problem_vs");
   matrixGraph_ = feiFactory_->createMatrixGraph(vectorSpace_,vectorSpace_,"problem_mg");

   // define a single id type: node
   nodeType_ = 0; 
   vectorSpace_->defineIDTypes(1,&nodeType_);
   edgeType_ = 1; 
   vectorSpace_->defineIDTypes(1,&edgeType_);

   communicator_ = Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(mpiComm)));
}

/** \brief Reset the indicies for this DOF manager.
  *
  * This method resets the indices and wipes out internal state. This method
  * does preserve the fields and the patterns added to the object. Also the
  * old connection manager is returned.
  *
  * \returns Old connection manager.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > DOFManager<LocalOrdinalT,GlobalOrdinalT>::resetIndices()
{
   Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr = connMngr_;

   connMngr_ = Teuchos::null;

   // wipe out FEI objects
   patternNum_.clear();
   feiFactory_ = Teuchos::null;
   vectorSpace_.reset();
   matrixGraph_.reset();

   ownedGIDHashTable_.clear(); 

   return connMngr;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::addField(const std::string & str,
                                                        const Teuchos::RCP<const FieldPattern> & pattern)
{
   std::vector<std::string> elementBlockIds;
   connMngr_->getElementBlockIds(elementBlockIds);

   // loop over blocks adding field pattern to each 
   for(std::size_t i=0;i<elementBlockIds.size();i++)
      addField(elementBlockIds[i],str,pattern);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::addField(const std::string & blockId,const std::string & str,
                                                        const Teuchos::RCP<const FieldPattern> & pattern)
{
   TEST_FOR_EXCEPTION(fieldsRegistered_,std::logic_error,
                      "DOFManager::addField: addField cannot be called after registerFields or"
                      "buildGlobalUnknowns has been called"); 

   fieldStringToPattern_[std::make_pair(blockId,str)] = pattern;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::registerFields() 
{
   numFields_ = 0;

   // test validity of the field order
   {
      // build a unique set of fields, so we can compare validate the ordered list
      std::set<std::string> fields;
      for(std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> >::const_iterator
          fieldItr=fieldStringToPattern_.begin(); fieldItr!=fieldStringToPattern_.end();++fieldItr) {
         std::string fieldName = fieldItr->first.second;
         fields.insert(fieldName);
      }

      // construct default field order if neccessary
      if(fieldOrder_.size()==0) {
         std::set<std::string>::const_iterator itr;
         for(itr=fields.begin();itr!=fields.end();itr++)
            fieldOrder_.push_back(*itr);
      }

      // check validity of field order: no repeats, and everything is accounted for
      bool validOrder = validFieldOrder(fieldOrder_,fields);
      if(!validOrder) {
         // for outputing
         std::stringstream ss;

         ss << "DOFManager::registerFields - Field order is invalid!\n";

         ss << "   fields = [ ";
         for(std::set<std::string>::const_iterator itr=fields.begin();
             itr!=fields.end();++itr)
            ss << "\"" << *itr << "\" ";
         ss << " ]\n";

         ss << "   fieldOrder = [ ";
         for(std::vector<std::string>::const_iterator itr=fieldOrder_.begin();
             itr!=fieldOrder_.end();++itr)
            ss << "\"" << *itr << "\" ";
         ss << " ]\n";

         TEST_FOR_EXCEPTION(!validOrder,std::logic_error,ss.str());
      }
   }

   // build field IDs
   for(std::size_t fo_index=0;fo_index<fieldOrder_.size();fo_index++) {
      std::string fieldName = fieldOrder_[fo_index];

      // field doesn't exist...add it
      int fieldNum = fo_index;
      int size = 1; // fields are always size 1
      vectorSpace_->defineFields(1,&fieldNum,&size);

      fieldStrToInt_[fieldName] = fieldNum;
      intToFieldStr_[fieldNum] = fieldName;
   }      
   numFields_ = fieldOrder_.size();

   // initialize blockToField_ vector to have at least empty sets
   // for each element block
   std::vector<std::string> elementBlockIds;
   getElementBlockIds(elementBlockIds);
   for(std::size_t ebi=0;ebi<elementBlockIds.size();ebi++)
      blockToField_.insert(std::make_pair(elementBlockIds[ebi],std::set<int>())); 

   // associate blocks with particular field ids
   for(std::map<std::pair<std::string,std::string>,Teuchos::RCP<const FieldPattern> >::const_iterator
       fieldItr=fieldStringToPattern_.begin(); fieldItr!=fieldStringToPattern_.end();++fieldItr) {
 
      std::string blockId = fieldItr->first.first;
      std::string fieldName = fieldItr->first.second;

      std::map<std::string,int>::const_iterator itr = fieldStrToInt_.find(fieldName);
      if(itr!=fieldStrToInt_.end()) {
         // field already exists!
         blockToField_[blockId].insert(itr->second); 
         fieldIntToPattern_[std::make_pair(blockId,itr->second)] = fieldItr->second;
      }
      else {
         // this statement should _never_ be executed. The reason is that
         // the fieldIntToPattern_ was filled before this function was run
         // directly from the fieldStringToPattern_ map. Possibly check the
         // order validator for letting something slip through!
         
         TEST_FOR_EXCEPTION(false,std::logic_error,
                            "DOFManager::registerFields - Impossible case discoverved!");
      }
   }

   fieldsRegistered_ = true;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int DOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldNum(const std::string & str) const
{
   std::map<std::string,int>::const_iterator itr = fieldStrToInt_.find(str);

   // return based on what was found
   if(itr==fieldStrToInt_.end()) {
      // incorrect field name
      TEST_FOR_EXCEPTION(true,std::logic_error,
                         "DOFManager::getFieldNum No field with the name \"" + str + "\" has been added");
   }
   else {
      return itr->second;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::setFieldOrder(const std::vector<std::string> & fieldOrder)
{
   fieldOrder_ = fieldOrder;
}

/** Get the field order used. Return the field strings.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldOrder(std::vector<std::string> & fieldOrder) const
{
   fieldOrder = fieldOrder_;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int DOFManager<LocalOrdinalT,GlobalOrdinalT>::getNumFields() const
{
   return vectorSpace_->getNumFields();
}

// build the global unknown numberings
//   1. this builds the pattens
//   2. initializes the connectivity
//   3. calls initComplete
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGlobalUnknowns(const Teuchos::RCP<const FieldPattern> & geomPattern)
{
   if(!fieldsRegistered_)
      registerFields();

   std::vector<std::string> fieldOrder;
   getFieldOrder(fieldOrder);

   Teuchos::RCP<const ConnManager<LocalOrdinalT,GlobalOrdinalT> > connMngr = connMngr_.getConst();
   geomPattern_ = geomPattern;

   // get element blocks
   std::vector<std::string> elementBlockIds;
   connMngr->getElementBlockIds(elementBlockIds);

   // setup connectivity mesh
   patternNum_.resize(connMngr->numElementBlocks()); 
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      std::size_t blockIndex = blockIdToIndex(blockId);

      // build the pattern
      bool patternBuilt = buildPattern(fieldOrder,blockId);

      if(patternBuilt) {
         // note that condition "patternBuilt==true" implies "fieldAggPattern_[blockId]!=Teuchos::null"

         // figure out what IDs are active for this pattern
         const std::vector<int> & numFieldsPerID = fieldAggPattern_[blockId]->numFieldsPerId();
         std::vector<int> activeIds;
         for(std::size_t i=0;i<numFieldsPerID.size();i++)
            if(numFieldsPerID[i]>0) 
               activeIds.push_back(i);
         std::vector<int> reduceConn(activeIds.size()); // which IDs to use
      
         // grab elements for this block
         const std::vector<LocalOrdinal> & elements = connMngr->getElementBlock(blockId);
   
         // build graph for this block
         matrixGraph_->initConnectivityBlock(blockIndex,elements.size(),patternNum_[blockIndex]);
         for(std::size_t e=0;e<elements.size();e++) {
            const GlobalOrdinal * conn = connMngr->getConnectivity(elements[e]);
            for(std::size_t i=0;i<activeIds.size();i++)
               reduceConn[i] = conn[activeIds[i]];
    
            matrixGraph_->initConnectivity(blockIndex,elements[e],&reduceConn[0]);
         }
      }
      // else: no fields on this block, don't try to do anything with it.
      //       basically no field has been added that is associated with this
      //       element block. This is OK, but we need to correctly ignore it.
   }
   matrixGraph_->initComplete();

   // build owned map
   std::vector<GlobalOrdinal> ownedIndices;
   getOwnedIndices(ownedIndices);
   ownedGIDHashTable_.insert(ownedIndices.begin(),ownedIndices.end());  
}

// build the global unknown numberings
//   1. this builds the pattens
//   2. initializes the connectivity
//   3. calls initComplete
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGlobalUnknowns()
{
   if(!fieldsRegistered_)
      registerFields();

   // build the pattern for the ID layout on the mesh
   std::vector<RCP<const FieldPattern> > patVector;
   RCP<GeometricAggFieldPattern> aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern);;
   std::map<std::pair<std::string,int>,Teuchos::RCP<const FieldPattern> >::iterator f2p_itr;
   for(f2p_itr=fieldIntToPattern_.begin();f2p_itr!=fieldIntToPattern_.end();f2p_itr++)
      patVector.push_back(f2p_itr->second);
   aggFieldPattern->buildPattern(patVector);

   // setup connectivity mesh
   connMngr_->buildConnectivity(*aggFieldPattern);

   // using new geometric pattern, build global unknowns
   buildGlobalUnknowns(aggFieldPattern);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOrderedBlock(const std::vector<std::string> & fieldOrder,
                                                               const std::string & blockId,
                                                               std::vector<int> & orderedBlock) const
{
   const std::set<int> & fieldSet = this->getFields(blockId);
   orderedBlock.clear();

   std::vector<std::string>::const_iterator itr;
   for(itr=fieldOrder.begin();itr!=fieldOrder.end();++itr) {
      int fieldNum = this->getFieldNum(*itr);

      // if field in in a particular block add it 
      if(fieldSet.find(fieldNum)!=fieldSet.end())
         orderedBlock.push_back(fieldNum);
   }
}

// build the pattern associated with this manager
template <typename LocalOrdinalT,typename GlobalOrdinalT>
bool DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildPattern(const std::vector<std::string> & fieldOrder,
                                                            const std::string & blockId)
{
   using Teuchos::rcp;
   using Teuchos::RCP;

   // use some generic field ordering if the current one is empty
   std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > blockPatterns;
   std::vector<int> orderedBlock;
   getOrderedBlock(fieldOrder,blockId,orderedBlock);

   // get a map of field patterns
   std::vector<int>::const_iterator itr;
   for(itr=orderedBlock.begin();itr!=orderedBlock.end();++itr) {
      Teuchos::RCP<const FieldPattern> fp = fieldIntToPattern_[std::make_pair(blockId,*itr)];
      blockPatterns.push_back(std::make_pair(*itr,fp));
   }

   std::size_t blockIndex = blockIdToIndex(blockId);
   if(blockPatterns.size()<=0) {
      patternNum_[blockIndex] = -1; // this should cause an error
                                    // if this is accidently visited
      return false;
   }
   
   // smash together all fields...do interlacing
   fieldAggPattern_[blockId] = rcp(new FieldAggPattern(blockPatterns));

   // build FEI pattern
   const std::vector<int> & fields = fieldAggPattern_[blockId]->fieldIds();
   const std::vector<int> & numFieldsPerID = fieldAggPattern_[blockId]->numFieldsPerId();

   std::vector<int> reduceNumFieldsPerID;
   for(std::size_t i=0;i<numFieldsPerID.size();i++)
      if(numFieldsPerID[i]>0) 
         reduceNumFieldsPerID.push_back(numFieldsPerID[i]);

   int idsPerElement  = reduceNumFieldsPerID.size();

   patternNum_[blockIndex] 
         = matrixGraph_->definePattern(idsPerElement,nodeType_,&reduceNumFieldsPerID[0],&fields[0]);

   return true;
}

// "Get" functions
/////////////////////////////////////////////////////////////////////

namespace {
// hide these functions

template <typename LocalOrdinalT,typename GlobalOrdinalT>
inline void getGIDsFromMatrixGraph(int blockIndex,int dof,LocalOrdinalT localElmtId,
                                   fei::MatrixGraph & mg,std::vector<GlobalOrdinalT> & gids)
{
   std::vector<int> indices(dof);

   // get elements indices
   int localSize = -1;
   mg.getConnectivityIndices(blockIndex,localElmtId,dof,&indices[0],localSize);
   
   // copy the indices
   gids.resize(dof);
   for(std::size_t i=0;i<indices.size();i++)
      gids[i] = (GlobalOrdinalT) indices[i];
}

template <typename LocalOrdinalT>
inline void getGIDsFromMatrixGraph(int blockIndex,int dof,LocalOrdinalT localElmtId,
                                   fei::MatrixGraph & mg,std::vector<int> & gids)
{
   // get elements indices
   gids.resize(dof);
   int localSize = -1;
   mg.getConnectivityIndices(blockIndex,localElmtId,dof,&gids[0],localSize);
}

}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids) const
{
   // get information about number of indicies
   std::string blockId = connMngr_->getBlockId(localElmtId);
   std::size_t blockIndex = blockIdToIndex(blockId);
   int dof = getElementBlockGIDCount(blockIndex);

   // getConnectivityNumIndices returns -1 if no block is found or
   // has been initialized. So if this DOFManager has no fields on
   // the block it should be ignored (hence dof>0)
   if(dof>0) {

      getGIDsFromMatrixGraph(blockIndex,dof,localElmtId,*matrixGraph_,gids);
   }
}

/** \brief Get a vector containg the orientation of the GIDs relative to the neighbors.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::
getElementOrientation(LocalOrdinalT localElmtId,std::vector<double> & gidsOrientation) const
{
   TEST_FOR_EXCEPTION(true,std::logic_error,"DOFManager::getElementOrientation not implemented yet!");
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::printFieldInformation(std::ostream & os) const
{
   os << "DOFManager Field Information: " << std::endl;
   
   if(fieldsRegistered_) {
      std::map<std::string,Teuchos::RCP<FieldAggPattern> >::const_iterator iter;
      for(iter=fieldAggPattern_.begin();iter!=fieldAggPattern_.end();++iter) {
         os << "Element Block = " << iter->first << std::endl; 
         iter->second->print(os);
   
         // output field information
         std::set<int>::const_iterator itr_fieldIds = blockToField_.find(iter->first)->second.begin(); 
         std::set<int>::const_iterator end_fieldIds = blockToField_.find(iter->first)->second.end(); 
         os << "   Field String to Field Id:\n";
         for( /*empty*/ ;itr_fieldIds!=end_fieldIds;++itr_fieldIds)
            os << "      \"" << getFieldString(*itr_fieldIds) << "\" is field ID " << *itr_fieldIds << std::endl;
      }
   }
   else {
      // fields are not registered
      os << "Fields not yet registered! Unknowns not built (call buildGlobalUnknowns)" << std::endl;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<const FieldPattern> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldPattern(const std::string & blockId, int fieldNum) const
{
   std::map<std::pair<std::string,int>,Teuchos::RCP<const FieldPattern> >::const_iterator itr;
   itr = fieldIntToPattern_.find(std::make_pair(blockId,fieldNum));

   if(itr==fieldIntToPattern_.end()) {
      // could not find requiested field pattern...return null
      return Teuchos::null;
   }

   return itr->second;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<const FieldPattern> DOFManager<LocalOrdinalT,GlobalOrdinalT>::
getFieldPattern(const std::string & blockId, const std::string & fieldName) const
{
   // return null even if field doesn't exist in manager
   int fieldNum = -1;
   try {
      fieldNum = getFieldNum(fieldName);
   }
   catch(const std::logic_error & le) {
      return Teuchos::null;
   }

   return getFieldPattern(blockId,fieldNum);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
std::size_t DOFManager<LocalOrdinalT,GlobalOrdinalT>::blockIdToIndex(const std::string & blockId) const
{
   // use lazy evaluation to build block indices
   if(blockIdToIndex_==Teuchos::null) {

      std::vector<std::string> elementBlockIds;
      connMngr_->getElementBlockIds(elementBlockIds);

      // build ID to Index map
      blockIdToIndex_ = Teuchos::rcp(new std::map<std::string,std::size_t>);
      for(std::size_t i=0;i<elementBlockIds.size();i++)
         (*blockIdToIndex_)[elementBlockIds[i]] = i;
   }
 
   return (*blockIdToIndex_)[blockId];
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::ownedIndices(const std::vector<GlobalOrdinalT> & indices,std::vector<bool> & isOwned) const
{
   // verify size is correct
   if(indices.size()!=isOwned.size())
      isOwned.resize(indices.size(),false);

   // use unordered set to check for ownership of the ID
   typename boost::unordered_set<GlobalOrdinal>::const_iterator endItr = ownedGIDHashTable_.end();
   for(std::size_t i=0;i<indices.size();i++) 
      isOwned[i] = (ownedGIDHashTable_.find(indices[i])!=endItr);
}

// These two functions are "helpers" for DOFManager::getOwnedIndices
///////////////////////////////////////////////////////////////////////////
template <typename OrdinalType> 
static void getOwnedIndices_T(const fei::SharedPtr<fei::VectorSpace> & vs,std::vector<OrdinalType> & indices) 
{
   int numIndices, ni;
   numIndices = vs->getNumIndices_Owned();
   indices.resize(numIndices);
   std::vector<int> int_Indices; // until FEI is templated

   // get the number of locally owned degrees of freedom...allocate space
   int_Indices.resize(numIndices);

   // get the global indices
   vs->getIndices_Owned(numIndices,&int_Indices[0],ni);

   for(std::size_t i=0;i<int_Indices.size();i++) 
      indices[i] = (OrdinalType) int_Indices[i];
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const
{
   getOwnedIndices_T<GlobalOrdinalT>(vectorSpace_,indices);
}

/** Check the validity of a field order. This is used internally
  * as a sanity check. Checks for no repeats, bogus fields, and all fields
  * being included.
  *
  * \param[in] fieldOrder_ut Field order vector under test (ut).
  *
  * \returns true if the vector is valid, false otherwise.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
bool DOFManager<LocalOrdinalT,GlobalOrdinalT>::validFieldOrder(const std::vector<std::string> & fieldOrder_ut,const std::set<std::string> & fields) const
{
   if(fields.size()!=fieldOrder_ut.size()) // something is wrong!
      return false;

   std::set<std::string> fieldOrderSet;

   // first check the size by shoving everything into a set
   std::vector<std::string>::const_iterator itr;
   for(itr=fieldOrder_ut.begin();itr!=fieldOrder_ut.end();++itr)
      fieldOrderSet.insert(*itr);

   if(fieldOrderSet.size()!=fieldOrder_ut.size()) // there are repeat fields!
      return false;

   // check to make sure each field is represented
   std::set<std::string>::const_iterator itr_ut = fieldOrderSet.begin();
   std::set<std::string>::const_iterator itr_src = fields.begin();
   while(itr_ut!=fieldOrderSet.end()) {
      if(*itr_ut!=*itr_src) 
         return false;

      itr_ut++;
      itr_src++;
   }

   return true;
}

///////////////////////////////////////////////////////////////////////////

// These two functions are "helpers" for DOFManager::getOwnedAndSharedIndices
///////////////////////////////////////////////////////////////////////////
template <typename OrdinalType> 
static void getOwnedAndSharedIndices_T(const fei::SharedPtr<fei::VectorSpace> & vs,std::vector<OrdinalType> & indices) 
{
   std::vector<int> int_Indices; // until FEI is templated

   // get the global indices
   vs->getIndices_SharedAndOwned(int_Indices);

   indices.resize(int_Indices.size());
   for(std::size_t i=0;i<int_Indices.size();i++) 
      indices[i] = (OrdinalType) int_Indices[i];
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const
{
   getOwnedAndSharedIndices_T<GlobalOrdinalT>(vectorSpace_,indices);
}
///////////////////////////////////////////////////////////////////////////

}
