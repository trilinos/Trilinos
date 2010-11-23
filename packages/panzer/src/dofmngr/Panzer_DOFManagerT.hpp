// #include "Panzer_DOFManager.hpp"

// FEI includes
#include "fei_Factory_Trilinos.hpp"

#include <map>

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class DOFManager
// ************************************************************

template <typename LocalOrdinalT,typename GlobalOrdinalT>
DOFManager<LocalOrdinalT,GlobalOrdinalT>::DOFManager()
   : numFields_(0)
{ }

template <typename LocalOrdinalT,typename GlobalOrdinalT>
DOFManager<LocalOrdinalT,GlobalOrdinalT>::DOFManager(Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm)
   : numFields_(0)
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
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::setConnManager(Teuchos::RCP<ConnManager<LocalOrdinalT,GlobalOrdinalT> > & connMngr,MPI_Comm mpiComm)
{
   // this kills any old connection manager as well as the old FEI objects
   resetIndices();

   connMngr_ = connMngr;

   // build fei components
   feiFactory_ = Teuchos::rcp(new Factory_Trilinos(mpiComm));

   // build fei components
   comm_ = Teuchos::rcp(new Epetra_MpiComm(mpiComm));
   vectorSpace_ = feiFactory_->createVectorSpace(mpiComm,"problem_vs");
   matrixGraph_ = feiFactory_->createMatrixGraph(vectorSpace_,vectorSpace_,"problem_mg");

   // define a single id type: node
   nodeType_ = 0; 
   vectorSpace_->defineIDTypes(1,&nodeType_);
   edgeType_ = 1; 
   vectorSpace_->defineIDTypes(1,&edgeType_);
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

   // wipe out maps and graphs
   map_ = Teuchos::null;
   overlappedMap_ = Teuchos::null;
   graph_ = Teuchos::null;
   overlappedGraph_ = Teuchos::null;

   // wipe out FEI objects
   patternNum_.clear();
   feiFactory_ = Teuchos::null;
   comm_ = Teuchos::null;
   vectorSpace_.reset();
   matrixGraph_.reset();

   return connMngr;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::addField(const std::string & str,const Teuchos::RCP<const FieldPattern> & pattern)
{
   std::vector<std::string> elementBlockIds;
   connMngr_->getElementBlockIds(elementBlockIds);

   // loop over blocks adding field pattern to each 
   for(std::size_t i=0;i<elementBlockIds.size();i++)
      addField(elementBlockIds[i],str,pattern);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::addField(const std::string & blockId,const std::string & str,const Teuchos::RCP<const FieldPattern> & pattern)
{
   std::map<std::string,int>::const_iterator itr = fieldStrToInt_.find(str);
   if(itr!=fieldStrToInt_.end()) {
      // field already exists!
      blockToField_[blockId].insert(itr->second); 
      fieldIntToPattern_[std::make_pair(blockId,itr->second)] = pattern;
   }
   else {
      // field doesn't exist...add it
      int fieldNum = numFields_;
      int size = 1; // fields are always size 1
      vectorSpace_->defineFields(1,&fieldNum,&size);
      fieldStrToInt_[str] = fieldNum;
      intToFieldStr_[fieldNum] = str;
      fieldIntToPattern_[std::make_pair(blockId,fieldNum)] = pattern;
      blockToField_[blockId].insert(fieldNum); 
      numFields_++;
   }
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
int DOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldNum(const std::string & str) const
{
   std::map<std::string,int>::const_iterator itr = fieldStrToInt_.find(str);

   // return based on what was found
   if(itr==fieldStrToInt_.end())
      return -1;
   else
      return itr->second;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::setFieldOrder(const std::vector<int> & fieldOrder)
{
   fieldOrder_.clear();
   fieldOrder_ = fieldOrder;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::setFieldOrder(const std::vector<std::string> & fieldOrder)
{
   // convert to vector of field IDs...call integer version of fieldOrder_
   std::vector<int> fieldOrderInt;
   std::vector<std::string>::const_iterator strItr; 
   for(strItr=fieldOrder.begin();strItr!=fieldOrder.end();++strItr) {
      fieldOrderInt.push_back(getFieldNum(*strItr));
   }

   setFieldOrder(fieldOrderInt);
}

/** Get the field order used. Return the field IDs.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldOrder(std::vector<int> & fieldOrder) const
{
   fieldOrder = fieldOrder_; // just assign the field order
}

/** Get the field order used. Return the field strings.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getFieldOrder(std::vector<std::string> & fieldOrder) const
{
   // converge fieldOrder_ into a vector of strings
   std::vector<int>::const_iterator intItr;
   for(intItr=fieldOrder_.begin();intItr!=fieldOrder_.end();++intItr) {
      fieldOrder.push_back(getFieldString(*intItr));
   }
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
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGlobalUnknowns()
{
   // build the pattern for the ID layout on the mesh
   std::vector<RCP<const FieldPattern> > patVector;
   RCP<GeometricAggFieldPattern> aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern);;
   std::map<std::pair<std::string,int>,Teuchos::RCP<const FieldPattern> >::iterator f2p_itr;
   for(f2p_itr=fieldIntToPattern_.begin();f2p_itr!=fieldIntToPattern_.end();f2p_itr++)
      patVector.push_back(f2p_itr->second);
   aggFieldPattern->buildPattern(patVector);

   // get element blocks
   std::vector<std::string> elementBlockIds;
   connMngr_->getElementBlockIds(elementBlockIds);

   // setup connectivity mesh
   connMngr_->buildConnectivity(*aggFieldPattern);
   patternNum_.resize(connMngr_->numElementBlocks()); 
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      std::size_t blockIndex = blockIdToIndex(blockId);

      // build the pattern
      buildPattern(blockId,aggFieldPattern);

      // figure out what IDs are active for this pattern
      const std::vector<int> & numFieldsPerID = fieldAggPattern_[blockId]->numFieldsPerId();
      std::vector<int> activeIds;
      for(std::size_t i=0;i<numFieldsPerID.size();i++)
         if(numFieldsPerID[i]>0) 
            activeIds.push_back(i);
      std::vector<int> reduceConn(activeIds.size()); // which IDs to use
   
      // grab elements for this block
      const std::vector<LocalOrdinal> & elements = connMngr_->getElementBlock(blockId);

      // build graph for this block
      matrixGraph_->initConnectivityBlock(blockIndex,elements.size(),patternNum_[blockIndex]);
      for(std::size_t e=0;e<elements.size();e++) {
         const GlobalOrdinal * conn = connMngr_->getConnectivity(elements[e]);
         for(std::size_t i=0;i<activeIds.size();i++)
            reduceConn[i] = conn[activeIds[i]];
 
         matrixGraph_->initConnectivity(blockIndex,elements[e],&reduceConn[0]);
      }
   }
   matrixGraph_->initComplete();
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildDefaultFieldOrder()
{
   std::vector<int> fieldOrder;

   // build field order int vector from the ordering of field names
   std::map<std::string,int>::const_iterator s2iItr;
   for(s2iItr=fieldStrToInt_.begin();s2iItr!=fieldStrToInt_.end();++s2iItr) {
      fieldOrder.push_back(s2iItr->second);
   }

   // set the field order
   setFieldOrder(fieldOrder);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
std::vector<int> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOrderedBlock(const std::string & blockId)
{
   const std::set<int> & fieldSet = blockToField_[blockId];
   std::vector<int> orderedBlock;

   std::vector<int>::const_iterator itr;
   for(itr=fieldOrder_.begin();itr!=fieldOrder_.end();++itr) {
      // if field in in a particular block add it 
      if(fieldSet.find(*itr)!=fieldSet.end())
         orderedBlock.push_back(*itr);
   }

   return orderedBlock;
}

// build the pattern associated with this manager
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildPattern(const std::string & blockId,const RCP<FieldPattern> & geomPattern)
{
   using Teuchos::rcp;
   using Teuchos::RCP;

   // use some generic field ordering if the current one is empty
   if(fieldOrder_.size()==0)
      buildDefaultFieldOrder();

   std::vector<int> orderedBlock = getOrderedBlock(blockId);
   std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > blockPatterns;

   // get a map of field patterns
   // std::set<int>::const_iterator itr;
   std::vector<int>::const_iterator itr;
   for(itr=orderedBlock.begin();itr!=orderedBlock.end();++itr) {
      Teuchos::RCP<const FieldPattern> fp = fieldIntToPattern_[std::make_pair(blockId,*itr)];
      blockPatterns.push_back(std::make_pair(*itr,fp));
   }

   // smash together all fields...do interlacing
   fieldAggPattern_[blockId] = rcp(new FieldAggPattern(blockPatterns));
   TEUCHOS_ASSERT(geomPattern->equals(*fieldAggPattern_[blockId]->getGeometricAggFieldPattern()));

   // build FEI pattern
   const std::vector<int> & fields = fieldAggPattern_[blockId]->fieldIds();
   const std::vector<int> & numFieldsPerID = fieldAggPattern_[blockId]->numFieldsPerId();

   std::vector<int> reduceNumFieldsPerID;
   for(std::size_t i=0;i<numFieldsPerID.size();i++)
      if(numFieldsPerID[i]>0) 
         reduceNumFieldsPerID.push_back(numFieldsPerID[i]);

   int idsPerSimplex  = reduceNumFieldsPerID.size();

   std::size_t blockIndex = blockIdToIndex(blockId);
   patternNum_[blockIndex] 
         = matrixGraph_->definePattern(idsPerSimplex,nodeType_,&reduceNumFieldsPerID[0],&fields[0]);
}

// "Get" functions
/////////////////////////////////////////////////////////////////////

// get the map from the matrix
template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_Map> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getMap() const
{
   if(map_==Teuchos::null) map_ = buildMap();

   return map_;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_Map> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOverlapMap() const
{
   if(overlappedMap_==Teuchos::null) overlappedMap_ = buildOverlapMap();

   return overlappedMap_;
}

// get the graph of the crs matrix
template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getGraph() const
{
   if(graph_==Teuchos::null) graph_ = buildGraph();

   return graph_;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOverlapGraph() const
{
   if(overlappedGraph_==Teuchos::null) overlappedGraph_ = buildOverlapGraph();

   return overlappedGraph_;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_Import> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOverlapImport() const
{
   return Teuchos::rcp(new Epetra_Import(*getOverlapMap(),*getMap()));
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_Export> DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOverlapExport() const
{
   return Teuchos::rcp(new Epetra_Export(*getOverlapMap(),*getMap()));
}

// "Build" functions
/////////////////////////////////////////////////////////////////////

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_Map> DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildMap() const
{
   Teuchos::RCP<Epetra_Map> map; // result
   int numIndices,ni;
   std::vector<int> indices;

   // get the number of locally owned degrees of freedom...allocate space
   numIndices = vectorSpace_->getNumIndices_Owned();
   indices.resize(numIndices);

   // get the global indices
   vectorSpace_->getIndices_Owned(numIndices,&indices[0],ni);

   TEUCHOS_ASSERT(ni==numIndices); // sanity check
 
   map = Teuchos::rcp(new Epetra_Map(-1,numIndices,&indices[0],0,*comm_));

   return map;
}

// build the overlapped map
template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_Map> DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildOverlapMap() const
{
   std::vector<int> indices;

   // get the global indices
   vectorSpace_->getIndices_SharedAndOwned(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));
}

// get the graph of the crs matrix
template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildGraph() const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the overlapped graph
   RCP<Epetra_Map> map = getMap();
   RCP<Epetra_CrsGraph> graph  = rcp(new Epetra_CrsGraph(Copy,*map,0));
   RCP<Epetra_CrsGraph> oGraph = getOverlapGraph();

   // perform the communication to finish building graph
   RCP<Epetra_Export> exporter = getOverlapExport();
   graph->Export( *oGraph, *exporter, Insert );
   graph->FillComplete();
   return graph;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> DOFManager<LocalOrdinalT,GlobalOrdinalT>::buildOverlapGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> map = getOverlapMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*map,0));

   std::vector<std::string> elementBlockIds;
   connMngr_->getElementBlockIds(elementBlockIds);

   // graph information about the mesh
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      std::size_t blockIndex = blockIdToIndex(blockId);

      // grab elements for this block
      const std::vector<LocalOrdinal> & elements = connMngr_->getElementBlock(blockId);

      // get information about number of indicies
      int rowDOF=-1, colDOF=-1;
      matrixGraph_->getConnectivityNumIndices(blockIndex,rowDOF,colDOF); 
      std::vector<int> rowIndices(rowDOF);
      std::vector<int> colIndices(colDOF);

      // loop over the elemnts
      for(std::size_t i=0;i<elements.size();i++) {
         int localRDOF = -1,localCDOF = -1;

         // get elements indices
         matrixGraph_->getConnectivityIndices(blockIndex,elements[i],
                                                         rowDOF,&rowIndices[0],localRDOF,
                                                         colDOF,&colIndices[0],localCDOF);
       
         // sanity check
         TEUCHOS_ASSERT(localRDOF==rowDOF && localCDOF==colDOF);

         // insert these indicies into the graph
         for(int j=0;j<rowDOF;j++)
            graph->InsertGlobalIndices(rowIndices[j],colDOF,&colIndices[0]);
      }
   }

   // finish filling the graph
   graph->FillComplete();

   return graph;
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getElementGIDs(LocalOrdinalT localElmtId,std::vector<GlobalOrdinalT> & gids) const
{
   // get information about number of indicies
   std::string blockId = connMngr_->getBlockId(localElmtId);
   std::size_t blockIndex = blockIdToIndex(blockId);
   int dof = matrixGraph_->getConnectivityNumIndices(blockIndex);
   std::vector<int> indices(dof);

   // get elements indices
   int localSize = -1;
   matrixGraph_->getConnectivityIndices(blockIndex,localElmtId,dof,&indices[0],localSize);

   // copy the indices
   gids.resize(dof);
   for(std::size_t i=0;i<indices.size();i++)
      gids[i] = (GlobalOrdinal) indices[i];
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::printFieldInformation(std::ostream & os) const
{
   os << "DOFManager Field Information: " << std::endl;
   
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

template < >
static void getOwnedIndices_T<int>(const fei::SharedPtr<fei::VectorSpace> & vs,std::vector<int> & indices) 
{
   int numIndices, ni;
   numIndices = vs->getNumIndices_Owned();
   indices.resize(numIndices);

   // directly write to int indices
   vs->getIndices_Owned(numIndices,&indices[0],ni);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOwnedIndices(std::vector<GlobalOrdinalT> & indices) const
{
   getOwnedIndices_T<GlobalOrdinalT>(vectorSpace_,indices);
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

template < >
static void getOwnedAndSharedIndices_T<int>(const fei::SharedPtr<fei::VectorSpace> & vs,std::vector<int> & indices) 
{
   // get the global indices
   vs->getIndices_SharedAndOwned(indices);
}

template <typename LocalOrdinalT,typename GlobalOrdinalT>
void DOFManager<LocalOrdinalT,GlobalOrdinalT>::getOwnedAndSharedIndices(std::vector<GlobalOrdinalT> & indices) const
{
   getOwnedAndSharedIndices_T<GlobalOrdinalT>(vectorSpace_,indices);
}
///////////////////////////////////////////////////////////////////////////

}
