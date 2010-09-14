#include "Panzer_DOFManager.hpp"

// FEI includes
#include "fei_Factory_Trilinos.hpp"

#include <map>

#include "Panzer_GeometricAggFieldPattern.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class DOFManager
// ************************************************************

DOFManager::DOFManager()
   : numFields_(0)
{ }

DOFManager::DOFManager(Teuchos::RCP<ConnManager<int,int> > & connMngr,MPI_Comm mpiComm)
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
void DOFManager::setConnManager(Teuchos::RCP<ConnManager<int,int> > & connMngr,MPI_Comm mpiComm)
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
Teuchos::RCP<ConnManager<int,int> > DOFManager::resetIndices()
{
   Teuchos::RCP<ConnManager<int,int> > connMngr = connMngr_;

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

void DOFManager::addField(const std::string & str,const Teuchos::RCP<const FieldPattern> & pattern)
{
   for(int i=0;i<connMngr_->numElementBlocks();i++)
      addField(i,str,pattern);
}

void DOFManager::addField(int blockId,const std::string & str,const Teuchos::RCP<const FieldPattern> & pattern)
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

int DOFManager::getFieldNum(const std::string & str) const
{
   std::map<std::string,int>::const_iterator itr = fieldStrToInt_.find(str);

   // return based on what was found
   if(itr==fieldStrToInt_.end())
      return -1;
   else
      return itr->second;
}

int DOFManager::getNumFields() const
{
   return vectorSpace_->getNumFields();
}

// build the global unknown numberings
//   1. this builds the pattens
//   2. initializes the connectivity
//   3. calls initComplete
void DOFManager::buildGlobalUnknowns()
{
   // build the pattern for the ID layout on the mesh
   std::vector<RCP<const FieldPattern> > patVector;
   RCP<GeometricAggFieldPattern> aggFieldPattern = Teuchos::rcp(new GeometricAggFieldPattern);;
   std::map<std::pair<int,int>,Teuchos::RCP<const FieldPattern> >::iterator f2p_itr;
   for(f2p_itr=fieldIntToPattern_.begin();f2p_itr!=fieldIntToPattern_.end();f2p_itr++)
      patVector.push_back(f2p_itr->second);
   aggFieldPattern->buildPattern(patVector);

   // setup connectivity mesh
   connMngr_->buildConnectivity(*aggFieldPattern);
   patternNum_.resize(connMngr_->numElementBlocks()); 
   for(int blockIndex=0;blockIndex<connMngr_->numElementBlocks();blockIndex++) {
      // build the pattern
      buildPattern(blockIndex,aggFieldPattern);

      // figure out what IDs are active for this pattern
      const std::vector<int> & numFieldsPerID = fieldAggPattern_[blockIndex]->numFieldsPerId();
      std::vector<int> activeIds;
      for(std::size_t i=0;i<numFieldsPerID.size();i++)
         if(numFieldsPerID[i]>0) 
            activeIds.push_back(i);
      std::vector<int> reduceConn(activeIds.size()); // which IDs to use
   
      // grab elements for this block
      const std::vector<int> & elements = connMngr_->getElementBlock(blockIndex);

      // build graph for this block
      matrixGraph_->initConnectivityBlock(blockIndex,elements.size(),patternNum_[blockIndex]);
      for(std::size_t e=0;e<elements.size();e++) {
         const int * conn = connMngr_->getConnectivity(elements[e]);
         for(std::size_t i=0;i<activeIds.size();i++)
            reduceConn[i] = conn[activeIds[i]];
 
         matrixGraph_->initConnectivity(blockIndex,elements[e],&reduceConn[0]);
      }
   }
   matrixGraph_->initComplete();
}

// build the pattern associated with this manager
void DOFManager::buildPattern(int blockIndex,const RCP<FieldPattern> & geomPattern)
{
   using Teuchos::rcp;
   using Teuchos::RCP;

   const std::set<int> & fieldSet = blockToField_[blockIndex];
   // std::map<int,Teuchos::RCP<const FieldPattern> > blockPatterns;
   std::vector<std::pair<int,Teuchos::RCP<const FieldPattern> > > blockPatterns;

   // get a map of field patterns
   std::set<int>::const_iterator itr;
   for(itr=fieldSet.begin();itr!=fieldSet.end();++itr) {
      blockPatterns.push_back(std::make_pair(*itr,fieldIntToPattern_[std::make_pair(blockIndex,*itr)])); 
   }

   // smash together all fields...do interlacing
   fieldAggPattern_[blockIndex] = rcp(new FieldAggPattern(blockPatterns));
   TEUCHOS_ASSERT(geomPattern->equals(*fieldAggPattern_[blockIndex]->getGeometricAggFieldPattern()));

   // std::cout << "Agg FEI (block = " << blockIndex << ")" << std::endl;
   // fieldAggPattern_[blockIndex]->print(std::cout);
 
   // build FEI pattern
   const std::vector<int> & fields = fieldAggPattern_[blockIndex]->fieldIds();
   const std::vector<int> & numFieldsPerID = fieldAggPattern_[blockIndex]->numFieldsPerId();

   std::vector<int> reduceNumFieldsPerID;
   for(std::size_t i=0;i<numFieldsPerID.size();i++)
      if(numFieldsPerID[i]>0) 
         reduceNumFieldsPerID.push_back(numFieldsPerID[i]);

   int idsPerSimplex  = reduceNumFieldsPerID.size();

   patternNum_[blockIndex] 
         = matrixGraph_->definePattern(idsPerSimplex,nodeType_,&reduceNumFieldsPerID[0],&fields[0]);
}

// "Get" functions
/////////////////////////////////////////////////////////////////////

// get the map from the matrix
const Teuchos::RCP<Epetra_Map> DOFManager::getMap() const
{
   if(map_==Teuchos::null) map_ = buildMap();

   return map_;
}

const Teuchos::RCP<Epetra_Map> DOFManager::getOverlapMap() const
{
   if(overlappedMap_==Teuchos::null) overlappedMap_ = buildOverlapMap();

   return overlappedMap_;
}

// get the graph of the crs matrix
const Teuchos::RCP<Epetra_CrsGraph> DOFManager::getGraph() const
{
   if(graph_==Teuchos::null) graph_ = buildGraph();

   return graph_;
}

const Teuchos::RCP<Epetra_CrsGraph> DOFManager::getOverlapGraph() const
{
   if(overlappedGraph_==Teuchos::null) overlappedGraph_ = buildOverlapGraph();

   return overlappedGraph_;
}

const Teuchos::RCP<Epetra_Import> DOFManager::getOverlapImport() const
{
   return Teuchos::rcp(new Epetra_Import(*getOverlapMap(),*getMap()));
}

const Teuchos::RCP<Epetra_Export> DOFManager::getOverlapExport() const
{
   return Teuchos::rcp(new Epetra_Export(*getOverlapMap(),*getMap()));
}

// "Build" functions
/////////////////////////////////////////////////////////////////////

const Teuchos::RCP<Epetra_Map> DOFManager::buildMap() const
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
const Teuchos::RCP<Epetra_Map> DOFManager::buildOverlapMap() const
{
   std::vector<int> indices;

   // get the global indices
   vectorSpace_->getIndices_SharedAndOwned(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));
}

// get the graph of the crs matrix
const Teuchos::RCP<Epetra_CrsGraph> DOFManager::buildGraph() const
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

const Teuchos::RCP<Epetra_CrsGraph> DOFManager::buildOverlapGraph() const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> map = getOverlapMap();
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*map,0));

   // graph information about the mesh
   for(int blockIndex=0;blockIndex<connMngr_->numElementBlocks();blockIndex++) {
      // grab elements for this block
      const std::vector<int> & elements = connMngr_->getElementBlock(blockIndex);

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

void DOFManager::getElementGIDs(int localElmtId,std::vector<int> & gids) const
{
   // get information about number of indicies
   int blockId = connMngr_->getBlockId(localElmtId);
   int dof = matrixGraph_->getConnectivityNumIndices(blockId);
   std::vector<int> indices(dof);

   // get elements indices
   int localSize = -1;
   matrixGraph_->getConnectivityIndices(blockId,localElmtId,dof,&indices[0],localSize);

   // copy the indices
   gids = indices;
}

void DOFManager::printFieldInformation(std::ostream & os) const
{
   os << "DOFManager Field Information: " << std::endl;
   
   std::map<int,Teuchos::RCP<FieldAggPattern> >::const_iterator iter;
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

Teuchos::RCP<const FieldPattern> DOFManager::getFieldPattern(int blockId, int fieldNum) const
{
   std::map<std::pair<int,int>,Teuchos::RCP<const FieldPattern> >::const_iterator itr;
   itr = fieldIntToPattern_.find(std::make_pair(blockId,fieldNum));

   if(itr==fieldIntToPattern_.end()) {
      // could not find requiested field pattern...return null
      return Teuchos::null;
   }

   return itr->second;
}

}
