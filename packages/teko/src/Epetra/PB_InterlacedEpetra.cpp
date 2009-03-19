#include "Epetra/PB_InterlacedEpetra.hpp"

#include <vector>

using Teuchos::RCP;
using Teuchos::rcp;

namespace PB {
namespace Epetra {

// this assumes that there are numGlobals with numVars each interlaced
// i.e. for numVars = 2 (u,v) then the vector is
//    [u_0,v_0,u_1,v_1,u_2,v_2, ..., u_(numGlobals-1),v_(numGlobals-1)]
void buildSubMaps(int numGlobals,int numVars,const Epetra_Comm & comm,std::vector<std::pair<int,RCP<Epetra_Map> > > & subMaps)
{
   std::vector<int> vars;
   
   // build vector describing the sub maps
   for(int i=0;i<numVars;i++) vars.push_back(1);

   // build all the submaps
   buildSubMaps(numGlobals,vars,comm,subMaps);
}

// build maps to make other conversions
void buildSubMaps(int numGlobals,const std::vector<int> & vars,const Epetra_Comm & comm,std::vector<std::pair<int,Teuchos::RCP<Epetra_Map> > > & subMaps)
{
   std::vector<int>::const_iterator varItr;

   // compute total number of variables...largely for fun and debugging
   int numGlobalVars = 0;
   for(varItr=vars.begin();varItr!=vars.end();++varItr)
      numGlobalVars += *varItr;

   // must be an even number of globals
   TEUCHOS_ASSERT((numGlobals%numGlobalVars)==0);

   Epetra_Map sampleMap(numGlobals/numGlobalVars,0,comm);
   int numBlocks  = sampleMap.NumMyElements();
   int minBlockID = sampleMap.MinMyGID();

   subMaps.clear();

   // index into local block in strided map
   int blockOffset = 0;
   for(varItr=vars.begin();varItr!=vars.end();++varItr) {
      int numLocalVars = *varItr;
      int numAllElmts = numLocalVars*numGlobals/numGlobalVars;

      // build a sample map for parallel decomposition
      int numMyElmts = numLocalVars * numBlocks;

      // create global arrays describing the as of yet uncreated maps
      std::vector<int> subGlobals;

      // loop over each block of variables
      for(int blockNum=0;blockNum<numBlocks;blockNum++) {

         // loop over each local variable in the block
         for(int local=0;local<numLocalVars;++local) {
            // global block number = minGID+blockNum 
            // block begin global id = numGlobalVars*(minGID+blockNum)
            // global id block offset = blockOffset+local
            subGlobals.push_back((minBlockID+blockNum)*numGlobalVars+blockOffset+local);
         }
      }

      // sanity check
      assert(numMyElmts==subGlobals.size());

      // create an actual map
      RCP<Epetra_Map> subMap = rcp(new Epetra_Map(numAllElmts,numMyElmts,&subGlobals[0],0,comm));
      subMaps.push_back(std::make_pair(numLocalVars,subMap));

      // std::cout << "all = " << numAllElmts << ", mine = " << numMyElmts << std::endl;
      // std::cout << "all:  min = " << subMap->MinAllGID() << ", max = " << subMap->MaxAllGID() << std::endl;
      // std::cout << "mine: min = " << subMap->MinMyGID()  << ", max = " << subMap->MaxMyGID()  << std::endl;

      // update the block offset
      blockOffset += numLocalVars;
   }
}

void buildExportImport(const Epetra_Map & baseMap, const std::vector<std::pair<int,RCP<Epetra_Map> > > & subMaps,
                       std::vector<RCP<Epetra_Export> > & subExport,
                       std::vector<RCP<Epetra_Import> > & subImport)
{
   std::vector<std::pair<int,RCP<Epetra_Map> > >::const_iterator mapItr;

   // build importers and exporters
   for(mapItr=subMaps.begin();mapItr!=subMaps.end();++mapItr) {
      // exctract basic map
      const Epetra_Map & map = *(mapItr->second);

      // add new elements to vectors
      subImport.push_back(rcp(new Epetra_Import(map,baseMap)));
      subExport.push_back(rcp(new Epetra_Export(map,baseMap)));
   }
}

void buildSubVectors(const std::vector<std::pair<int,RCP<Epetra_Map> > > & subMaps,std::vector<RCP<Epetra_MultiVector> > & subVectors,int count)
{
   std::vector<std::pair<int,RCP<Epetra_Map> > >::const_iterator mapItr;

   // build vectors, importers and exporters
   for(mapItr=subMaps.begin();mapItr!=subMaps.end();++mapItr) {
      // exctract basic map
      const Epetra_Map & map = *(mapItr->second);

      // add new elements to vectors
      // subVectors.push_back(rcp(new Epetra_Vector(map)));
      subVectors.push_back(rcp(new Epetra_MultiVector(map,count)));
   }
}

// build a single subblock Epetra_CrsMatrix
RCP<Epetra_CrsMatrix> buildSubBlock(int i,int j,const Epetra_CrsMatrix & A,const std::vector<std::pair<int,RCP<Epetra_Map> > > & subMaps)
{
   // get the number of variables families
   int numVarFamily = subMaps.size();

   TEUCHOS_ASSERT(i>=0 && i<numVarFamily);
   TEUCHOS_ASSERT(j>=0 && j<numVarFamily);

   const Epetra_Map & rowMap = *subMaps[i].second;
   const Epetra_Map & colMap = *subMaps[j].second;
   int rowFamilyCnt = subMaps[i].first;
   int colFamilyCnt = subMaps[j].first;

   // compute the number of global variables
   // and the row and column block offset
   int numGlobalVars = 0;
   int rowBlockOffset = 0;
   int colBlockOffset = 0;
   for(int k=0;k<numVarFamily;k++) {
      numGlobalVars += subMaps[k].first;
 
      // compute block offsets
      if(k<i) rowBlockOffset += subMaps[k].first;
      if(k<j) colBlockOffset += subMaps[k].first;
   }

   // get the number of nodal blocks in the matrix
   int numBlocks = A.NumGlobalRows()/numGlobalVars;

   // copy all global rows to here
   Epetra_Import import(rowMap,A.RowMap());
   Epetra_CrsMatrix localA(Copy,rowMap,0);
   localA.Import(A,import,Insert);

   RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy,rowMap,0));

   // get entry information
   int numMyRows = rowMap.NumMyElements();
   int maxNumEntries = A.GlobalMaxNumEntries();

   // for extraction
   int indicies[maxNumEntries];
   double values[maxNumEntries];

   // for insertion
   int colIndicies[maxNumEntries];
   double colValues[maxNumEntries];

   // std::cout << "building (i,j) = ( " << i << ", " << j << " )" << std::endl;

   // insert each row into subblock
   // let FillComplete handle column distribution
   for(int localRow=0;localRow<numMyRows;localRow++) {
      int numEntries = -1; 
      int globalRow = rowMap.GID(localRow);

      TEUCHOS_ASSERT(globalRow>=0);

      // extract a global row copy
      int err = localA.ExtractGlobalRowCopy(globalRow, maxNumEntries, numEntries, values, indicies);
      TEUCHOS_ASSERT(err==0);

      int numOwnedCols = 0;
      for(int localCol=0;localCol<numEntries;localCol++) {
         int globalCol = indicies[localCol];

         // determinate which block this column ID is in
         int block = globalCol / numGlobalVars;
         
         bool inFamily = true; 
 
         // test the beginning of the block
         inFamily &= (block*numGlobalVars+colBlockOffset <= globalCol);
         inFamily &= ((block*numGlobalVars+colBlockOffset+colFamilyCnt) > globalCol);

         // is this column in the variable family
         // if(globalCol % numVars-j==0) {
         if(inFamily) {
            colIndicies[numOwnedCols] = indicies[localCol];
            colValues[numOwnedCols] = values[localCol];

            numOwnedCols++;

       //      std::cout << globalCol << " ";
         }
      }

      // insert it into the new matrix
      mat->InsertGlobalValues(globalRow,numOwnedCols,colValues,colIndicies);
   }
   // std::cout << std::endl;

   // fill it and automagically optimize the storage
   mat->FillComplete(colMap,rowMap);

   return mat;
}


// collect subvectors into a single global vector
void many2one(Epetra_MultiVector & one, const std::vector<RCP<const Epetra_MultiVector> > & many,
                                   const std::vector<RCP<Epetra_Export> > & subExport)
{
   // std::vector<RCP<const Epetra_Vector> >::const_iterator vecItr;
   std::vector<RCP<const Epetra_MultiVector> >::const_iterator vecItr;
   std::vector<RCP<Epetra_Export> >::const_iterator expItr;

   // using Exporters fill the empty vector from the sub-vectors
   for(vecItr=many.begin(),expItr=subExport.begin();
       vecItr!=many.end();++vecItr,++expItr) {
      one.Export(**vecItr,**expItr,Insert);
   }
}

// distribute one global vector into a many subvectors
void one2many(std::vector<RCP<Epetra_MultiVector> > & many,const Epetra_MultiVector & single,
                                   const std::vector<RCP<Epetra_Import> > & subImport)
{
   // std::vector<RCP<Epetra_Vector> >::const_iterator vecItr;
   std::vector<RCP<Epetra_MultiVector> >::const_iterator vecItr;
   std::vector<RCP<Epetra_Import> >::const_iterator impItr;

   // using Importers fill the sub vectors from the mama vector
   for(vecItr=many.begin(),impItr=subImport.begin();
       vecItr!=many.end();++vecItr,++impItr) {
      (*vecItr)->Import(single,**impItr,Insert);
   }
}

} // end namespace Epetra 
} // end namespace PB
