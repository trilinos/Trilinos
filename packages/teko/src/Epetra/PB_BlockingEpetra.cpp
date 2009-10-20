#include "PB_BlockingEpetra.hpp"

using Teuchos::RCP;
using Teuchos::rcp;

namespace PB {
namespace Epetra {
namespace Blocking {

/** Build maps to make other conversions. This function builds a map 
  * using a vector of global ids local to this processor.  It also builds
  * a seperate map that (globally) starts from zero. For instance if the following
  * GIDs are passed in PID = 0, GID = [0 2 4 6 8] and PID = 1, GID = [10 12 14]
  * the the two maps created are
  *    Global Map = [(PID=0,GID=[0 2 4 6 8]), (PID=1,GID=[10 12 14])]
  *    Contiguous Map = [(PID=0,GID=[0 1 2 3 4]), (PID=1,GID=[5 6 7])]
  *
  * \param[in] gid Local global IDs to use
  * \param[in] comm Communicator to use in construction of the maps
  *
  * \returns A pair of maps: (Global Map, Contiguous Map)
  */
const MapPair buildSubMap(const std::vector< int > & gid, const Epetra_Comm &comm)
{
   Teuchos::RCP<Epetra_Map> gidMap = rcp(new Epetra_Map(-1,gid.size(),&gid[0],0,comm));
   Teuchos::RCP<Epetra_Map> contigMap = rcp(new Epetra_Map(-1,gid.size(),0,comm));

   return std::make_pair(gidMap,contigMap); 
}

/** Build the Export/Import objects that take the single vector global map and 
  * build individual sub maps.
  *
  * \param[in] baseMap Single global vector map.
  * \param[in] maps Pair of maps containing the global sub map and the contiguous sub map.
  *                 These come directly from <code>buildSubMap</code>.
  *
  * \returns A pair containing pointers to the Import/Export objects.
  */
const ImExPair buildExportImport(const Epetra_Map & baseMap,const MapPair & maps)
{
   return std::make_pair(rcp(new Epetra_Import(*maps.first,baseMap)),
                         rcp(new Epetra_Export(*maps.first,baseMap)));
}

/** Using a list of map pairs created by <code>buildSubMap</code>, buidl the corresponding
  * multi-vector objects.
  *
  * \param[in] maps Map pairs created by <code>buildSubMap</code>
  * \param[in,out] vectors Vector objects created using the Contiguous maps
  * \param[in] count Number of multivectors to build.
  */
void buildSubVectors(std::vector<MapPair> & maps,
                     std::vector<RCP<Epetra_MultiVector> > & vectors,int count)
{
   std::vector<MapPair>::const_iterator mapItr;
   
   // loop over all maps
   for(mapItr=maps.begin();mapItr!=maps.end();mapItr++) {
      // add new elements to vectors
      RCP<Epetra_MultiVector> mv = rcp(new Epetra_MultiVector(*(*mapItr).second,count));
      vectors.push_back(mv);
   } 
}

/** Copy the contents of a global vector into many sub-vectors created by <code>buildSubVectors</code>.
  *
  * \param[in,out] many Sub-vector to be filled by this operation created by <code>buildSubVectors</code>.
  * \param[in] one The source vector.
  * \param[in] subImport A list of import objects to use in copying.
  */
void one2many(std::vector<RCP<Epetra_MultiVector> > & many, const Epetra_MultiVector & single,
                                                            const std::vector<RCP<Epetra_Import> > & subImport)
{
   // std::vector<RCP<Epetra_Vector> >::const_iterator vecItr;
   std::vector<RCP<Epetra_MultiVector> >::const_iterator vecItr;
   std::vector<RCP<Epetra_Import> >::const_iterator impItr;

   // using Importers fill the sub vectors from the mama vector
   for(vecItr=many.begin(),impItr=subImport.begin();
       vecItr!=many.end();++vecItr,++impItr) {
      // for ease of access to the destination
      RCP<Epetra_MultiVector> destVec = *vecItr;

      // extract the map with global indicies from the current vector
      const Epetra_BlockMap & globalMap = (*impItr)->TargetMap();

      // build the import vector as a view on the destination
      Epetra_MultiVector importVector(View,globalMap,destVec->Values(),destVec->Stride(),destVec->NumVectors());

      // perform the import
      importVector.Import(single,**impItr,Insert);
   }
}

/** Copy the contents of many sub vectors (created from a contigous sub maps) to a single global
  * vector. This should have the map used to create the Export/Import objects in the <code>buildExportImport</code>
  * function. If more then one sub vector contains values for a particular GID in the single vector
  * then the value in the final vector will be that of the last sub vector (in the list <code>many</code>).
  *
  * \param[in,out] one The single vector to be filled by this operation.
  * \param[in] many Sub-vectors created by <code>buildSubVectors</code> used to fill <code>one</code>.
  * \param[in] subExport A list of export objects to use in copying.
  */
void many2one(Epetra_MultiVector & one, const std::vector<RCP<const Epetra_MultiVector> > & many,
                                        const std::vector<RCP<Epetra_Export> > & subExport)
{
   // std::vector<RCP<const Epetra_Vector> >::const_iterator vecItr;
   std::vector<RCP<const Epetra_MultiVector> >::const_iterator vecItr;
   std::vector<RCP<Epetra_Export> >::const_iterator expItr;

   // using Exporters fill the empty vector from the sub-vectors
   for(vecItr=many.begin(),expItr=subExport.begin();
       vecItr!=many.end();++vecItr,++expItr) {

      // for ease of access to the source
      RCP<const Epetra_MultiVector> srcVec = *vecItr;

      // extract the map with global indicies from the current vector
      const Epetra_BlockMap & globalMap = (*expItr)->SourceMap();

      // build the export vector as a view of the destination
      Epetra_MultiVector exportVector(View,globalMap,srcVec->Values(),srcVec->Stride(),srcVec->NumVectors());
      one.Export(exportVector,**expItr,Insert);
   }
}

// build a single subblock Epetra_CrsMatrix
RCP<Epetra_CrsMatrix> buildSubBlock(int i,int j,const Epetra_CrsMatrix & A,const std::vector<MapPair> & subMaps)
{
   // get the number of variables families
   int numVarFamily = subMaps.size();

   TEUCHOS_ASSERT(i>=0 && i<numVarFamily);
   TEUCHOS_ASSERT(j>=0 && j<numVarFamily);

   const Epetra_Map & gRowMap = *subMaps[i].first;
   const Epetra_Map & gColMap = *subMaps[j].first;
   const Epetra_Map & rowMap = *subMaps[i].second;
   const Epetra_Map & colMap = *subMaps[j].second;

   RCP<Epetra_CrsMatrix> mat = rcp(new Epetra_CrsMatrix(Copy,rowMap,0));

   // get entry information
   int numMyRows = rowMap.NumMyElements();
   int maxNumEntries = A.GlobalMaxNumEntries();

   // for extraction
   std::vector<int> indices(maxNumEntries);
   std::vector<double> values(maxNumEntries);

   // for insertion
   std::vector<int> colIndices(maxNumEntries);
   std::vector<double> colValues(maxNumEntries);

   // insert each row into subblock
   // let FillComplete handle column distribution
   for(int localRow=0;localRow<numMyRows;localRow++) {
      int numEntries = -1;
      int globalRow = gRowMap.GID(localRow);
      int contigRow = rowMap.GID(localRow);

      // extract a global row copy
      int err = A.ExtractGlobalRowCopy(globalRow, maxNumEntries, numEntries, &values[0], &indices[0]);
      TEUCHOS_ASSERT(err==0);

      int numOwnedCols = 0;
      for(int localCol=0;localCol<numEntries;localCol++) {
         // if global id is not owned by this column
         int gLID = gColMap.LID(indices[localCol]);
         if(gLID==-1) continue;

         // get the contiguous column ID
         int contigCol = colMap.GID(gLID);
         TEUCHOS_ASSERT(contigCol>=0);

         colIndices[numOwnedCols] = contigCol;
         colValues[numOwnedCols] = values[localCol];
         numOwnedCols++;
      }

      // insert it into the new matrix
      mat->InsertGlobalValues(contigRow,numOwnedCols,&colValues[0],&colIndices[0]);
   }

   // fill it and automagically optimize the storage
   mat->FillComplete(colMap,rowMap);

   return mat;
}

// build a single subblock Epetra_CrsMatrix
void rebuildSubBlock(int i,int j,const Epetra_CrsMatrix & A,const std::vector<MapPair> & subMaps,Epetra_CrsMatrix & mat)
{
   // get the number of variables families
   int numVarFamily = subMaps.size();

   TEUCHOS_ASSERT(i>=0 && i<numVarFamily);
   TEUCHOS_ASSERT(j>=0 && j<numVarFamily);
   TEUCHOS_ASSERT(mat.Filled());

   const Epetra_Map & gRowMap = *subMaps[i].first;
   const Epetra_Map & gColMap = *subMaps[j].first;
   const Epetra_Map & rowMap = *subMaps[i].second;
   const Epetra_Map & colMap = *subMaps[j].second;

   mat.PutScalar(0.0);

   // get entry information
   int numMyRows = rowMap.NumMyElements();
   int maxNumEntries = A.GlobalMaxNumEntries();

   // for extraction
   std::vector<int> indices(maxNumEntries);
   std::vector<double> values(maxNumEntries);

   // for insertion
   std::vector<int> colIndices(maxNumEntries);
   std::vector<double> colValues(maxNumEntries);

   // insert each row into subblock
   // let FillComplete handle column distribution
   for(int localRow=0;localRow<numMyRows;localRow++) {
      int numEntries = -1;
      int globalRow = gRowMap.GID(localRow);
      int contigRow = rowMap.GID(localRow);

      // extract a global row copy
      int err = A.ExtractGlobalRowCopy(globalRow, maxNumEntries, numEntries, &values[0], &indices[0]);
      TEUCHOS_ASSERT(err==0);

      int numOwnedCols = 0;
      for(int localCol=0;localCol<numEntries;localCol++) {
         // if global id is not owned by this column
         int gLID = gColMap.LID(indices[localCol]);
         if(gLID==-1) continue;

         // get the contiguous column ID
         int contigCol = colMap.GID(gLID);
         TEUCHOS_ASSERT(contigCol>=0);

         colIndices[numOwnedCols] = contigCol;
         colValues[numOwnedCols] = values[localCol];
         numOwnedCols++;
      }

      // insert it into the new matrix
      mat.SumIntoGlobalValues(contigRow,numOwnedCols,&colValues[0],&colIndices[0]);
   }
}

} // end Blocking
} // end Epetra
} // end PB
