//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER

#ifndef EPETRAEXT_MMHELPERS_H
#define EPETRAEXT_MMHELPERS_H

#include "EpetraExt_ConfigDefs.h"
#include "Epetra_DistObject.h"
#include "Epetra_Map.h"
#include "Teuchos_RCP.hpp"

#include <vector>
#include <set>
#include <map>


class Epetra_CrsMatrix;
class Epetra_Import;
class Epetra_Distributor;

namespace EpetraExt {
class LightweightCrsMatrix;

// Turn on to enable the special MMM timers
//#define ENABLE_MMM_TIMINGS

//#define HAVE_EPETRAEXT_DEBUG // for extra sanity checks

// ==============================================================
//struct that holds views of the contents of a CrsMatrix. These
//contents may be a mixture of local and remote rows of the
//actual matrix.
class CrsMatrixStruct {
public:
  CrsMatrixStruct();

  virtual ~CrsMatrixStruct();

  void deleteContents();

  int numRows;

  // The following class members get used in the transpose modes of the MMM
  // but not in the A*B mode.  
  int* numEntriesPerRow;
  int** indices;
  double** values;
  bool* remote;
  int numRemote;
  const Epetra_BlockMap* importColMap;

  // Maps and matrices (all modes)
  const Epetra_Map* origRowMap;
  const Epetra_Map* rowMap;
  const Epetra_Map* colMap;
  const Epetra_Map* domainMap;
  LightweightCrsMatrix* importMatrix;
  const Epetra_CrsMatrix *origMatrix;

  // The following class members are only used for A*B mode
  std::vector<int> targetMapToOrigRow;
  std::vector<int> targetMapToImportRow;
};

int dumpCrsMatrixStruct(const CrsMatrixStruct& M);

// ==============================================================
class CrsWrapper {
 public:
  virtual ~CrsWrapper(){}

  virtual const Epetra_Map& RowMap() const = 0;

  virtual bool Filled() = 0;

  virtual int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices) = 0;

  virtual int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices) = 0;
};

// ==============================================================
class CrsWrapper_Epetra_CrsMatrix : public CrsWrapper {
 public:
  CrsWrapper_Epetra_CrsMatrix(Epetra_CrsMatrix& epetracrsmatrix);
  virtual ~CrsWrapper_Epetra_CrsMatrix();

  const Epetra_Map& RowMap() const;

  bool Filled();

  int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);

 private:
  Epetra_CrsMatrix& ecrsmat_;
};

// ==============================================================
class CrsWrapper_GraphBuilder : public CrsWrapper {
 public:
  CrsWrapper_GraphBuilder(const Epetra_Map& emap);
  virtual ~CrsWrapper_GraphBuilder();

  const Epetra_Map& RowMap() const {return rowmap_; }

  bool Filled();

  int InsertGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);
  int SumIntoGlobalValues(int GlobalRow, int NumEntries, double* Values, int* Indices);

  std::map<int,std::set<int>*>& get_graph();

  int get_max_row_length() { return max_row_length_; }

 private:
  std::map<int,std::set<int>*> graph_;
  const Epetra_Map& rowmap_;
  int max_row_length_;
};

// ==============================================================
void insert_matrix_locations(CrsWrapper_GraphBuilder& graphbuilder,
                              Epetra_CrsMatrix& C);

void pack_outgoing_rows(const Epetra_CrsMatrix& mtx,
                        const std::vector<int>& proc_col_ranges,
                        std::vector<int>& send_rows,
                        std::vector<int>& rows_per_send_proc);

std::pair<int,int> get_col_range(const Epetra_Map& emap);

std::pair<int,int> get_col_range(const Epetra_CrsMatrix& mtx);

// ==============================================================
class LightweightMapData : Epetra_Data {
  friend class LightweightMap; 
 public:
  LightweightMapData();
  ~LightweightMapData();
  int IndexBase_;
  std::vector<int> MyGlobalElements_; 
  Epetra_HashTable<int> * LIDHash_;

  // For "copy" constructor only...
  Epetra_Map * CopyMap_;
  
 };

class LightweightMap {
 public:
  LightweightMap();
  LightweightMap(int NumGlobalElements,int NumMyElements, const int * MyGlobalElements, int IndexBase, bool GenerateHash=true);
  LightweightMap(const Epetra_Map & Map);
  LightweightMap(const LightweightMap & Map);
  ~LightweightMap();

  LightweightMap & operator=(const LightweightMap & map);
  int LID(int GID) const;
  int GID(int LID) const;      
  int NumMyElements() const;
  int* MyGlobalElements() const; 
  int IndexBase() const {return Data_->IndexBase_;}

  int MinLID() const;
  int MaxLID() const;

 private:
  void CleanupData();
  LightweightMapData *Data_;
  //Epetra_BlockMapData* Data_;
};


// ==============================================================
class RemoteOnlyImport {
 public:
  RemoteOnlyImport(const Epetra_Import & Importer, LightweightMap & RemoteOnlyTargetMap);
  ~RemoteOnlyImport();

  int NumSameIDs() {return 0;}

  int NumPermuteIDs() {return 0;} 

  int NumRemoteIDs() {return NumRemoteIDs_;}

  int NumExportIDs() {return NumExportIDs_;}

  int* ExportLIDs() {return ExportLIDs_;}

  int* ExportPIDs() {return ExportPIDs_;}

  int* RemoteLIDs() {return RemoteLIDs_;}  

  int* PermuteToLIDs() {return 0;}

  int* PermuteFromLIDs() {return 0;}

  int NumSend() {return NumSend_;}
  
  Epetra_Distributor & Distributor() {return *Distor_;}  

  const Epetra_BlockMap & SourceMap() const {return *SourceMap_;}
  const LightweightMap & TargetMap() const {return *TargetMap_;}

 private:
  int NumSend_;
  int NumRemoteIDs_;
  int NumExportIDs_;
  int * ExportLIDs_;
  int * ExportPIDs_;
  int * RemoteLIDs_;
  Epetra_Distributor* Distor_;
  const Epetra_BlockMap* SourceMap_;
  const LightweightMap  *TargetMap_;
};
   
// ==============================================================
class LightweightCrsMatrix {
 public:
  LightweightCrsMatrix(const Epetra_CrsMatrix & A, RemoteOnlyImport & RowImporter);
  LightweightCrsMatrix(const Epetra_CrsMatrix & A, Epetra_Import & RowImporter);
  ~LightweightCrsMatrix();

  // Standard crs data structures
  std::vector<int>    rowptr_;
  std::vector<int>    colind_;
  std::vector<double> vals_;

  // Colind in LL-GID space (if needed)
  std::vector<long long>   colind_LL_;

  // Epetra Maps
  bool                     use_lw;
  LightweightMap           *RowMapLW_;
  Epetra_BlockMap          *RowMapEP_;
  LightweightMap           ColMap_;
  Epetra_Map               DomainMap_;


  // List of owning PIDs (from the DomainMap) as ordered by entries in the column map.
  std::vector<int>    ColMapOwningPIDs_;

  // For building the final importer for C
  std::vector<int>    ExportLIDs_;
  std::vector<int>    ExportPIDs_;

 private: 

  template <typename ImportType>
  void Construct(const Epetra_CrsMatrix & A, ImportType & RowImporter);

  // Templated versions of MakeColMapAndReindex (to prevent code duplication)
  template <class GO>
  int MakeColMapAndReindex(std::vector<int> owningPIDs,std::vector<GO> Gcolind);

  template <int>
  int MakeColMapAndReindex(std::vector<int> owningPIDs,std::vector<int> Gcolind);

  template <long long>
  int MakeColMapAndReindex(std::vector<int> owningPIDs,std::vector<long long> Gcolind);


  template<typename ImportType>
  int PackAndPrepareReverseComm(const Epetra_CrsMatrix & SourceMatrix, ImportType & RowImporter,
				std::vector<int> &ReverseSendSizes, std::vector<int> &ReverseSendBuffer);

  template<typename ImportType>
  int MakeExportLists(const Epetra_CrsMatrix & SourceMatrix, ImportType & RowImporter,
		      std::vector<int> &ReverseRecvSizes, const int *ReverseRecvBuffer,
		      std::vector<int> & ExportPIDs, std::vector<int> & ExportLIDs);

};

}//namespace EpetraExt

#endif

