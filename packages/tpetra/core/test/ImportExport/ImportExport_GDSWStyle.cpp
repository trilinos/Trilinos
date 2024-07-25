// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Core.hpp>
#include <Teuchos_LocalTestingHelpers.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include <iterator>
#include <numeric>
#include <sstream>
#include <ostream>
#include <fstream>
#include <unistd.h>
#include <sys/resource.h>

#include <Tpetra_Distributor.hpp>
#include "Teuchos_FancyOStream.hpp"
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Tpetra_Import_Util.hpp"
#include "Tpetra_TestingXMLUtilities.hpp"
#include "GDSW_Proxy.hpp"



namespace {

  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::outArg;
  using Tpetra::global_size_t;
  using std::vector;
  using std::sort;
  using Teuchos::arrayViewFromVector;
  using Teuchos::broadcast;
  using Teuchos::OrdinalTraits;
  using Teuchos::tuple;
  using Teuchos::Range1D;
  using Tpetra::Map;
  using Tpetra::Import;
  using Tpetra::Export;
  using Teuchos::ScalarTraits;
  using Teuchos::Comm;
  using Teuchos::Array;
  using Tpetra::REPLACE;
  using Tpetra::ADD;
  using std::ostream_iterator;
  using std::endl;

  using Tpetra::createContigMapWithNode;



namespace { //anonymous


template<class ScalarType, class LocalOrdinalType, class GlobalOrdinalType, class NodeType>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarType, LocalOrdinalType, GlobalOrdinalType, NodeType> >
createLaplace1D (const Teuchos::RCP<const Tpetra::Map<LocalOrdinalType, GlobalOrdinalType, NodeType> >& rowMap,
                    Teuchos::FancyOStream& out,
                    const bool dbg)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::as;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcpFromRef;
  typedef ScalarType ST;
  typedef LocalOrdinalType LO;
  typedef GlobalOrdinalType GO;
  typedef NodeType NT;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> matrix_type;
  typedef typename matrix_type::local_matrix_device_type LMT;

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm ();
  //#define OLD_AND_MEMORY_HOGGING
#ifndef OLD_AND_MEMORY_HOGGING
  ST ONE = Teuchos::ScalarTraits<ST>::one();

  ST TWO = ONE+ONE;

  // Generate column map
  int rank = comm->getRank();
  int numProcs = comm->getSize();
  GO GO_INVALID = Teuchos::OrdinalTraits<GO>::invalid();

  size_t numRows = rowMap->getLocalNumElements();
  LO nnz = 3*numRows;


  // Make_the column map
  int size = numRows;
  if(rank !=0) size++;
  if(rank !=numProcs-1) size++; 
  Teuchos::Array<GO> col_ids(size);
  for(LO i=0; i<(LO)numRows; i++)
    col_ids[i] = rowMap->getGlobalElement(i);
  LO ct = numRows;

  LO initial_row_length;
  LO final_row_length;
  if(rank!=0) {
    col_ids[ct] = rowMap->getGlobalElement(0) -1;
    initial_row_length = 3;
    ct++;
  }
  else {
    initial_row_length = 2;
    nnz--;
  }

  if(rank!=numProcs-1) {
    col_ids[ct] = rowMap->getGlobalElement(numRows-1) +1;
    final_row_length = 3;
    ct++;
  }
  else {
    final_row_length = 2;
    nnz--;    
  }

  RCP<const map_type> colMap = rcp(new map_type(GO_INVALID,col_ids(),rowMap->getIndexBase(),comm));
  size_t numCols = colMap->getLocalNumElements();


  // Fill the matrix
  typename LMT::values_type::non_const_type  values("values",nnz);
  typename LMT::index_type::non_const_type   colind("colind",nnz);
  typename LMT::row_map_type::non_const_type rowptr("rowptr",numRows+1);

  Kokkos::parallel_for("matrix fill", numRows,KOKKOS_LAMBDA(const LO& row) {     
      if(row == 0) {
        // First row on proc
        LO row_start = 0;
        LO row_stop = initial_row_length;
        rowptr[0] = row_start;
        rowptr[1] = row_stop;
        if(initial_row_length == 2) {
          colind[row_start  ] = row;
          colind[row_start+1] = row+1;
          values[row_start  ] =  TWO;
          values[row_start+1] = -ONE;
        }
        else {
          colind[row_start  ] = numRows;
          colind[row_start+1] = row;
          colind[row_start+2] = row+1;
          values[row_start  ] = -ONE;
          values[row_start+1] =  TWO;
          values[row_start+2] = -ONE;        
        }
      }
      else if (row == (LO)numRows -1) {
        // Last row on proc
        LO row_start = (row-1)*3 + initial_row_length;
        LO row_stop  = row_start+final_row_length;
        rowptr[row+1] = row_stop;
        if(final_row_length == 2) {
          colind[row_start  ] = row-1;
          colind[row_start+1] = row;
          values[row_start  ] = -ONE;
          values[row_start+1] =  TWO;

        }
        else {
          colind[row_start  ] = row-1;
          colind[row_start+1] = row;
          colind[row_start+2] = numRows;
          values[row_start  ] = -ONE;
          values[row_start+1] =  TWO;
          values[row_start+2] = -ONE;        
        }
      }
      else {
        // All other rows
        LO row_start = (row-1)*3 + initial_row_length;
        LO row_stop  = row_start+3;
        rowptr[row+1] = row_stop;
        colind[row_start  ] = row-1;
        colind[row_start+1] = row;
        colind[row_start+2] = row+1;
        values[row_start  ] = -ONE;
        values[row_start+1] =  TWO;
        values[row_start+2] = -ONE;        
      }    
    });

  // Put matrix together
  LMT A_lcl("local",numRows,numCols,nnz,values,rowptr,colind);

  RCP<matrix_type> A = rcp(new matrix_type(A_lcl,rowMap,colMap));

#else
  const size_t myNumElts = rowMap->getLocalNumElements ();

  RCP<matrix_type> A = rcp(new matrix_type(rowMap, 3));

  Array<GO> ind (3);
  Array<ST> val (3);
  for (size_t myRow = 0; myRow < myNumElts; ++myRow) {
    const GO globalRow = rowMap->getGlobalElement (myRow);
    if (globalRow == rowMap->getMinAllGlobalIndex ()) {
      val[0] = as<ST> (2);
      val[1] = as<ST> (-1);
      ind[0] = globalRow;
      ind[1] = globalRow + 1;
      A->insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
    }
    else if (globalRow == rowMap->getMaxAllGlobalIndex ()) {
      val[0] = as<ST> (-1);
      val[1] = as<ST> (2);
      ind[0] = globalRow - 1;
      ind[1] = globalRow;
      A->insertGlobalValues (globalRow, ind.view (0, 2), val.view (0, 2));
    }
    else {
      val[0] = as<ST> (-1);
      val[1] = as<ST> (2);
      val[2] = as<ST> (-1);
      ind[0] = globalRow - 1;
      ind[1] = globalRow;
      ind[2] = globalRow + 1;
      A->insertGlobalValues (globalRow, ind.view (0, 3), val.view (0, 3));
    }
  }
  
  A->fillComplete (rowMap, rowMap);
#endif
  return A;
}


// Input matrices must be fill complete, and all four of their Maps
// (row, column, domain, and range) must be the same.
template<class CrsMatrixType>
bool
compareCrsMatrix (const CrsMatrixType& A_orig, const CrsMatrixType& A)
{
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::RCP;
  using Teuchos::reduceAll;
  using Teuchos::REDUCE_MIN;
  typedef typename CrsMatrixType::global_ordinal_type GO;
  typedef typename ArrayView<const GO>::size_type size_type;
  typedef typename CrsMatrixType::nonconst_global_inds_host_view_type gids_type;
  typedef typename CrsMatrixType::nonconst_values_host_view_type vals_type;

  int localEqual = 1;

  //
  // Are my local matrices equal?
  //
  gids_type indOrig, ind;
  vals_type valOrig, val;
  size_t numEntriesOrig = 0;
  size_t numEntries = 0;

  ArrayView<const GO> localElts = A.getRowMap ()->getLocalElementList ();
  const size_type numLocalElts = localElts.size ();
  for (size_type i = 0; i < numLocalElts; ++i) {
    const GO globalRow = localElts[i];
    numEntriesOrig = A_orig.getNumEntriesInGlobalRow (globalRow);
    numEntries = A.getNumEntriesInGlobalRow (globalRow);

    if (numEntriesOrig != numEntries) {
      localEqual = 0;
      break;
    }
    Kokkos::resize(indOrig,numEntriesOrig);
    Kokkos::resize(valOrig,numEntriesOrig);
    A_orig.getGlobalRowCopy (globalRow, indOrig, valOrig, numEntriesOrig);
    Kokkos::resize(ind,numEntries);
    Kokkos::resize(val,numEntries);
    A.getGlobalRowCopy (globalRow, ind, val, numEntries);

    // Global row entries are not necessarily sorted.  Sort them so
    Tpetra::sort2 (indOrig, indOrig.extent(0), valOrig);
    Tpetra::sort2 (ind, ind.extent(0), val);

    for (size_t k = 0; k < numEntries; ++k) {
      // Values should be _exactly_ equal.
      if (indOrig[k] != ind[k] || valOrig[k] != val[k]) {
        localEqual = 0;
        break;
      }
    }
  }

  RCP<const Comm<int> > comm = A.getRowMap ()->getComm ();
  int globalEqual = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, 1, &localEqual, &globalEqual);
  return globalEqual == 1;
}



// return current memory usage in kilobytes
size_t get_memory_usage_now()
{
  size_t memory = 0; 


#ifdef PROC_STAT
  unsigned long rss_pages = 0;
  // Read memory size data from /proc/pid/stat
  // see "man proc" for details.
  std::ifstream proc_stat("/proc/self/stat");
  if (proc_stat)
  {
    std::string buf;
    proc_stat
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf >> buf >> buf
       >> buf >> buf >> buf
       >> rss_pages;
    proc_stat.close();
  }
  memory = rss_pages * sysconf( _SC_PAGESIZE);

#else  
  // darwin reports rusage.ru_maxrss in bytes
#if defined(__APPLE__) || defined(__MACH__)
#define RU_MAXRSS_UNITS 1024
#else
#define RU_MAXRSS_UNITS 1
#endif

  struct rusage sys_resources;
  getrusage(RUSAGE_SELF, &sys_resources);
  memory = (unsigned long)sys_resources.ru_maxrss / RU_MAXRSS_UNITS;
#endif  


  /* Success */
  return memory;
}


template <class crs_matrix_type, class map_type>
RCP<crs_matrix_type> Filter(const RCP<crs_matrix_type> & A,const RCP<const map_type> &filterMap){
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using indices_type = typename crs_matrix_type::local_inds_host_view_type;
  using values_type = typename crs_matrix_type::values_host_view_type;
  using LO = typename crs_matrix_type::local_ordinal_type;
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  auto rowMap = A->getRowMap();
  auto colMap = A->getColMap();
  size_t nrows = rowMap->getLocalNumElements();

  // Count
  Teuchos::Array<size_t> count(filterMap->getLocalNumElements(),0);
  for(size_t i=0; i<nrows; i++) {
    indices_type  indices;
    values_type   values;
    A->getLocalRowView(i,indices,values);
    LO frow = filterMap->getLocalElement(rowMap->getGlobalElement(i));
    for(LO j=0; j<(LO)indices.size(); j++) {
      if (filterMap->getLocalElement(colMap->getGlobalElement(indices[j])) != LO_INVALID)
        count[frow]++;
    }
  }


  // Alloc
  RCP<crs_matrix_type> B = rcp(new crs_matrix_type(filterMap,filterMap,count()));

  // Fill
  for(size_t i=0; i<nrows; i++) {
    indices_type  indices;
    values_type   values;
    A->getLocalRowView(i,indices,values);
    LO frow = filterMap->getLocalElement(rowMap->getGlobalElement(i));
    for(LO j=0; j<(LO)indices.size(); j++) {
      LO col = filterMap->getLocalElement(colMap->getGlobalElement(indices[j]));
      if(col != LO_INVALID) {
        B->insertLocalValues(frow,1,&values[j],&col);
      }
    }   
  }
  B->fillComplete(rowMap,rowMap);
  
  return B;
}


}//anonymous namespace


  template<class SC, class LO, class GO, class NT>
  void GDSWStyle_Test0(int run_case, bool watchr_output, bool & success ) {
    Teuchos::FancyOStream  out(Teuchos::rcpFromRef(std::cout));

    using map_type = Map<LO, GO, NT>;
    using import_type = Tpetra::Import<LO, GO, NT>;
    using crs_matrix_type = Tpetra::CrsMatrix<SC,LO,GO,NT>;
    bool debug=true;

    const Tpetra::global_size_t INVALID =  Teuchos::OrdinalTraits<GO>::invalid ();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    TpetraFunctions<SC,LO,GO,NT> tFunctions;

    // diagnostic setup
    Teuchos::Array<std::string> names;  names.reserve(50);
    Teuchos::Array<double> memory;      memory.reserve(50);

    // create map
    RCP<const map_type > rowMap1to1 =  createContigMapWithNode<LO, GO, NT> (INVALID, 1000000, comm);
    size_t mem0 = get_memory_usage_now();

    // Create matrix
    RCP<crs_matrix_type> A = createLaplace1D<SC,LO,GO,NT>(rowMap1to1,out,debug);

    // Grab the column map to create a proxy for the input region
    RCP<const map_type> regionMap = A->getColMap();
    size_t mem1 = get_memory_usage_now();
    RCP<crs_matrix_type> outputMatrix1,outputMatrix2,outputMatrix3,outputMatrix4,outputMatrix5,outputMatrix6,outputMatrix7;

    // 1) Tpetra-based code duplication what is in GDSW
    import_type Importer(rowMap1to1, regionMap);
    size_t mem2a = get_memory_usage_now();
    if(run_case == -1 || run_case == 1) {
      outputMatrix1 = rcp( new crs_matrix_type(regionMap, regionMap, 0) );
      outputMatrix1->doImport(*A, Importer, Tpetra::INSERT);
      outputMatrix1->fillComplete(rowMap1to1, rowMap1to1);
    }
    size_t mem2b = get_memory_usage_now();

    // 2) GDSW Proxy code
    if(run_case == -1 || run_case == 2) {
      tFunctions.importSquareMatrix(A, regionMap, outputMatrix2);
    }
    size_t mem3 = get_memory_usage_now();

    // 3) Import-and-FillComplete code
    // NOTE: This will *NOT* generate the same matrix because we're not specifying the column map
    // This is here only for memory use comparisons
    if(run_case == -1 || run_case == 3 || run_case == 5) 
      outputMatrix3 = Tpetra::importAndFillCompleteCrsMatrix<crs_matrix_type>(A,Importer,rowMap1to1,rowMap1to1);

    size_t mem4 = get_memory_usage_now();

    // 4) Import-based GDSW style
    if(run_case == -1 || run_case == 4)
      tFunctions.importSquareMatrixFromImporter(A, Teuchos::rcpFromRef(Importer), outputMatrix4);
    size_t mem5 = get_memory_usage_now();


    // 5) Locally filtered matrix
    if(run_case == -1 || run_case == 5)     
      outputMatrix5 = Filter(outputMatrix3,regionMap);   
    size_t mem6 = get_memory_usage_now();


    // 6) Import-based GDSW style, V2
    if(run_case == -1 || run_case == 6) 
      tFunctions.importSquareMatrixFromImporter2(A, Teuchos::rcpFromRef(Importer), outputMatrix6);
    size_t mem7 = get_memory_usage_now();

    // 7) Import-based GDSW style, V3
    if(run_case == -1 || run_case == 7) 
      tFunctions.importSquareMatrixFromImporter3(A, Teuchos::rcpFromRef(Importer), outputMatrix7);
    size_t mem8 = get_memory_usage_now();

    /***********************************************************************************/
    //std::cout<<"Breakdown Mem 0/1/2a/2b/3/4 = "<<mem0<<"/"<<mem1<<"/"<<mem2a<<"/"<<mem2b<<"/"<<mem3<<"/"<<mem4<<std::endl;

    std::cout<<"Orig matrix storage                           = "<< (mem1-mem0) <<std::endl;
    std::cout<<"Importer storage                              = "<< (mem2a-mem1) <<std::endl;
    if(run_case == -1 || run_case == 1) {
      std::cout<<"1) Tpetra Importer+Import+FC storage          = "<< (mem2b-mem2a) <<std::endl;
      names.push_back("Tpetra baseline");
      memory.push_back(mem2b-mem2a);
    }

    if(run_case == -1 || run_case == 2) {
      std::cout<<"2) GDSW-proxy storage                         = "<< (mem3-mem2b) <<std::endl;
      names.push_back("GDSW proxy");
      memory.push_back(mem3-mem2b);
    }

    if(run_case == -1 || run_case == 3 || run_case == 5) {
      std::cout<<"3) importAndFillComplete storage              = "<< (mem4-mem3) <<std::endl;
      names.push_back("importAndFillComplete");
      memory.push_back(mem4-mem3);
    }
    if(run_case == -1 || run_case == 4) {
      std::cout<<"4) Import-based GDSW-proxy                    = "<< (mem5-mem4) <<std::endl;
      names.push_back("Import-based GDSW-proxy");
      memory.push_back(mem5-mem4);
    }

    if(run_case == -1 || run_case == 5) {
      std::cout<<"5) Locally filtered IACF                      = "<< (mem6-mem5) <<std::endl;
      names.push_back("Locally filtered IACF");
      memory.push_back(mem6-mem5);
    }

    if(run_case == -1 || run_case == 6) {
      std::cout<<"6) V2 Import-based GDSW-proxy                 = "<< (mem7-mem6) <<std::endl;
      names.push_back("V2 Import-based GDSW");
      memory.push_back(mem2b-mem2a);
    }

    if(run_case == -1 || run_case == 7) {
      std::cout<<"7) V3 Import-based GDSW-proxy                 = "<< (mem8-mem7) <<std::endl;
      names.push_back("V3 Import-based GDSW");
      memory.push_back(mem2b-mem2a);
    }

    // Compare the output matrices
    bool result;
    if(run_case == -1) {
      result = compareCrsMatrix(*outputMatrix1,*outputMatrix2);
      if(!result) {
        out<<"*** Tpetra-based matrix ***"<<std::endl;
        outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
        out<<"*** GDSW-proxy matrix ***"<<std::endl;
        outputMatrix2->describe(out,Teuchos::VERB_EXTREME);
      }
      TEST_EQUALITY( result, true );
    }

    // Compare the output matrices
    if(run_case == -1) {
      result = compareCrsMatrix(*outputMatrix1,*outputMatrix2);
      if(!result) {
        out<<"*** Tpetra-based matrix ***"<<std::endl;
        outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
        out<<"*** GDSW-proxy matrix ***"<<std::endl;
        outputMatrix2->describe(out,Teuchos::VERB_EXTREME);
      }
      TEST_EQUALITY( result, true );
    }


    // Compare the output matrices
    if(run_case == -1) { 
      result = compareCrsMatrix(*outputMatrix1,*outputMatrix4);
      if(!result) {
        out<<"*** Tpetra-based matrix ***"<<std::endl;
        outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
        out<<"*** Import-based GDSW-proxy matrix ***"<<std::endl;
        outputMatrix4->describe(out,Teuchos::VERB_EXTREME);
      }
      TEST_EQUALITY( result, true );
    }      

    // Compare the output matrices
    if(run_case == -1) {
      result = compareCrsMatrix(*outputMatrix1,*outputMatrix5);
      if(!result) {
        out<<"*** Tpetra-based matrix ***"<<std::endl;
        outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
        out<<"*** Locally Filtered IAFC ***"<<std::endl;
        outputMatrix5->describe(out,Teuchos::VERB_EXTREME);
      }
      TEST_EQUALITY( result, true );
    }

    // Compare the output matrices
    if(run_case == -1) {
      result = compareCrsMatrix(*outputMatrix1,*outputMatrix6);
      if(!result) {
        out<<"*** Tpetra-based matrix ***"<<std::endl;
        outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
        out<<"*** V2 Import-based GDSW-proxy matrix ***"<<std::endl;
        outputMatrix6->describe(out,Teuchos::VERB_EXTREME);
      }
      TEST_EQUALITY( result, true );
    }

    // Compare the output matrices
    if(run_case == -1) {
      result = compareCrsMatrix(*outputMatrix1,*outputMatrix7);
      if(!result) {
        out<<"*** Tpetra-based matrix ***"<<std::endl;
        outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
        out<<"*** V3 Import-based GDSW-proxy matrix ***"<<std::endl;
        outputMatrix7->describe(out,Teuchos::VERB_EXTREME);
      }
      TEST_EQUALITY( result, true );
    }

    if(watchr_output) {
      out<<"Beginning Watchr output w/ "<<names.size()<<" records"<<std::endl;
      Tpetra::TestingXMLUtilities<double> XML;
      std::string xmlOut = XML.reportWatchrXML("memory","KB","GDSW Memory " + std::to_string(comm->getSize()) + " ranks",names,memory,comm);
      if(xmlOut.size() != 0)
      out<<"XML output in: "<<xmlOut<<std::endl;
      
    }
                          
  }

}// anonymous namespace

  //
  // INSTANTIATIONS
  //


int main(int narg, char *arg[]) 
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  using SC = Tpetra::Details::DefaultTypes::scalar_type; 
  using LO = Tpetra::Map<>::local_ordinal_type;
  using GO = Tpetra::Map<>::global_ordinal_type;
  using NT = Tpetra::Map<>::node_type;


  
  Teuchos::CommandLineProcessor cmdp(false,true);
  int run_case = -1; cmdp.setOption("case", &run_case, "Which case to run");
  bool watchr_output = false;  cmdp.setOption("watchr-output","no-watchr-output", &watchr_output, "Output memory data for watchr");

  if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }
  bool success = true;
  GDSWStyle_Test0<SC,LO,GO,NT>(run_case,watchr_output,success);


  if(!comm->getRank()) {
    if(success) std::cout<<"End Result: TEST PASSED"<<std::endl;
    else std::cout<<"End Result: TEST FAILED"<<std::endl;
  }
}

