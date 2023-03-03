/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/

#include <Tpetra_TestingUtilities.hpp>
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


// Yes, I'm including the CPP file on purpose.  Don't hate.
#include "GDSW_Proxy.hpp"
#include "GDSW_Proxy.cpp"

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
  typedef Tpetra::global_size_t GST;
  typedef Tpetra::Map<LO, GO, NT> map_type;
  typedef Tpetra::Export<LO, GO, NT> export_type;
  typedef Tpetra::CrsMatrix<ST, LO, GO, NT> matrix_type;

  RCP<const Teuchos::Comm<int> > comm = rowMap->getComm ();

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



// return current memory usage in bytes
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
  memory = (unsigned long)sys_resources.ru_maxrss * RU_MAXRSS_UNITS;
#endif  


  /* Success */
  return memory;
}


}//anonymous namespace


  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( GDSWStyle, Test0, SC, LO, GO, NT ) {
    using map_type = Map<LO, GO, NT>;
    using import_type = Tpetra::Import<LO, GO, NT>;
    using crs_matrix_type = Tpetra::CrsMatrix<SC,LO,GO,NT>;
    bool debug=true;

    const Tpetra::global_size_t INVALID =  Teuchos::OrdinalTraits<GO>::invalid ();
    RCP<const Comm<int> > comm = Tpetra::getDefaultComm();

    // create map
    RCP<const map_type > rowMap1to1 =  createContigMapWithNode<LO, GO, NT> (INVALID, 10000, comm);
    size_t mem0 = get_memory_usage_now();


    // Create matrix
    RCP<crs_matrix_type> A = createLaplace1D<SC,LO,GO,NT>(rowMap1to1,out,debug);

    // Grab the column map to create a proxy for the input region
    RCP<const map_type> regionMap = A->getColMap();

    size_t mem1 = get_memory_usage_now();


    // 1) Tpetra-based code duplication what is in GDSW
    import_type Importer(rowMap1to1, regionMap);
    size_t mem2a = get_memory_usage_now();
    RCP<crs_matrix_type> outputMatrix1 = rcp( new crs_matrix_type(regionMap, regionMap, 0) );
    outputMatrix1->doImport(*A, Importer, Tpetra::INSERT);
    outputMatrix1->fillComplete(rowMap1to1, rowMap1to1);
    size_t mem2b = get_memory_usage_now();


    // 2) GDSW Proxy code
    RCP<crs_matrix_type> outputMatrix2;
    TpetraFunctions<SC,LO,GO,NT> tFunctions;
    tFunctions.importSquareMatrix(A, regionMap, outputMatrix2);
    size_t mem3 = get_memory_usage_now();

    // 3) Import-and-FillComplete code
    // NOTE: This will *NOT* generate the same matrix because we're not specifying the column map
    // This is here only for memory use comparisons
    RCP<crs_matrix_type> outputMatrix3 = Tpetra::importAndFillCompleteCrsMatrix<crs_matrix_type>(A,Importer,rowMap1to1,rowMap1to1);

    size_t mem4 = get_memory_usage_now();

    // 4) Import-based GDSW style
    RCP<crs_matrix_type> outputMatrix4;
    tFunctions.importSquareMatrixFromImporter(A, Teuchos::rcpFromRef(Importer), outputMatrix4);
    size_t mem5 = get_memory_usage_now();


    /***********************************************************************************/
    //std::cout<<"Breakdown Mem 0/1/2a/2b/3/4 = "<<mem0<<"/"<<mem1<<"/"<<mem2a<<"/"<<mem2b<<"/"<<mem3<<"/"<<mem4<<std::endl;

    std::cout<<"Orig matrix storage                           = "<< (mem1-mem0) <<std::endl;
    std::cout<<"Importer storage                              = "<< (mem2a-mem1) <<std::endl;
    std::cout<<"1) Tpetra Importer+Import+FC storage          = "<< (mem2b-mem1) <<std::endl;
    std::cout<<"2) GDSW-proxy storage                         = "<< (mem3-mem2b) <<std::endl;
    std::cout<<"3) importAndFillComplete storage              = "<< (mem4-mem3) <<std::endl;
    std::cout<<"4) Import-based GDSW-proxy                    = "<< (mem5-mem4) <<std::endl;

    // Compare the output matrices
    bool result = compareCrsMatrix(*outputMatrix1,*outputMatrix2);
    if(!result) {
      out<<"*** Tpetra-based matrix ***"<<std::endl;
      outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
      out<<"*** GDSW-proxy matrix ***"<<std::endl;
      outputMatrix2->describe(out,Teuchos::VERB_EXTREME);
    }
    TEST_EQUALITY( result, true );

    
    // Compare the output matrices
    result = compareCrsMatrix(*outputMatrix1,*outputMatrix4);
    if(!result) {
      out<<"*** Tpetra-based matrix ***"<<std::endl;
      outputMatrix1->describe(out,Teuchos::VERB_EXTREME);
      out<<"*** Import-based GDSW-proxy matrix ***"<<std::endl;
      outputMatrix4->describe(out,Teuchos::VERB_EXTREME);
    }
    TEST_EQUALITY( result, true );

  }


  //
  // INSTANTIATIONS
  //


#define UNIT_TEST_4( SCALAR, LO, GO, NT )                                \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( GDSWStyle, Test0, SCALAR, LO, GO, NT ) 


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN( UNIT_TEST_4 )

} // namespace (anonymous)


