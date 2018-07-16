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
#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include "TpetraExt_MatrixMatrix.hpp"
#include "Tpetra_MatrixIO.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrixMultiplyOp.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "MatrixMarket_Tpetra.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include <cmath>
#include "KokkosSparse_spgemm.hpp"
#include "Tpetra_RowMatrixTransposer.hpp"
#include "Tpetra_Import_Util2.hpp"

namespace {
  static const double defaultEpsilon = 1e-10;
  bool verbose = false;
  std::string matnamesFile;

  using Tpetra::MatrixMarket::Reader;
  using Tpetra::CrsMatrix;
  using Tpetra::CrsMatrixMultiplyOp;
  using Tpetra::global_size_t;
  using Tpetra::Map;
  using Tpetra::RowMatrixTransposer;
  using Tpetra::Vector;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::FancyOStream;
  using Teuchos::null;
  using Teuchos::outArg;
  using Teuchos::ParameterList;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::REDUCE_MIN;
  using Teuchos::reduceAll;
  using Teuchos::tuple;
  using std::endl;

TEUCHOS_STATIC_SETUP(){
  Teuchos::CommandLineProcessor &clp = Teuchos::UnitTestRepository::getCLP();
  clp.setOption("matnames-file", &matnamesFile,
    "A file containing a list of matricies we'll import", true);
  clp.setOption("v", "not-verbose", &verbose,
    "Whether or not to use verbose output");
}

template<class SC, class LO, class GO,class NT>
RCP<CrsMatrix<SC, LO, GO, NT> >
getIdentityMatrix (Teuchos::FancyOStream& out,
                   const Tpetra::global_size_t globalNumRows,
                   const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> Matrix_t;
  typedef Tpetra::Map<LO, GO, NT> Map_t;

  Teuchos::OSTab tab0 (out);
  out << "getIdentityMatrix" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Create row Map" << endl;
  RCP<const Map_t> identityRowMap =
    Tpetra::createUniformContigMapWithNode<LO, GO, NT> (globalNumRows, comm);

  out << "Create CrsMatrix" << endl;
  RCP<Matrix_t> identityMatrix =
    Tpetra::createCrsMatrix<SC, LO, GO, NT> (identityRowMap, 1);

  out << "Fill CrsMatrix" << endl;
  Teuchos::ArrayView<const GO> gblRows = identityRowMap->getNodeElementList ();
  for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
    Teuchos::Array<GO> col (1, *it);
    Teuchos::Array<SC> val (1, Teuchos::ScalarTraits<SC>::one ());
    identityMatrix->insertGlobalValues (*it, col (), val ());
  }

  out << "Call fillComplete" << endl;
  identityMatrix->fillComplete ();

  out << "Done!" << endl;
  return identityMatrix;
}

template<class SC, class LO, class GO,class NT>
RCP<CrsMatrix<SC, LO, GO, NT> >
getIdentityMatrixWithMap (Teuchos::FancyOStream& out,
			  Teuchos::RCP<const Tpetra::Map<LO,GO,NT> >& identityRowMap,
			  const Teuchos::RCP<const Teuchos::Comm<int> >& comm)
{
  using Teuchos::RCP;
  using std::endl;
  typedef Tpetra::CrsMatrix<SC, LO, GO, NT> Matrix_t;

  Teuchos::OSTab tab0 (out);
  out << "getIdentityMatrix" << endl;
  Teuchos::OSTab tab1 (out);

  out << "Create CrsMatrix" << endl;
  RCP<Matrix_t> identityMatrix =
    Tpetra::createCrsMatrix<SC, LO, GO, NT> (identityRowMap, 1);

  out << "Fill CrsMatrix" << endl;
  Teuchos::ArrayView<const GO> gblRows = identityRowMap->getNodeElementList ();
  for (auto it = gblRows.begin (); it != gblRows.end (); ++it) {
    Teuchos::Array<GO> col (1, *it);
    Teuchos::Array<SC> val (1, Teuchos::ScalarTraits<SC>::one ());
    identityMatrix->insertGlobalValues (*it, col (), val ());
  }

  out << "Call fillComplete" << endl;
  identityMatrix->fillComplete ();

  out << "Done!" << endl;
  return identityMatrix;
}



typedef struct add_test_results_struct{
  double correctNorm;
  double computedNorm;
  double epsilon;
} add_test_results;

typedef struct mult_test_results_struct{  
  mult_test_results_struct():epsilon(1e10),cNorm(-1.0),compNorm(-1.0){}

  double epsilon;
  double cNorm;
  double compNorm;
} mult_test_results;


template<class Matrix_t>
add_test_results regular_add_test(
    const std::string& name,
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    bool BT,
    RCP<Matrix_t > C,
    RCP<const Comm<int> > comm)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  add_test_results toReturn;
  toReturn.correctNorm = C->getFrobeniusNorm ();

  RCP<const Map_t > rowmap = AT ? A->getDomainMap() : A->getRowMap();
  RCP<Matrix_t> computedC = rcp( new Matrix_t(rowmap, 1));

  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::Add(*A, AT, one, *B, BT, one, computedC);
  computedC->fillComplete(A->getDomainMap(), A->getRangeMap());
  toReturn.computedNorm = computedC->getFrobeniusNorm ();
  toReturn.epsilon = fabs(toReturn.correctNorm - toReturn.computedNorm);

#if 0
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_calculated.mtx", computedC);
  Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile(
    name+"_real.mtx", C);
#endif


  return toReturn;
}

/// \brief Test the three-argument (A, B, C) version of CrsMatrix add,
///   where the output argument C is null on input.
///
/// \tparam Matrix_t A specialization of Tpetra::CrsMatrix.
template<class Matrix_t>
add_test_results
null_add_test (const Matrix_t& A,
               const Matrix_t& B,
               const bool AT,
               const bool BT,
               const Matrix_t& C,
               Teuchos::FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type scalar_type;
  typedef typename Matrix_t::local_ordinal_type local_ordinal_type;
  typedef typename Matrix_t::global_ordinal_type global_ordinal_type;
  typedef typename Matrix_t::node_type NT;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, NT> map_type;
  typedef Tpetra::Export<local_ordinal_type, global_ordinal_type, NT> export_type;
  const scalar_type one = STS::one ();

  out << "  Computing Frobenius norm of the expected result C" << endl;
  add_test_results toReturn;
  toReturn.correctNorm = C.getFrobeniusNorm ();

  out << "  Calling 3-arg add" << endl;
  RCP<const map_type> domainMap = BT ? B.getRangeMap () : B.getDomainMap ();
  RCP<const map_type> rangeMap = BT ? B.getDomainMap () : B.getRangeMap ();
  RCP<Matrix_t> C_computed =
    Tpetra::MatrixMatrix::add (one, AT, A, one, BT, B,
                               domainMap, rangeMap, Teuchos::null);
  TEUCHOS_TEST_FOR_EXCEPTION(
    C_computed.is_null (), std::logic_error, "3-arg add returned null.");

  RCP<Matrix_t> C_exported;
  if (! C_computed->getRowMap ()->isSameAs (* (C.getRowMap ()))) {
    // Export C_computed to C's row Map, so we can compare the two.
    export_type exp (C_computed->getRowMap (), C.getRowMap ());
    C_exported =
      Tpetra::exportAndFillCompleteCrsMatrix<Matrix_t> (C_computed, exp,
                                                             C.getDomainMap (),
                                                             C.getRangeMap ());
  } else {
    C_exported = C_computed;
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getRowMap ()->isSameAs (* (C.getRowMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different row Maps.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getDomainMap ()->isSameAs (* (C.getDomainMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different domain Maps.");
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! C_exported->getRangeMap ()->isSameAs (* (C.getRangeMap ())),
    std::logic_error,
    "Sorry, C_exported and C have different range Maps.");

  toReturn.computedNorm = C_exported->getFrobeniusNorm ();
  toReturn.epsilon = STS::magnitude (toReturn.correctNorm - toReturn.computedNorm);
  return toReturn;
}


template<class Matrix_t>
add_test_results add_into_test(
    RCP<Matrix_t > A,
    RCP<Matrix_t > B,
    bool AT,
    RCP<Matrix_t > C,
    RCP<const Comm<int> > comm)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;

  add_test_results toReturn;
  toReturn.correctNorm = C->getFrobeniusNorm ();

  RCP<const Map_t > rowmap =
    AT ? A->getDomainMap () : A->getRowMap ();
  RCP<Matrix_t> computedC = rcp (new Matrix_t (rowmap, 1));
  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::Add (*A, AT, one, *B, one);
  B->fillComplete ();
  toReturn.computedNorm = B->getFrobeniusNorm ();
  toReturn.epsilon = fabs (toReturn.correctNorm - toReturn.computedNorm);

  return toReturn;
}

template<class Matrix_t>
mult_test_results multiply_test_kernel(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  bool AT,
  bool BT,
  RCP<Matrix_t> C,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NO;
  typedef RowMatrixTransposer<SC,LO,GO,NO>  transposer_type;
  typedef Map<LO,GO,NO> Map_t;
  RCP<const Map_t> map = C->getRowMap();
  LO LO_INVALID = Teuchos::OrdinalTraits<LO>::invalid();

  // SC one = Teuchos::ScalarTraits<SC>::one();
  mult_test_results results;

  // Transposes
  RCP<Matrix_t> Aeff,  Beff;
  if(AT) {
    transposer_type transposer(A);
    Aeff = transposer.createTranspose();   
  }
  else
    Aeff=A;
  if(BT) {
    transposer_type transposer(B);
    Beff = transposer.createTranspose();
  }
  else
    Beff=B;


  // Now let's handle incompatibilities between the cols of Aeff and the rows of Beff, by copying and rearranging Beff if needed
  if(!Aeff->getGraph()->getColMap()->isSameAs(*Beff->getGraph()->getColMap())) {
    Teuchos::ArrayRCP<const size_t> Be1_rowptr;
    Teuchos::ArrayRCP<const LO> Be1_colind;
    Teuchos::ArrayRCP<const SC> Be1_vals;
    Beff->getAllValues(Be1_rowptr,Be1_colind,Be1_vals);

    RCP<const Map_t> Be1_rowmap = Beff->getGraph()->getRowMap();
    RCP<const Map_t> Be2_rowmap = Aeff->getGraph()->getColMap();

    // Do rowptr
    Teuchos::ArrayRCP<size_t> Be2_rowptr(Be2_rowmap->getNodeNumElements()+1);
    Be2_rowptr[0]=0;
    for(size_t i=0; i<Be2_rowmap->getNodeNumElements(); i++) {
      LO lrid_1 = Be1_rowmap->getLocalElement(Be2_rowmap->getGlobalElement(i));
      if(lrid_1 == LO_INVALID)
	Be2_rowptr[i+1] = Be2_rowptr[i];
      else {
	Be2_rowptr[i+1] = Be2_rowptr[i] + Be1_rowptr[lrid_1+1] - Be1_rowptr[lrid_1];
      }
    }
    int nnz_2 =  Be2_rowptr[Be2_rowmap->getNodeNumElements()];

    Teuchos::ArrayRCP<LO> Be2_colind(nnz_2);
    Teuchos::ArrayRCP<SC> Be2_vals(nnz_2);

    // Copy colind/vals
    for(size_t i=0; i<Be2_rowmap->getNodeNumElements(); i++) {
      LO lrid_1 = Be1_rowmap->getLocalElement(Be2_rowmap->getGlobalElement(i));      
      // NOTE: Invalid rows will be length zero, so this will always work
      for(size_t j=Be2_rowptr[i]; j<Be2_rowptr[i+1]; j++) {
	size_t j1 = j-Be2_rowptr[i] + Be1_rowptr[lrid_1];
	Be2_colind[j] = Be1_colind[j1];
	Be2_vals[j]   = Be1_vals[j1];
      }      
    }

    // Make new Beff
    RCP<Matrix_t> Beff2 = rcp(new Matrix_t(Be2_rowmap,Beff->getGraph()->getColMap(),Be2_rowptr,Be2_colind,Be2_vals));
			      //    Beff2->setAllValues(
    Beff2->expertStaticFillComplete(Beff->getGraph()->getDomainMap(),Beff->getGraph()->getRangeMap());//Not efficient w.r.t. import, but that's OK
    Beff=Beff2;
  }

  // Extract Kokkos CrsMatrices
  typedef typename Tpetra::CrsMatrix<SC,LO,GO,NO>::local_matrix_type KCRS;
  typedef typename KCRS::device_type device_t;
  typedef typename KCRS::StaticCrsGraphType graph_t;
  typedef typename graph_t::row_map_type::non_const_type lno_view_t;
  typedef typename graph_t::entries_type::non_const_type lno_nnz_view_t;
  typedef typename KCRS::values_type::non_const_type scalar_view_t;
  typedef KokkosKernels::Experimental::KokkosKernelsHandle<typename lno_view_t::const_value_type,typename lno_nnz_view_t::const_value_type, typename scalar_view_t::const_value_type, 
          typename device_t::execution_space, typename device_t::memory_space,typename device_t::memory_space > KernelHandle;

  // Grab the  Kokkos::SparseCrsMatrix-es
  const KCRS & Ak = Aeff->getLocalMatrix();
  const KCRS & Bk = Beff->getLocalMatrix();

  // Setup
  // As a note "SPGEMM_MKL" will *NOT* pass all of these tests
  //  std::vector<std::string> ALGORITHMS={"SPGEMM_MKL","SPGEMM_KK_MEMSPEED","SPGEMM_KK_SPEED","SPGEMM_KK_MEMORY"};
  std::vector<std::string> ALGORITHMS={"SPGEMM_KK_MEMORY"};

  for(int alg = 0; alg < (int)ALGORITHMS.size(); alg++) {
    std::string myalg = ALGORITHMS[alg];
    printf("Testing algorithm %s on %s\n",myalg.c_str(),name.c_str());
    
    KokkosSparse::SPGEMMAlgorithm alg_enum = KokkosSparse::StringToSPGEMMAlgorithm(myalg);
    typename KernelHandle::nnz_lno_t AnumRows = Ak.numRows();
    typename KernelHandle::nnz_lno_t BnumRows = Bk.numRows();
    typename KernelHandle::nnz_lno_t BnumCols = Bk.numCols();  
    lno_view_t      row_mapC ("non_const_lnow_row", AnumRows + 1);
    lno_nnz_view_t  entriesC;
    scalar_view_t   valuesC;
    KernelHandle kh;
    kh.create_spgemm_handle(alg_enum);
    
    // Symbolic
    KokkosSparse::Experimental::spgemm_symbolic(&kh,AnumRows,BnumRows,BnumCols,Ak.graph.row_map,Ak.graph.entries,false,Bk.graph.row_map,Bk.graph.entries,false,row_mapC);
    size_t c_nnz_size = kh.get_spgemm_handle()->get_c_nnz();
    //    printf("DEBUG: c_nnz_size = %d\n",c_nnz_size); // This isn't relevant for MKL
    if (c_nnz_size){
      entriesC = lno_nnz_view_t (Kokkos::ViewAllocateWithoutInitializing("entriesC"), c_nnz_size);
      valuesC = scalar_view_t (Kokkos::ViewAllocateWithoutInitializing("valuesC"), c_nnz_size);
    }
    
    // Numeric
    KokkosSparse::Experimental::spgemm_numeric(&kh,AnumRows,BnumRows,BnumCols,Ak.graph.row_map,Ak.graph.entries,Ak.values,false,Bk.graph.row_map,Bk.graph.entries,Bk.values,false,row_mapC,entriesC,valuesC);
    kh.destroy_spgemm_handle();
    
    // Sort
    Tpetra::Import_Util::sortCrsEntries(row_mapC, entriesC, valuesC);

    // Compare the returned arrays with that of actual C
    Teuchos::ArrayRCP<const size_t> Real_rowptr;
    Teuchos::ArrayRCP<const LO> Real_colind;
    Teuchos::ArrayRCP<const SC> Real_vals;
    C->getAllValues(Real_rowptr,Real_colind,Real_vals);

    // Check number of rows
    if((size_t)Real_rowptr.size() != (size_t)row_mapC.size()) throw std::runtime_error("mult_test_results multiply_test_kernel: rowmap size mismatch");
  
    // Check row sizes
    bool has_mismatch=false;
    for(size_t i=0; i<(size_t) Real_rowptr.size(); i++) {
      if(Real_rowptr()[i] !=row_mapC[i]) {has_mismatch=true;break;}
    }
    if(has_mismatch) {
#if 0     
      if(AT) {
	std::string fname(name + "_AT.out");
	Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile (fname,Aeff,"AT","AT");
      }
      if(BT) {
	std::string fname(name + "_BT.out");
	Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile (fname,Beff,"BT","BT");
      }

      printf("Real rowgids = ");
      for(size_t i=0; i<(size_t) C->getGraph()->getRowMap()->getNodeNumElements(); i++) 
	printf("%d ",(int)C->getGraph()->getRowMap()->getGlobalElement(i));
      printf("\n");
      printf("Aeff rowgids = ");
      for(size_t i=0; i<(size_t) Aeff->getGraph()->getRowMap()->getNodeNumElements(); i++) 
	printf("%d ",(int)Aeff->getGraph()->getRowMap()->getGlobalElement(i));
      printf("\n");

      printf("Real rowptr = ");
      for(size_t i=0; i<(size_t) Real_rowptr.size(); i++) 
	printf("%d ",(int)Real_rowptr()[i]);
      printf("\nCalc rowptr = ");
      for(size_t i=0; i<(size_t) row_mapC.size(); i++) 
	printf("%d ",(int)row_mapC[i]);
      printf("\n");
      
      printf("Real colind = ");
      for(size_t i=0; i<(size_t) Real_colind.size(); i++) 
	printf("%d(%6.1f) ",(int)Real_colind()[i],Real_vals()[i]);
      printf("\nCalc colind = ");
      for(size_t i=0; i<(size_t) entriesC.size(); i++) 
	printf("%d(%6.1f) ",(int)entriesC[i],valuesC[i]);
      printf("\n");
#endif
      throw std::runtime_error("mult_test_results multiply_test_kernel: rowmap entries mismatch");
    } 

    // Check graphs & entries (sorted by GID)
    results.cNorm=0.0;
    results.compNorm=0.0;
    results.epsilon=0.0;
    for(size_t i=0; i<(size_t) Real_rowptr.size()-1; i++) {      
      size_t nnz = Real_rowptr()[i+1] - Real_rowptr()[i];
      if(nnz==0) continue;
      std::vector<std::pair<LO,SC> > R_sorted(nnz), C_sorted(nnz);
      for(size_t j=0; j<nnz; j++)  {
	R_sorted[j].first  = C->getColMap()->getGlobalElement(Real_colind()[Real_rowptr()[i]+j]);
	R_sorted[j].second = Real_vals()[Real_rowptr()[i]+j];
	C_sorted[j].first  = Beff->getColMap()->getGlobalElement(entriesC[Real_rowptr()[i]+j]);
	C_sorted[j].second = valuesC[Real_rowptr()[i]+j];
      }
      std::sort(R_sorted.begin(),R_sorted.end());
      std::sort(C_sorted.begin(),C_sorted.end());
      for(size_t j=0; j<nnz; j++) {
	if(R_sorted[j].first !=C_sorted[j].first) {has_mismatch=true;break;}
	results.cNorm+=std::abs(R_sorted[j].second);
	results.compNorm+=std::abs(C_sorted[j].second);
	results.epsilon+=std::abs(C_sorted[j].second-R_sorted[j].second);
      }
    }
    if(has_mismatch) {
#if 0     
      printf("Real colmap = ");
      for(size_t i=0; i<(size_t) C->getGraph()->getColMap()->getNodeNumElements(); i++) 
	printf("%d ",(int)C->getGraph()->getColMap()->getGlobalElement(i));
      printf("\n");
      printf("   B colmap = ");
      for(size_t i=0; i<(size_t) Beff->getGraph()->getColMap()->getNodeNumElements(); i++) 
	printf("%d ",(int)Beff->getGraph()->getColMap()->getGlobalElement(i));
      printf("\n");
 
      if(AT) {
	std::string fname(name + "_AT.out");
	Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile (fname,Aeff,"AT","AT");
      }
      if(BT) {
	std::string fname(name + "_BT.out");
	Tpetra::MatrixMarket::Writer<Matrix_t>::writeSparseFile (fname,Beff,"BT","BT");
      }
      printf("Real rowptr = ");
      for(size_t i=0; i<(size_t) Real_rowptr.size(); i++) 
	printf("%d ",(int)Real_rowptr()[i]);
      printf("\nCalc rowptr = ");
      for(size_t i=0; i<(size_t) row_mapC.size(); i++) 
	printf("%d ",(int)row_mapC[i]);
      printf("\n");
      
      printf("Real colind = ");
      for(size_t i=0; i<(size_t) Real_colind.size(); i++) 
	printf("%d(%6.1f) ",(int)Real_colind()[i],Real_vals()[i]);
      printf("\nCalc colind = ");
      for(size_t i=0; i<(size_t) entriesC.size(); i++) 
	printf("%d(%6.1f) ",(int)entriesC[i],valuesC[i]);
      printf("\n");
#endif
      throw std::runtime_error("mult_test_results multiply_test_kernel: colmap entries mismatch");
    } 
    
    if(results.cNorm>1e-10) results.epsilon = results.epsilon / results.cNorm;

    if(results.epsilon>1e-10)
      throw std::runtime_error("mult_test_results multiply_test_kernel: values mismatch");

  }


  return results;
}


template<class Matrix_t>
mult_test_results multiply_reuse_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  bool AT,
  bool BT,
  RCP<Matrix_t> C,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Map<LO,GO,NT> Map_t;
  typedef Vector<SC,LO,GO,NT> Vector_t;

  RCP<const Map_t> map = C->getRowMap();

  // Scaling vectors
  Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);

  RCP<Vector_t> leftScaling  = rcp( new Vector_t(C->getRangeMap()) );
  leftScaling->randomize();
  leftScaling->norm2(norms);
  leftScaling->scale(1.0/norms[0]);

  RCP<Vector_t> rightScaling = rcp( new Vector_t(C->getDomainMap()) );
  rightScaling->randomize();
  rightScaling->norm2(norms);
  rightScaling->scale(1.0/norms[0]);

  // computedC1 = leftScaling * (op(A)*op(B)) * rightScaling
  RCP<Matrix_t> computedC1 = rcp( new Matrix_t(map, 1));
  Tpetra::MatrixMatrix::Multiply(*A, AT, *B, BT, *computedC1, false/*call_FillCompleteOnResult*/);
  computedC1->fillComplete(C->getDomainMap(), C->getRangeMap());
  computedC1->leftScale (*leftScaling);
  computedC1->rightScale(*rightScaling);

  // NOTE (mfh 05 Jun 2016) This may even be null.  It exists at this
  // point only for the syntax.
  RCP<NT> node = map->getNode ();

  // As = leftScaling * op(A) =
  //   leftScaling * A, if AT=false
  //   A*leftScaling,   if AT=true
  RCP<Matrix_t> As = A->clone(node);
  if (AT == false) As->leftScale (*leftScaling);
  else             As->rightScale(*leftScaling);

  // Bs = op(B) * rightScaling =
  //   B * rightScaling, if BT=false
  //   rightScaling*B,   if BT=true
  RCP<Matrix_t> Bs = B->clone(node);
  if (BT == false) Bs->rightScale(*rightScaling);
  else             Bs->leftScale (*rightScaling);

  // computedC2 = op(As) * op(Bs)
  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(C->getCrsGraph()) );
  computedC2->fillComplete(C->getDomainMap(), C->getRangeMap());

  computedC2->resumeFill();
  Tpetra::MatrixMatrix::Multiply(*As, AT, *Bs, BT, *computedC2, true/*call_FillCompleteOnResult*/);

  // diffMatrix = computedC2 - computedC1
  SC one = Teuchos::ScalarTraits<SC>::one();
  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(C->getRowMap());
  Tpetra::MatrixMatrix::Add(*computedC1, false, -one, *computedC2, false, one, diffMatrix);
  diffMatrix->fillComplete(C->getDomainMap(), C->getRangeMap());

  mult_test_results results;
  results.cNorm    = C->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  return results;
}


template<class Matrix_t>
mult_test_results jacobi_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Vector<SC,LO,GO,NT> Vector_t;
  typedef Map<LO,GO,NT> Map_t;
  RCP<const Map_t> map = A->getRowMap();

  SC omega=Teuchos::ScalarTraits<SC>::one();
  Vector_t Dinv(B->getRowMap());
  Dinv.putScalar(1.0);

  // Jacobi version
  RCP<Matrix_t> C = rcp(new Matrix_t(B->getRowMap(),0));
  Tpetra::MatrixMatrix::Jacobi(omega,Dinv,*A,*B,*C);

  // Multiply + Add version
  Dinv.putScalar(omega);
  RCP<Matrix_t> AB = rcp(new Matrix_t(B->getRowMap(),0));
  RCP<Matrix_t> C_check = rcp(new Matrix_t(B->getRowMap(),0));
  Tpetra::MatrixMatrix::Multiply(*A,false,*B,false,*AB);
  AB->leftScale(Dinv);
  SC one = Teuchos::ScalarTraits<SC>::one();
  Tpetra::MatrixMatrix::Add(*AB,false,-one,*B,false,one,C_check);

  // Check the difference
  Tpetra::MatrixMatrix::Add(*C, false, -one, *C_check, one);
  C_check->fillComplete(B->getDomainMap(),B->getRangeMap());

  // Error Check
  double compNorm = C_check->getFrobeniusNorm();
  double cNorm = C->getFrobeniusNorm();
  mult_test_results results;
  results.epsilon  = compNorm/cNorm;
  results.cNorm    = cNorm;
  results.compNorm = compNorm;
  return results;
}


template<class Matrix_t>
mult_test_results jacobi_reuse_test(
  const std::string& name,
  RCP<Matrix_t> A,
  RCP<Matrix_t> B,
  RCP<const Comm<int> > comm,
  FancyOStream& out)
{
  typedef typename Matrix_t::scalar_type SC;
  typedef typename Matrix_t::local_ordinal_type LO;
  typedef typename Matrix_t::global_ordinal_type GO;
  typedef typename Matrix_t::node_type NT;
  typedef Vector<SC,LO,GO,NT> Vector_t;
  typedef Map<LO,GO,NT> Map_t;

  RCP<const Map_t> map = A->getRowMap();

  // Scaling vectors
  Teuchos::Array<typename Teuchos::ScalarTraits<SC>::magnitudeType> norms(1);
  RCP<Vector_t> rightScaling = rcp( new Vector_t(B->getDomainMap()) );
  rightScaling->randomize();
  rightScaling->norm2(norms);
  rightScaling->scale(1.0/norms[0]);
  SC one = Teuchos::ScalarTraits<SC>::one();

  SC omega = one;
  Vector_t Dinv(B->getRowMap());
  Dinv.putScalar(1.0);

  // Jacobi version
  RCP<Matrix_t> computedC1 = rcp(new Matrix_t(B->getRowMap(), 0));
  Tpetra::MatrixMatrix::Jacobi(omega, Dinv, *A, *B, *computedC1);
  computedC1->rightScale(*rightScaling);

  // Bs = B * rightScaling

  // NOTE (mfh 05 Jun 2016) This may even be null.  It exists at this
  // point only for the syntax.
  RCP<NT> node = map->getNode ();
  RCP<Matrix_t> Bs = B->clone(node);
  Bs->rightScale(*rightScaling);

  // computedC2 = (I - Dinv*A)*Bs
  RCP<Matrix_t> computedC2 = rcp( new Matrix_t(computedC1->getCrsGraph()) );
  Tpetra::MatrixMatrix::Jacobi(omega, Dinv, *A, *Bs, *computedC2);

  // diffMatrix = computedC2 - computedC1

  RCP<Matrix_t> diffMatrix = Tpetra::createCrsMatrix<SC,LO,GO,NT>(computedC1->getRowMap());
  Tpetra::MatrixMatrix::Add(*computedC1, false, -one, *computedC2, false, one, diffMatrix);
  diffMatrix->fillComplete(computedC1->getDomainMap(), computedC1->getRangeMap());

  mult_test_results results;
  results.cNorm    = computedC1->getFrobeniusNorm ();
  results.compNorm = diffMatrix->getFrobeniusNorm ();
  results.epsilon  = results.compNorm/results.cNorm;

  return results;
}


TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Tpetra_MatMatKernels, operations_test,SC,LO, GO, NT)  {
  RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
  
  // NOTE: The matrix reader doesn't read real matrices into a complex data type, so we just swap down to MT here
  typedef typename Teuchos::ScalarTraits<SC>::magnitudeType MT;
  typedef CrsMatrix<MT,LO,GO,NT> Matrix_t;
  const int myRank = comm->getRank ();
  //const int numProcs = comm->getSize();

  // Create an output stream that prints immediately, rather than
  // waiting until the test finishes.  This lets us get debug output
  // even if an unexpected exception or crash happens.  Unfortunately,
  // Teuchos::FancyOStream does not implement operator=, so we can't
  // replace 'out' temporarily.
  RCP<FancyOStream> newOutP = (myRank == 0) ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::getFancyOStream (rcp (new Teuchos::oblackholestream ()));
  FancyOStream& newOut = *newOutP;

  newOut << "Tpetra sparse matrix-matrix multiply: operations_test" << endl;
  Teuchos::OSTab tab1 (newOut);

  newOut << "Get parameters from XML file" << endl;
  Teuchos::RCP<Teuchos::ParameterList> matrixSystems =
    Teuchos::getParametersFromXmlFile(matnamesFile);

  for (Teuchos::ParameterList::ConstIterator it = matrixSystems->begin();
       it != matrixSystems->end();
       ++it) {
    TEUCHOS_TEST_FOR_EXCEPTION(!it->second.isList(), std::runtime_error,
      "All top-level items in the list of matrix file names must be "
      "ParameterList instances.  " << endl <<
      "Bad tag's name: " << it->first <<
      "Type name: " << it->second.getAny().typeName() <<
      endl << endl);
 
    ParameterList currentSystem = matrixSystems->sublist (it->first);
    std::string name = currentSystem.name();
    std::string A_file = currentSystem.get<std::string> ("A");
    std::string B_file = currentSystem.get<std::string> ("B");
    std::string C_file = currentSystem.get<std::string> ("C");
    const bool AT = currentSystem.get<bool> ("TransA");
    const bool BT = currentSystem.get<bool> ("TransB");
    double epsilon = currentSystem.get<double> ("epsilon", defaultEpsilon);
    std::string op = currentSystem.get<std::string> ("op");

    RCP<Matrix_t> A, B, C;
    if (A_file != "")
      A = Reader<Matrix_t>::readSparseFile (A_file, comm);
    if (B_file != "")
      B = Reader<Matrix_t>::readSparseFile (B_file, comm);
    if (C_file != "")
      C = Reader<Matrix_t>::readSparseFile (C_file, comm);

    TEUCHOS_TEST_FOR_EXCEPTION(op != "multiply" && op != "add" && op != "RAP", std::runtime_error,
                               "Unrecognized Matrix Operation: " << op);

    if (op == "multiply") {
      if (verbose)
        newOut << "Running multiply test (kernel) for " << currentSystem.name() << endl;
      mult_test_results results;	
      try {
	results = multiply_test_kernel(name, A, B, AT, BT, C, comm, newOut);
      }
      catch(const std::exception & err) {
	printf("FAILED: %s\n",err.what());
      }
      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);
    }    
      // NOTE: Re-enable these tests when we get kernels for these things
      /*
      if (verbose)
        newOut << "Running multiply reuse test for " << currentSystem.name() << endl;

      results = multiply_reuse_test(name, A, B, AT, BT, C, comm, newOut);

      if (verbose) {
        newOut << "Results:"     << endl;
        newOut << "\tEpsilon: "  << results.epsilon  << endl;
        newOut << "\tcNorm: "    << results.cNorm    << endl;
        newOut << "\tcompNorm: " << results.compNorm << endl;
      }
      TEST_COMPARE(results.epsilon, <, epsilon);


      // Do we try Jacobi?
      if (AT == false && BT == false && A->getDomainMap()->isSameAs(*A->getRangeMap())) {
        if (verbose)
          newOut << "Running jacobi test for " << currentSystem.name() << endl;

        results = jacobi_test(name, A, B, comm, newOut);
        if (verbose) {
          newOut << "Results:"     << endl;
          newOut << "\tEpsilon: "  << results.epsilon  << endl;
          newOut << "\tcNorm: "    << results.cNorm    << endl;
          newOut << "\tcompNorm: " << results.compNorm << endl;
        }
        TEST_COMPARE(results.epsilon, <, epsilon)

        if (verbose)
          newOut << "Running jacobi reuse test for " << currentSystem.name() << endl;

        results = jacobi_reuse_test(name, A, B, comm, newOut);

        if (verbose) {
          newOut << "Results:"     << endl;
          newOut << "\tEpsilon: "  << results.epsilon  << endl;
          newOut << "\tcNorm: "    << results.cNorm    << endl;
          newOut << "\tcompNorm: " << results.compNorm << endl;
        }
        TEST_COMPARE(results.epsilon, <, epsilon)
      }
    }
    else if (op == "add") {
      if (verbose)
        newOut << "Running 3-argument add test (nonnull C on input) for "
               << currentSystem.name() << endl;

      add_test_results results = regular_add_test(name+"_add",A, B, AT, BT, C, comm);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Regular Add Test Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      // FIXME (mfh 09 May 2013) This test doesn't currently pass.  I
      // don't think anyone ever exercised the case where C is null on
      // input before.  I'm disabling this test for now until I have a
      // chance to fix that case.
      if (verbose)
        newOut << "Running 3-argument add test (null C on input) for "
               << currentSystem.name() << endl;

      TEUCHOS_TEST_FOR_EXCEPTION(A.is_null (), std::logic_error,
                                 "Before null_add_test: A is null");
      TEUCHOS_TEST_FOR_EXCEPTION(B.is_null (), std::logic_error,
                                 "Before null_add_test: B is null");
      TEUCHOS_TEST_FOR_EXCEPTION(C.is_null (), std::logic_error,
                                 "Before null_add_test: C is null");

      results = null_add_test<Matrix_t> (*A, *B, AT, BT, *C, newOut);

      TEST_COMPARE(results.epsilon, <, epsilon);
      newOut << "Null Add Test Results: " << endl;
      newOut << "\tCorrect Norm: " << results.correctNorm << endl;
      newOut << "\tComputed norm: " << results.computedNorm << endl;
      newOut << "\tEpsilon: " << results.epsilon << endl;

      B = Reader<Matrix_t >::readSparseFile(B_file, comm, false);

      if (! BT) {
        if (verbose)
          newOut << "Running 2-argument add test for "
                 << currentSystem.name() << endl;

        results = add_into_test(A,B,AT,C,comm);
        TEST_COMPARE(results.epsilon, <, epsilon)
        newOut << "Add Into Test Results: " << endl;
        newOut << "\tCorrect Norm: " << results.correctNorm << endl;
        newOut << "\tComputed norm: " << results.computedNorm << endl;
        newOut << "\tEpsilon: " << results.epsilon << endl;
      }
    }   
      */
  } 

  const int lclSuccess = success ? 1 : 0;
  int gblSuccess = 0;
  reduceAll<int, int> (*comm, REDUCE_MIN, lclSuccess, outArg (gblSuccess));
  newOut << "We made it through operations_test on all processes!" << endl;
  if (gblSuccess != 1) {
    newOut << "FAILED on at least one process!" << endl;
  }
  TEST_EQUALITY_CONST( gblSuccess, 1 );
}



#define UNIT_TEST_GROUP_SC_LO_GO_NO( SC, LO, GO, NT )			\
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Tpetra_MatMatKernels, operations_test,SC, LO, GO, NT)


  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_SLGN_NO_ORDINAL_SCALAR( UNIT_TEST_GROUP_SC_LO_GO_NO )


  } //namespace Tpetra

