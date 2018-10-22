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

#include <Tpetra_ConfigDefs.hpp>
#include <Tpetra_TestingUtilities.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <map>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_as.hpp>
#include <Teuchos_Tuple.hpp>
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Distributor.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Util.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_Details_gathervPrint.hpp"



// Macro that marks a function as "possibly unused," in order to
// suppress build warnings.
#if ! defined(TRILINOS_UNUSED_FUNCTION)
#  if defined(__GNUC__) || (defined(__INTEL_COMPILER) && !defined(_MSC_VER))
#    define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#  elif defined(__clang__)
#    if __has_attribute(unused)
#      define TRILINOS_UNUSED_FUNCTION __attribute__((__unused__))
#    else
#      define TRILINOS_UNUSED_FUNCTION
#    endif // Clang has 'unused' attribute
#  elif defined(__IBMCPP__)
// IBM's C++ compiler for Blue Gene/Q (V12.1) implements 'used' but not 'unused'.
//
// http://pic.dhe.ibm.com/infocenter/compbg/v121v141/index.jsp
#    define TRILINOS_UNUSED_FUNCTION
#  else // some other compiler
#    define TRILINOS_UNUSED_FUNCTION
#  endif
#endif // ! defined(TRILINOS_UNUSED_FUNCTION)


namespace {

  template<class T1, class T2>
  void vector_check(size_t N, T1 & v1, T2 & v2) {
    int myRank = v1.getMap()->getComm()->getRank();
    for(size_t i=0; i<N; i++)
      if(v1.getDataNonConst(0)[i] != v2.getDataNonConst(0)[i]) {
        std::stringstream oss;
        oss<<"["<<myRank<<"] Error: Mismatch on unknown "<<i<<" "<<v1.getDataNonConst(0)[i]<<" != "<<v2.getDataNonConst(0)[i]<<std::endl;
        throw std::runtime_error(oss.str());
      }
  }


  // Teuchos using statements
   using Tpetra::TestingUtilities::getDefaultComm;

  using std::endl;
  using std::copy;
  using std::ostream_iterator;
  using std::string;

  using Teuchos::RCP;
  using Teuchos::ArrayRCP;
  using Teuchos::rcp;
  using Teuchos::null;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::Comm;
  using Teuchos::Range1D;
  using Teuchos::Tuple;
  using Teuchos::as;
  using Teuchos::OrdinalTraits;
  using Teuchos::ScalarTraits;
  using Teuchos::arrayView;
  using Teuchos::tuple;
  using Teuchos::NO_TRANS;
  using Teuchos::TRANS;
  using Teuchos::CONJ_TRANS;
  using Teuchos::VERB_DEFAULT;
  using Teuchos::VERB_NONE;
  using Teuchos::VERB_LOW;
  using Teuchos::VERB_MEDIUM;
  using Teuchos::VERB_HIGH;
  using Teuchos::VERB_EXTREME;

  using Tpetra::Map;
  using Tpetra::MultiVector;
  using Tpetra::global_size_t;
  using Tpetra::GloballyDistributed;
  typedef Tpetra::global_size_t GST;

  
  using Tpetra::createContigMapWithNode;
  using Tpetra::createLocalMapWithNode;

  double errorTolSlack = 1.0e+2;



   template <class Scalar>
  typename Teuchos::ScalarTraits<Scalar>::magnitudeType testingTol() { return Teuchos::ScalarTraits<Scalar>::eps(); }
  template <>
  TRILINOS_UNUSED_FUNCTION int testingTol<int>() { return 0; }
  template <>
  TRILINOS_UNUSED_FUNCTION long testingTol<long>() { return 0; }
  template <>
  TRILINOS_UNUSED_FUNCTION long long testingTol<long long>() { return 0; }


  //
  // UNIT TEST SERVICE FUNCTIONS
  //
  



  //
  // UNIT TESTS
  //

  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, doImport, LO, GO, Scalar, NO )
  {

    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();
    const RCP<const Comm<int> > comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();
    const int numProcs = comm->getSize();
    if(numProcs==1) return;


    // Prepare for verbose output, if applicable.
    //    Teuchos::EVerbosityLevel verbLevel = verbose ? VERB_EXTREME : VERB_NONE;
    Teuchos::EVerbosityLevel verbLevel = VERB_EXTREME;
    const bool doPrint = includesVerbLevel (verbLevel, VERB_EXTREME, true);
    if (doPrint) {
      out << "FEMultiVector unit test" << std::endl;
    }
    Teuchos::OSTab tab1 (out); // Add one tab level

    std::ostringstream err;
    int lclErr = 0;

    try {
      Teuchos::OSTab tab2 (out);
      const int num_local_elements = 3;

      // create Map
      RCP<const Tpetra::Map<LO, GO, NO> > map =
        rcp( new Tpetra::Map<LO,GO,NO>(INVALID, num_local_elements, 0, comm));

      // create CrsGraph object
      RCP<Tpetra::CrsGraph<LO, GO, NO> > graph =
             rcp (new Tpetra::CrsGraph<LO, GO, NO> (map, 3, Tpetra::DynamicProfile));

      // Create a simple tridiagonal source graph.
      Array<GO> entry(1);
      for (size_t i = 0; i < map->getNodeNumElements (); i++) {
        const GO globalrow = map->getGlobalElement (i);
        entry[0] = globalrow;
        graph->insertGlobalIndices (globalrow, entry());
        if (myRank!=0) {
          entry[0] = globalrow-1;
          graph->insertGlobalIndices (globalrow, entry());
        }
        if (myRank!=numProcs-1) {
          entry[0] = globalrow+1;
          graph->insertGlobalIndices (globalrow, entry());
        }       
      }
      graph->fillComplete();


      Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
      Scalar ONE = Teuchos::ScalarTraits<Scalar>::one();
      RCP<const Tpetra::Map<LO,GO,NO> > domainMap = graph->getDomainMap();
      RCP<const Tpetra::Map<LO,GO,NO> > columnMap = graph->getColMap();
      RCP<const Tpetra::Import<LO,GO,NO> > importer = graph->getImporter();
      Tpetra::MultiVector<Scalar,LO,GO,NO> Vdomain(domainMap,1), Vcolumn(columnMap,1);
      Tpetra::FEMultiVector<Scalar,LO,GO,NO> Vfe(domainMap,importer,1);
      size_t Ndomain = domainMap->getNodeNumElements();
      size_t Ncolumn = domainMap->getNodeNumElements();

      if(importer.is_null()) throw std::runtime_error("No valid importer");

      // 1) Test column -> domain (no off-proc addition)
      Vcolumn.putScalar(ZERO);
      for(size_t i=0; i<Ndomain; i++)
        Vcolumn.getDataNonConst(0)[i] = domainMap->getGlobalElement(i);
      Vdomain.doExport(Vcolumn,*importer,Tpetra::ADD);


      Vfe.beginFill();
      Vfe.putScalar(ZERO);
      for(size_t i=0; i<Ndomain; i++)
        Vfe.getDataNonConst(0)[i] = domainMap->getGlobalElement(i);
      Vfe.endFill();
      vector_check(Ndomain,Vfe,Vdomain);

      // 2) Test column -> domain (with off-proc addition)
      Vdomain.putScalar(ZERO);
      Vcolumn.putScalar(ONE);
      Vdomain.doExport(Vcolumn,*importer,Tpetra::ADD);

      Vfe.putScalar(ZERO);
      Vfe.beginFill();
      Vfe.putScalar(ONE);
      Vfe.endFill();
      vector_check(Ncolumn,Vfe,Vdomain);
    } catch (std::exception& e) {
      err << "Proc " << myRank << ": " << e.what () << std::endl;
      lclErr = 1;
    }


    int gblErr = 0;
    Teuchos::reduceAll<int, int> (*comm, Teuchos::REDUCE_MAX, lclErr, Teuchos::outArg (gblErr));
    TEST_EQUALITY_CONST( gblErr, 0 );
    if (gblErr != 0) {
      Tpetra::Details::gathervPrint (out, err.str (), *comm);
      out << "Above test failed; aborting further tests" << std::endl;
      return;
    }
  }





// ===============================================================================
  template<class Scalar, class LO, class GO, class Node>
  bool test_offsetview(Teuchos::FancyOStream &out, int femv_type) {
    bool success=true;
    typedef typename ScalarTraits<Scalar>::magnitudeType Mag;
    typedef Tpetra::FEMultiVector<Scalar,LO,GO,Node> FEMV;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    const global_size_t INVALID = OrdinalTraits<global_size_t>::invalid();
    const Scalar S0 = ScalarTraits<Scalar>::zero();
    const Mag M0 = ScalarTraits<Mag>::zero();
    const Mag tol = errorTolSlack * testingTol<Scalar>();
    // get a comm and node
    RCP<const Comm<int> > comm = getDefaultComm();
    // create a Map
    const size_t numLocal1 = 3;
    const size_t numLocal2 = 4;
    const size_t numLocal = numLocal1 + numLocal2;
    const size_t numVectors = 6;
    Array<size_t> even(tuple<size_t>(1,3,5));
    Array<size_t>  odd(tuple<size_t>(0,2,4));
    TEUCHOS_TEST_FOR_EXCEPTION( even.size() != odd.size(), std::logic_error, "Test setup assumption violated.");
    RCP<const Map<LO,GO,Node> > fullMap = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal,comm);
    RCP<const Map<LO,GO,Node> >    map1 = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal1,comm);
    RCP<const Map<LO,GO,Node> >    map2 = createContigMapWithNode<LO,GO,Node>(INVALID,numLocal2,comm);

    Array<GO> sharedList(numLocal-1), nonsharedList(numLocal-1);
    for(size_t i = 0; i<(size_t)sharedList.size(); i++) {
      sharedList[i]    = fullMap->getGlobalElement(i);
      nonsharedList[i] = fullMap->getGlobalElement(i+1);
    }
    RCP<const Map<LO,GO,Node> > sharedMap    = rcp(new Map<LO,GO,Node>(INVALID,sharedList(),0,comm));
    RCP<const Map<LO,GO,Node> > nonsharedMap = rcp(new Map<LO,GO,Node>(INVALID,nonsharedList(),0,comm));
    
    RCP<const Tpetra::Import<LO,GO,Node> > importer;
    // FEMV Type 1: Importer w/ shared memory
    if(femv_type == 1)
      importer = rcp(new Tpetra::Import<LO,GO,Node>(fullMap,sharedMap));
    // FEMV Type 2: Importer w/ non-shared memory
    else if(femv_type == 2)
      importer = rcp(new Tpetra::Import<LO,GO,Node>(fullMap,nonsharedMap));
    // FEMV Type 3: No importer
    
    // Use the FEMV in source mode
    RCP<FEMV> A = rcp(new FEMV(fullMap,importer,numVectors,false));
    A->beginFill(); A->endFill();
    
    {
      // config source multivector
      RCP<MV> A1 = A->offsetViewNonConst(map1, 0);
      RCP<MV> A2 = A->offsetViewNonConst(map2, numLocal1);
      TEST_EQUALITY( A1->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A1->getNumVectors(), numVectors );
      TEST_EQUALITY( A2->getNumVectors(), numVectors );
      Array<Mag>  A_befr(numVectors),
                 A1_befr(numVectors),
                 A2_befr(numVectors),
                  A_aft1(numVectors),
                 A1_aft1(numVectors),
                 A2_aft1(numVectors),
                  A_aft2(numVectors),
                 A1_aft2(numVectors),
                 A2_aft2(numVectors);
      // compute norms of A, A1 and A2
      A->randomize();
      A->norm2(A_befr());
      A1->norm2(A1_befr());
      A2->norm2(A2_befr());
      // set A1 = zeros, compute norms of A, A1 and A2
      A1->putScalar(S0);
      A->norm2(A_aft1());
      A1->norm2(A1_aft1());
      A2->norm2(A2_aft1());
      // set A2 = zeros, compute norms of A, A1 and A2
      A2->putScalar(S0);
      A->norm2(A_aft2());
      A1->norm2(A1_aft2());
      A2->norm2(A2_aft2());
      // change to A1 should not affect A2
      // change to A2 should not affect A1
      // change to A1 or A2 should change A
      // A should be zero after setting A1 to zero and A2 to zero
      for (size_t i=0; i<numVectors; ++i) {
        TEST_EQUALITY_CONST( A_aft1[i] < A_befr[i] + tol, true ); // shrunk as A1 = 0
        TEST_EQUALITY_CONST( A_aft2[i] < A_aft1[i] + tol, true ); // shrunk as A2 = 0
        TEST_EQUALITY_CONST( A_aft2[i] , M0 );                    // ... to zero
        TEST_EQUALITY_CONST( A1_aft1[i] , M0 );                   // was set to zero
        TEST_EQUALITY_CONST( A1_aft2[i] , M0 );                   // should not have been changed
        TEST_FLOATING_EQUALITY( A2_befr[i], A2_aft1[i], tol);     // should not have been changed
        TEST_EQUALITY_CONST( A2_aft2[i] , M0 );                   // was set to zero
      }
    }
    {
      // non-contig source multivector
      RCP<MV> A1e = A->subViewNonConst(even)->offsetViewNonConst(map1, 0);
      RCP<MV> A2e = A->subViewNonConst(even)->offsetViewNonConst(map2, numLocal1);
      RCP<MV> A1o = A->subViewNonConst(odd)->offsetViewNonConst(map1, 0);
      RCP<MV> A2o = A->subViewNonConst(odd)->offsetViewNonConst(map2, numLocal1);
      TEST_EQUALITY( A1e->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A1o->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2e->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A2o->getLocalLength(), numLocal2 );
      const size_t numSubVecs = (size_t)even.size();
      TEST_EQUALITY( A1e->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A2e->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A1o->getNumVectors(), numSubVecs );
      TEST_EQUALITY( A2o->getNumVectors(), numSubVecs );
      A->randomize();
      Array<Mag> b1(numSubVecs), b2(numSubVecs), b3(numSubVecs), bw(numVectors); // before putScalar(): unchanged 1, 2, 3; whole
      Array<Mag> a1(numSubVecs), a2(numSubVecs), a3(numSubVecs), aw(numVectors); // after putScalar(): ...
      Array<Mag> changed(numSubVecs), zeros(numSubVecs,M0);
      for (int i=0; i<4; ++i) {
        ArrayView<RCP<MV> > allMVs; // (changed,three unchanged)
        switch (i) {
        case 0:
          allMVs = tuple<RCP<MV> >(A1e,A2e,A1o,A2o); break;
        case 1:
          allMVs = tuple<RCP<MV> >(A2e,A1o,A2o,A1e); break;
        case 2:
          allMVs = tuple<RCP<MV> >(A1o,A2o,A1e,A2e); break;
        case 3:
          allMVs = tuple<RCP<MV> >(A2o,A1e,A2e,A1o); break;
        }
        allMVs[1]->norm2(b1()); allMVs[2]->norm2(b2()); allMVs[3]->norm2(b3());
        A->norm2(bw());
        allMVs[0]->putScalar(S0);
        allMVs[0]->norm2(changed());
        allMVs[1]->norm2(a1()); allMVs[2]->norm2(a2()); allMVs[3]->norm2(a3());
        A->norm2(aw());
        TEST_COMPARE_FLOATING_ARRAYS(b1,a1,tol);
        TEST_COMPARE_FLOATING_ARRAYS(b2,a2,tol);
        TEST_COMPARE_FLOATING_ARRAYS(b3,a3,tol);
        TEST_COMPARE_ARRAYS(changed(), zeros());
        for (size_t ii = 0; ii < numVectors; ++ii) {
          TEST_EQUALITY_CONST( aw[ii] < bw[ii] + tol, true ); // shrunk
        }
      }
    }
    {
      RCP<const MV> A1 = A->offsetView(map1, 0);
      RCP<const MV> A2 = A->offsetView(map2, numLocal1);
      TEST_EQUALITY( A1->getLocalLength(), numLocal1 );
      TEST_EQUALITY( A2->getLocalLength(), numLocal2 );
      TEST_EQUALITY( A1->getNumVectors(), numVectors );
      TEST_EQUALITY( A2->getNumVectors(), numVectors );
      Array<Mag>  A_bef(numVectors),
                 A1_bef(numVectors),
                 A2_bef(numVectors),
                  A_aft(numVectors),
                 A1_aft(numVectors),
                 A2_aft(numVectors);
      // compute norms of A, A1 and A2
      A->randomize();
      A->norm2(A_bef());
      A1->norm2(A1_bef());
      A2->norm2(A2_bef());
      A->putScalar(S0);
      A->norm2(A_aft());
      A1->norm2(A1_aft());
      A2->norm2(A2_aft());
      for (size_t i=0; i<numVectors; ++i) {
        TEST_EQUALITY_CONST( A_bef[i] < A1_bef[i] + A2_bef[i] + tol, true );
        TEST_EQUALITY_CONST( A_aft[i], S0 );
        TEST_EQUALITY_CONST( A1_aft[i], S0 );
        TEST_EQUALITY_CONST( A2_aft[i], S0 );
      }
    }
    return success;
  }




// ===============================================================================
 ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, OffsetView_Type1, LO , GO , Scalar , Node )
  {
    if(getDefaultComm()->getSize() ==1) return;
    success = test_offsetview<Scalar,LO,GO,Node>(out,1);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, OffsetView_Type2, LO , GO , Scalar , Node )
  {
    if(getDefaultComm()->getSize() ==1) return;
    success = test_offsetview<Scalar,LO,GO,Node>(out,2);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, OffsetView_Type3, LO , GO , Scalar , Node )
  {
    success = test_offsetview<Scalar,LO,GO,Node>(out,3);
  }



  
// ===============================================================================
  // This unit test exercises the following situation: Given a
  // Tpetra::MultiVector X, partition it into row blocks [X1; X2]
  // (Matlab notation) using offsetView (or offsetViewNonConst).  The
  // sum of the local number of rows in X1 and X2 equals the local
  // number of rows in X, but either X1 or X2 might have a zero number
  // of local rows.  We exercise each of the latter cases, in two
  // separate tests.  Repeat both cases for offsetView (const X1 and
  // X2) and offsetViewNonConst (nonconst X1 and X2).
  //
  // The most interesting thing this test exercises is that neither of
  // the above cases should throw exceptions.  This was not originally
  // true for the case where X2 has zero local rows.  Thanks to
  // Deaglan Halligan for pointing this out (on 23 Oct 2013).
  template<class Scalar, class LO, class GO, class Node>
  bool test_offsetviewzerolength(Teuchos::FancyOStream &out, int femv_type) {
    bool success = true;
    typedef Tpetra::global_size_t GST;
    typedef Tpetra::FEMultiVector<Scalar,LO,GO,Node> FEMV;
    typedef Tpetra::MultiVector<Scalar,LO,GO,Node> MV;
    typedef Tpetra::Map<LO, GO, Node> map_type;
    const GST INVALID = Teuchos::OrdinalTraits<GST>::invalid ();

    // Get a communicator and Kokkos node instance.
    RCP<const Comm<int> > comm = getDefaultComm ();

    // Create a Map with a nonzero number of entries on each process.
    const size_t numLocalEntries = 10;
    const GO indexBase = 0;
    RCP<const map_type> map =
      rcp (new map_type (INVALID, numLocalEntries, indexBase, comm));

    // Create a MultiVector X using that Map.  Give it some number of
    // columns (vectors) other than 1, just to exercise the most
    // general case.
    const size_t numVecs = 3;

    Array<GO> sharedList(numLocalEntries-1), nonsharedList(numLocalEntries-1);
    for(size_t i = 0; i<(size_t)sharedList.size(); i++) {
      sharedList[i]    = map->getGlobalElement(i);
      nonsharedList[i] = map->getGlobalElement(i+1);
    }
    RCP<const Map<LO,GO,Node> > sharedMap    = rcp(new Map<LO,GO,Node>(INVALID,sharedList(),0,comm));
    RCP<const Map<LO,GO,Node> > nonsharedMap = rcp(new Map<LO,GO,Node>(INVALID,nonsharedList(),0,comm));
    
    RCP<const Tpetra::Import<LO,GO,Node> > importer;
    // FEMV Type 1: Importer w/ shared memory
    if(femv_type == 1)
      importer = rcp(new Tpetra::Import<LO,GO,Node>(map,sharedMap));
    // FEMV Type 2: Importer w/ non-shared memory
    else if(femv_type == 2)
      importer = rcp(new Tpetra::Import<LO,GO,Node>(map,nonsharedMap));
    // FEMV Type 3: No importer
    
    // Use the FEMV in source mode
    RCP<FEMV> X = rcp(new FEMV(map,importer,numVecs));
    X->beginFill(); X->endFill();
    

    // Make sure that X has the right (local) dimensions.
    TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
    TEST_EQUALITY( X->getNumVectors (), numVecs );

    // Create a Map with zero entries on every process.
    RCP<const map_type> mapZero =
      rcp (new map_type (INVALID, 0, indexBase, comm));

    // Case 1: X1 has the same local number of rows as X, and X2 has
    // zero local rows.  Thus, X2 will be a zero-length view of X,
    // starting at the end of the local part of X (so the offset is
    // numLocalEntries).
    {
      RCP<const MV> X1;
      RCP<const MV> X2;
      try {
        X1 = X->offsetView (map, 0);
        X2 = X->offsetView (mapZero, numLocalEntries);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X2 has zero "
          "local rows." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY( X1->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X1->getNumVectors (), numVecs );
      TEST_EQUALITY_CONST( X2->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X2->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->template getLocalView<Kokkos::HostSpace> ();
      auto X1_local = X1->template getLocalView<Kokkos::HostSpace> ();
      auto X2_local = X2->template getLocalView<Kokkos::HostSpace> ();

      // Make sure the pointers match.  It doesn't really matter to
      // what X2_local points, as long as it has zero rows.
      TEST_EQUALITY( X1_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY( X1_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY_CONST( X2_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1 = Teuchos::null;
        X2 = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Nonconst version of Case 1.
    {
      RCP<MV> X1_nonconst;
      RCP<MV> X2_nonconst;
      try {
        X1_nonconst = X->offsetViewNonConst (map, 0);
        X2_nonconst = X->offsetViewNonConst (mapZero, numLocalEntries);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X2 has zero "
          "local rows, and X1 and X2 are nonconst." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY( X1_nonconst->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X1_nonconst->getNumVectors (), numVecs );
      TEST_EQUALITY_CONST( X2_nonconst->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X2_nonconst->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X1_local = X1_nonconst->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X2_local = X2_nonconst->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();

      // Make sure the pointers match.  It doesn't really matter to
      // what X2_local points, as long as it has zero rows.
      TEST_EQUALITY( X1_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY( X1_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY_CONST( X2_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1_nonconst = Teuchos::null;
        X2_nonconst = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Case 2: X1 has zero rows, and X2 has the same local number of
    // rows as X.  Thus, X1 will be a zero-length view of X, starting
    // at the beginning of the local part of X (so the offset is 0).
    {
      RCP<const MV> X1;
      RCP<const MV> X2;
      try {
        X1 = X->offsetView (mapZero, 0);
        X2 = X->offsetView (map, 0);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X1 has zero "
          "local rows." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY_CONST( X1->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X1->getNumVectors (), numVecs );
      TEST_EQUALITY( X2->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X2->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X1_local = X1->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X2_local = X2->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      // Make sure the pointers match.  It doesn't really matter to
      // what X1_local points, as long as it has zero rows.
      TEST_EQUALITY( X2_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY_CONST( X1_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY( X2_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1 = Teuchos::null;
        X2 = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }

    // Nonconst version of Case 2.
    {
      RCP<MV> X1_nonconst;
      RCP<MV> X2_nonconst;
      try {
        X1_nonconst = X->offsetViewNonConst (mapZero, 0);
        X2_nonconst = X->offsetViewNonConst (map, 0);
      } catch (...) {
        out << "The following case failed: X = [X1; X2] where X1 has zero "
          "local rows, and X1 and X2 are nonconst." << std::endl;
        throw;
      }
      // Make sure that offsetView() didn't change X's dimensions.
      TEST_EQUALITY( X->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X->getNumVectors (), numVecs );

      // Make sure that X1 and X2 have the right (local) dimensions.
      TEST_EQUALITY_CONST( X1_nonconst->getLocalLength (), static_cast<size_t> (0) );
      TEST_EQUALITY( X1_nonconst->getNumVectors (), numVecs );
      TEST_EQUALITY( X2_nonconst->getLocalLength (), numLocalEntries );
      TEST_EQUALITY( X2_nonconst->getNumVectors (), numVecs );

      // Make sure the pointers are the same, by extracting the
      // Kokkos::DualView objects.  Get the host pointer, just in case
      // MV allocation favors host space for initial allocations and
      // defers device allocations.

      auto X_local = X->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X1_local = X1_nonconst->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();
      auto X2_local = X2_nonconst->template getLocalView<typename MV::dual_view_type::t_host::memory_space> ();

      // Make sure the pointers match.  It doesn't really matter to
      // what X1_local points, as long as it has zero rows.
      TEST_EQUALITY( X2_local.data (), X_local.data () );

      // Make sure the local dimensions of X1 are correct.
      TEST_EQUALITY_CONST( X1_local.extent (0), static_cast<size_t> (0) );
      TEST_EQUALITY( X1_local.extent (1), X_local.extent (1) );

      // Make sure the local dimensions of X2 are correct.
      TEST_EQUALITY( X2_local.extent (0), X_local.extent (0) );
      TEST_EQUALITY( X2_local.extent (1), X_local.extent (1) );

      // Make sure that nothing bad happens on deallocation.
      try {
        X1_nonconst = Teuchos::null;
        X2_nonconst = Teuchos::null;
      } catch (...) {
        out << "Failed to deallocate X1 or X2." << std::endl;
        throw;
      }
    }
    return success;
  }

// ===============================================================================
 ////
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, OffsetViewZeroLength_Type1, LO , GO , Scalar , Node )
  {
    if(getDefaultComm()->getSize() ==1) return;
    success = test_offsetviewzerolength<Scalar,LO,GO,Node>(out,1);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, OffsetViewZeroLength_Type2, LO , GO , Scalar , Node )
  {
    if(getDefaultComm()->getSize() ==1) return;
    success = test_offsetviewzerolength<Scalar,LO,GO,Node>(out,2);
  }
  TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL( FEMultiVector, OffsetViewZeroLength_Type3, LO , GO , Scalar , Node )
  {
    success = test_offsetviewzerolength<Scalar,LO,GO,Node>(out,3);
  }




// ===============================================================================


#define UNIT_TEST_GROUP( SC, LO, GO, NO  )                         \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, doImport, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, OffsetView_Type1, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, OffsetView_Type2, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, OffsetView_Type3, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, OffsetViewZeroLength_Type1, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, OffsetViewZeroLength_Type2, LO, GO, SC, NO ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT( FEMultiVector, OffsetViewZeroLength_Type3, LO, GO, SC, NO )

 

  TPETRA_ETI_MANGLING_TYPEDEFS()



  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )
}
