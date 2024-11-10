// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Creates vectors with different maps; tests results of export into them
// Documents behavior of Tpetra::ADD in Tpetra::MultiVector for many common
// (and a few less common) use cases.
// Analogous tests for Epetra are in epetra/test/MultiVector/Bug7758.cpp

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_Array.hpp"

#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"

namespace {


//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7758, DefaultToDefault, Scalar,LO,GO,Node)
{
  // This case demonstrates that owned entries shared between the source and
  // target map are copied from the source vector into the target (during
  // copyAndPermute).  Each entry of the resulting target vector has value
  // srcScalar.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar tgtScalar = 100. * (me+1);
  const Scalar srcScalar = 2.;
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos

  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // Create vectors; see what the result is with CombineMode=ADD

  vector_t defaultVecTgt(defaultMap);
  defaultVecTgt.putScalar(tgtScalar);

  vector_t defaultVecSrc(defaultMap);
  defaultVecSrc.putScalar(srcScalar);

  // Export Default-to-default:  should be a copy of src to tgt

  Tpetra::Export<LO,GO,Node> defaultToDefault(defaultMap, defaultMap);
  defaultVecTgt.doExport(defaultVecSrc, defaultToDefault, Tpetra::ADD);

  // Check result; all vector entries should be srcScalar
  auto data = defaultVecTgt.getLocalViewHost(Tpetra::Access::ReadOnly);

  for (size_t i = 0; i < defaultVecTgt.getLocalLength(); i++)
    if (data(i,0) != srcScalar) ierr++;
  if (ierr > 0) 
    std::cout << "TEST FAILED:  DEFAULT-TO-DEFAULT TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7758, CyclicToDefault, Scalar,LO,GO,Node)
{
  // This case demonstrates that owned entries shared between the source and
  // target map are copied from the source vector into the target (during
  // copyAndPermute).  Owned entries that are not in the source map
  // are NOT reset; their initial values persist.  Then received shared 
  // entries are added to the owned entries.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar tgtScalar = 100. * (me+1);
  const Scalar srcScalar = 2.;
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Default one-to-one linear block map in Trilinos

  Teuchos::RCP<const map_t> defaultMap = 
           rcp(new map_t(nGlobalEntries, 0, comm));

  // One-to-one cyclic map:  deal out entries like cards

  int nMyEntries = 0;
  for (size_t i = 0; i < nGlobalEntries; i++) {
    if (int(i % np) == me) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> cyclicMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Create vectors; see what the result is with CombineMode=ADD

  vector_t defaultVecTgt(defaultMap);
  defaultVecTgt.putScalar(tgtScalar);

  vector_t cyclicVecSrc(cyclicMap);
  cyclicVecSrc.putScalar(srcScalar);

  // Export Cyclic-to-default

  Tpetra::Export<LO,GO,Node> cyclicToDefault(cyclicMap, defaultMap);
  defaultVecTgt.doExport(cyclicVecSrc, cyclicToDefault, Tpetra::ADD);

  // Check result

  auto invalid = Teuchos::OrdinalTraits<LO>::invalid();
  auto data = defaultVecTgt.getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i = 0; i < defaultVecTgt.getLocalLength(); i++)
    if (cyclicMap->getLocalElement(defaultMap->getGlobalElement(i)) != invalid){
      // element is in both cyclic (source) and default (target) map;
      // initial target value was overwritten
      if (data(i,0) != srcScalar) ierr++;
    }
    else {
      // element is in only default (target) map;
      // initial target value was not overwritten
      if (data(i,0) != tgtScalar + srcScalar) ierr++;
    }
  if (ierr > 0) 
    std::cout << "TEST FAILED:  CYCLIC-TO-DEFAULT TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7758, OverlapToDefault, Scalar,LO,GO,Node)
{
  // This case demonstrates that owned entries shared between the source and
  // target map are copied from the source vector into the target (during
  // copyAndPermute).  Owned entries that are not in the source map
  // are NOT reset; their initial values persist.  Then received shared 
  // entries are added to the owned entries.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

    using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;

    const size_t nGlobalEntries = 8 * np;
    const Scalar tgtScalar = 100. * (me+1);
    const Scalar srcScalar = 2.;
    Teuchos::Array<GO> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    Teuchos::RCP<const map_t> defaultMap = 
             rcp(new map_t(nGlobalEntries, 0, comm));

    // Overlap map; some entries are stored on two procs
    int nMyEntries = 0;
    for (size_t i = 0; i < defaultMap->getLocalNumElements()/2; i++) {
      myEntries[nMyEntries++] = defaultMap->getGlobalElement(i);
    }
    for (size_t i = 0; i < defaultMap->getLocalNumElements(); i++) {
      myEntries[nMyEntries++] =
        (defaultMap->getMaxGlobalIndex() + 1 + i) % nGlobalEntries;
    }
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> overlapMap = 
             rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));
  
    // Create vectors; see what the result is with CombineMode=ADD

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.putScalar(tgtScalar);

    vector_t overlapVecSrc(overlapMap);
    overlapVecSrc.putScalar(srcScalar);

    // Export Overlap-to-default

    Tpetra::Export<LO,GO,Node> overlapToDefault(overlapMap, defaultMap);
    defaultVecTgt.doExport(overlapVecSrc, overlapToDefault, Tpetra::ADD);

    auto data = defaultVecTgt.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (size_t i = 0; i < defaultVecTgt.getLocalLength()/2; i++) {
      // overlapped; initial target values were overwritten
      if (data(i,0) != srcScalar + srcScalar) ierr++;  
    }
    for (size_t i = defaultVecTgt.getLocalLength()/2;
             i < defaultVecTgt.getLocalLength(); i++) {
      // not overlapped; initial target values were not overwritten
      if (data(i,0) != tgtScalar + srcScalar) ierr++;  
    }
    if (ierr > 0) 
      std::cout << "TEST FAILED:  OVERLAP-TO-DEFAULT TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7758, OddEvenToSerial, Scalar,LO,GO,Node)
{
  // Test case showing behavior when target map is all on processor zero.
  // In the source map, even numbered entries are on even numbered processors;
  // odd numbered entreis are on odd numbered processors.
  // In copyAndPermute, even numbered entries are copied from processor zero's
  // source vector to the target vector, and odd numbered entries are unchanged.
  // Then received values are added to the target vector.  The result is that
  // odd entries include the initial target values in their sum, while the 
  // even entries are not included.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

  using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
  using map_t = Tpetra::Map<LO,GO,Node>;

  const size_t nGlobalEntries = 8 * np;
  const Scalar tgtScalar = 100. * (me+1);
  const Scalar srcScalar = 2.;
  Teuchos::Array<GO> myEntries(nGlobalEntries); 

  // Odd entries given to odd procs; even entries given to even procs
  int nMyEntries = 0;
  for (size_t i = 0; i < nGlobalEntries; i++) {
    if (int(i % 2) == (me % 2)) {
      myEntries[nMyEntries++] = i;
    }
  }

  Tpetra::global_size_t dummy =
          Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  Teuchos::RCP<const map_t> oddEvenMap = 
           rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));

  // Map with all entries on one processor

  dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
  size_t nSerialEntries = (me == 0 ? nGlobalEntries : 0);
  Teuchos::RCP<const map_t> serialMap = 
           rcp(new map_t(dummy, nSerialEntries, 0, comm));

  // Create vectors; see what the result is with CombineMode=ADD

  vector_t oddEvenVecSrc(oddEvenMap);
  oddEvenVecSrc.putScalar(srcScalar);

  vector_t serialVecTgt(serialMap);
  serialVecTgt.putScalar(tgtScalar);

  // Export oddEven-to-serial

  Tpetra::Export<LO,GO,Node> oddEvenToSerial(oddEvenMap, serialMap);
  serialVecTgt.doExport(oddEvenVecSrc, oddEvenToSerial, Tpetra::ADD);

  // Check result

  auto data = serialVecTgt.getLocalViewHost(Tpetra::Access::ReadOnly);
  for (size_t i = 0; i < serialVecTgt.getLocalLength(); i++) {
    Scalar nCopies = Scalar(((np+1) / 2) - ((i % 2 == 1) && (np % 2 == 1)));
    if (data(i,0) != (i%2 ? tgtScalar : Scalar(0)) + srcScalar * nCopies)
      ierr++;
  }
  if (ierr > 0) 
    std::cout << "TEST FAILED:  ODDEVEN-TO-SERIAL TEST HAD " << ierr 
              << " FAILURES ON RANK " << me << std::endl;

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7758, SupersetToDefault, Scalar,LO,GO,Node)
{
  // This use case is similar to matrix assembly case in which user 
  // has a map of owned entries and a map of owned+shared entries, with the
  // owned+shared map being a superset of the owned map.  In this case, 
  // the owned values in the owned+shared vector are copied into the owned
  // vector in copyAndPermute; then received shared entries are 
  // added to the owned entries' values.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

    using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;

    const size_t nGlobalEntries = 8 * np;
    const Scalar tgtScalar = 100. * (me+1);
    const Scalar srcScalar = 2.;
    Teuchos::Array<GO> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    Teuchos::RCP<const map_t> defaultMap = 
             rcp(new map_t(nGlobalEntries, 0, comm));

    // Superset map; some entries are stored on two procs
    int nMyEntries = 0;
    for (size_t i = 0; i < defaultMap->getLocalNumElements(); i++) {
      myEntries[nMyEntries++] = defaultMap->getGlobalElement(i);
    }
    for (size_t i = 0; i < defaultMap->getLocalNumElements()/2; i++) {
      myEntries[nMyEntries++] =
        (defaultMap->getMaxGlobalIndex() + 1 + i) % nGlobalEntries;
    }
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> supersetMap =
             rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));
  
    // Create vectors; see what the result is with CombineMode=ADD

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.putScalar(tgtScalar);

    vector_t supersetVecSrc(supersetMap);
    supersetVecSrc.putScalar(srcScalar);

    // Export Superset-to-default

    Tpetra::Export<LO,GO,Node> supersetToDefault(supersetMap, defaultMap);
    defaultVecTgt.doExport(supersetVecSrc, supersetToDefault, Tpetra::ADD);

    auto data = defaultVecTgt.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (size_t i = 0; i < defaultVecTgt.getLocalLength()/2; i++)
      if (data(i,0) != srcScalar+srcScalar) ierr++;  // overlapped
    for (size_t i = defaultVecTgt.getLocalLength()/2;
             i < defaultVecTgt.getLocalLength(); i++)
      if (data(i,0) != srcScalar) ierr++;  // not overlapped
    if (ierr > 0) 
      std::cout << "TEST FAILED:  SUPERSET-TO-DEFAULT TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

//////////////////////////////////////////////////////////////////////////////
TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(Bug7758, NoSamesToDefault, Scalar,LO,GO,Node)
{
  // This use case is similar to matrix assembly case in which user 
  // has a map of owned entries and a map of shared entries, with no
  // overlap between the maps.  In this case, received shared entries are 
  // added to the owned entries' values, as copyAndPermute is never called.
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int np = comm->getSize();
  int ierr = 0;

  if (np > 1) {  // Need more than one proc to avoid duplicate entries in maps
    Teuchos::FancyOStream foo(Teuchos::rcp(&std::cout,false));

    using vector_t = Tpetra::Vector<Scalar,LO,GO,Node>;
    using map_t = Tpetra::Map<LO,GO,Node>;

    const size_t nGlobalEntries = 8 * np;
    const Scalar tgtScalar = 100. * (me+1);
    const Scalar srcScalar = 2.;
    Teuchos::Array<GO> myEntries(nGlobalEntries); 

    // Default one-to-one linear block map in Trilinos

    Teuchos::RCP<const map_t> defaultMap = 
             rcp(new map_t(nGlobalEntries, 0, comm));

    // Map with no sames or permutes
    int nMyEntries = 0;
    for (size_t i = 0; i < defaultMap->getLocalNumElements(); i++) {
      myEntries[nMyEntries++] =
        (defaultMap->getMaxGlobalIndex() + 1 + i) % nGlobalEntries;
    }
  
    Tpetra::global_size_t dummy = 
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    Teuchos::RCP<const map_t> noSamesMap = 
             rcp(new map_t(dummy, myEntries(0,nMyEntries), 0, comm));
  
    // Create vectors; see what the result is with CombineMode=ADD

    vector_t defaultVecTgt(defaultMap);
    defaultVecTgt.putScalar(tgtScalar);

    vector_t noSamesVecSrc(noSamesMap);
    noSamesVecSrc.putScalar(srcScalar);

    // Export noSames-to-default

    Tpetra::Export<LO,GO,Node> noSamesToDefault(noSamesMap, defaultMap);
    defaultVecTgt.doExport(noSamesVecSrc, noSamesToDefault, Tpetra::ADD);

    auto data = defaultVecTgt.getLocalViewHost(Tpetra::Access::ReadOnly);
    for (size_t i = 0; i < defaultVecTgt.getLocalLength(); i++)
      if (data(i,0) != tgtScalar + srcScalar) ierr++;  
    if (ierr > 0) 
      std::cout << "TEST FAILED:  NOSAMES-TO-DEFAULT TEST HAD " << ierr 
                << " FAILURES ON RANK " << me << std::endl;
  }

  int gerr;
  Teuchos::reduceAll<int,int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  TEST_ASSERT(gerr == 0);
}

#define UNIT_TEST_GROUP( SCALAR, LO, GO, NODE ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7758, DefaultToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7758, CyclicToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7758, OverlapToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7758, OddEvenToSerial, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7758, SupersetToDefault, SCALAR, LO, GO, NODE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(Bug7758, NoSamesToDefault, SCALAR, LO, GO, NODE)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace (anonymous)

