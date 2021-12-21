#include <vector>
#include <iostream>
#include <numeric>
#include "Tpetra_ConfigDefs.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Teuchos_TimeMonitor.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "TpetraCore_ETIHelperMacros.h"


namespace {


template <typename LO, typename GO, typename NT>
Teuchos::RCP<const Tpetra::Map<LO,GO,NT> >
GetDofMap(std::vector<GO>& globalIDs) {
  Tpetra::global_size_t m_IGO = globalIDs.size();
  auto comm = Tpetra::getDefaultComm();
  return Teuchos::rcp(new Tpetra::Map<LO,GO,NT>(m_IGO,
                                      Teuchos::ArrayView<GO>(globalIDs),
                                      0, comm));
}

template <typename SC, typename LO, typename GO, typename NT>
Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NT> > 
GetDummyCrsMatrix(int numRow, int numCol) {
  assert(numRow>=numCol);
  std::vector<GO> rowIds(numRow);
  std::vector<GO> colIds(numCol);
  std::iota(rowIds.begin(), rowIds.end(), 0);
  std::iota(colIds.begin(), colIds.end(), 0);

  auto rowMap = GetDofMap<LO,GO,NT>(rowIds);
  auto colMap = GetDofMap<LO,GO,NT>(colIds);
  auto domainMap = GetDofMap<LO,GO,NT>(colIds);

  Teuchos::ArrayRCP<size_t> count(numRow);
  for (LO i=0; i<numRow; i++) count[i] = 1;

  Teuchos::RCP<Tpetra::CrsMatrix<SC,LO,GO,NT> > T =
    Teuchos::rcp(new Tpetra::CrsMatrix<SC,LO,GO,NT>(rowMap, colMap, count()));

  std::vector<LO> colIdxT(numRow);
  LO numRepeat = numRow-numCol;
  LO rowInd = 0;
  LO colInd = 0;
  for (LO ind=0; ind<numRepeat; ++ind) {
    colIdxT[rowInd++]=colInd;
    colIdxT[rowInd++]=colInd;
    ++colInd;
  }
  for (; rowInd<numRow; ++rowInd) {
    colIdxT[rowInd]=colInd;
    ++colInd;
  }
  Teuchos::Array<SC> dummyVal(1); dummyVal[0] = SC(1.);
  for (LO i = 0; i<numRow; ++i) {
    LO ii = rowIds[i];
    T->insertLocalValues(i, Teuchos::ArrayView<LO>(&colIdxT[ii], count[i]),
                         dummyVal());
  }
  T->fillComplete(domainMap, rowMap);
  return T;
}

TEUCHOS_UNIT_TEST_TEMPLATE_4_DECL(CrsMatrix, Bug9944, SC, LO, GO, NT)
{
  // This replicates the memory access of a matrix multiplication step,
  // which accounts for roughly 80% of the full analysis runtime in debug builds.
  const LO numRow = 3800;
  const LO numCol = 3000;
  const LO repeat = 4000;
  auto T = GetDummyCrsMatrix<SC, LO, GO, NT>(numRow, numCol);

#ifdef TPETRA_ENABLE_DEPRECATED_CODE
  {
    // Loop over getLocalRowView using deprecated Teuchos::ArrayView interface
    // This use case was provided by #9944
    Teuchos::ArrayView<const LO> Indices;
    Teuchos::ArrayView<const SC> Values;
    GO sumi = 0.;
    SC sumv = 0.;

    auto timer = 
         Teuchos::TimeMonitor::getNewTimer("getLocalRowView "
                                           "Teuchos::ArrayView");
    {
      Teuchos::TimeMonitor tt(*timer);
      for (LO i=0; i<repeat; ++i) {
        for (LO row=0; row<numRow; ++row) {
          T->getLocalRowView(row, Indices, Values); // This line slows down more than 100x
          sumi += Indices[0];
          sumv += Values[0];
        }
      }
    }
    std::cout << "sumi = " << sumi << std::endl;
    std::cout << "sumv = " << sumv << std::endl;
  }
#endif
  {
    // Same loop using the new Kokkos::View interface
    typename Tpetra::CrsMatrix<SC,LO,GO,NT>::local_inds_host_view_type Indices;
    typename Tpetra::CrsMatrix<SC,LO,GO,NT>::values_host_view_type Values;

    GO sumi = 0.;
    typename Tpetra::CrsMatrix<SC,LO,GO,NT>::impl_scalar_type sumv = 0.;

    auto timer = 
         Teuchos::TimeMonitor::getNewTimer("getLocalRowView Kokkos::View");
    {
      Teuchos::TimeMonitor tt(*timer);
      for (LO i=0; i<repeat; ++i) {
        for (LO row=0; row<numRow; ++row) {
          T->getLocalRowView(row, Indices, Values); 
          sumi += Indices[0];
          sumv += Values[0];
        }
      }
    }
    std::cout << "sumi = " << sumi << std::endl;
    std::cout << "sumv = " << sumv << std::endl;
  }
  {
    // A more efficient implementation:
    // get the views once, then index into them
    GO sumi = 0;
    typename Tpetra::CrsMatrix<SC,LO,GO,NT>::impl_scalar_type sumv = 0.;

    auto lclMatrix = T->getLocalMatrixHost();
    auto allIndices = lclMatrix.graph.entries;
    auto allValues = lclMatrix.values;
    auto timer = 
         Teuchos::TimeMonitor::getNewTimer("getLocalMatrixHost "
                                           "Kokkos::subview");
    {
      Teuchos::TimeMonitor tt(*timer);
      for (LO i=0; i<repeat; ++i) {
        for (LO row=0; row<numRow; ++row) {
          auto offset = lclMatrix.graph.row_map[row];
          auto numEntries = lclMatrix.graph.row_map[row+1] - offset;
          auto Indices =
               Kokkos::subview(allIndices,
                               Kokkos::pair<int,int>(offset,offset+numEntries));
          auto Values =
               Kokkos::subview(allValues,
                               Kokkos::pair<int,int>(offset,offset+numEntries));
          sumi += Indices[0];
          sumv += Values[0];
        }
      }
    }
    std::cout << "sumi = " << sumi << std::endl;
    std::cout << "sumv = " << sumv << std::endl;
  }

  Teuchos::TimeMonitor::summarize();
  Teuchos::TimeMonitor::zeroOutTimers();
  std::cout << "TEST PASSED" << std::endl;
}

#define UNIT_TEST_GROUP( SC, LO, GO, NT ) \
  TEUCHOS_UNIT_TEST_TEMPLATE_4_INSTANT(CrsMatrix, Bug9944, SC, LO, GO, NT)

  TPETRA_ETI_MANGLING_TYPEDEFS()

  TPETRA_INSTANTIATE_TESTMV( UNIT_TEST_GROUP )

} // namespace
