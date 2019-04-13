#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <string>
#include <sstream>
#include <iostream>

// Struct representing a vertex label.
// We define our own "addition" for these labels.
// Later, we'll create a Tpetra::FEMultiVector of these labels.
class vtxLabel {
public:
  int label;
  // Constructors
  vtxLabel(int y) { label = y; }
  vtxLabel() { label = 0; }

  // vtxLabel assignment
  vtxLabel& operator=(const vtxLabel& other)  { 
    label = other.label; return *this;
  }

  // Added to satisfy Kokkos compilation (kokkos_atomic_fetch)
  vtxLabel& operator=(volatile const vtxLabel& other)  { 
    label = other.label; return *this;
  }

  // int assignment
  vtxLabel& operator=(const int& other) { 
    label = other; return *this;
  }
  // += overload
  vtxLabel& operator+=(const vtxLabel& other) { 
    label = -label - other.label; return *this;
  }
  // addition overload
  friend vtxLabel operator+(const vtxLabel& lhs, const vtxLabel& rhs) {
    vtxLabel result(-lhs.label + -rhs.label);
    return result;
  }
  // vtxLabel equality overload
  friend bool operator==(const vtxLabel& lhs, const vtxLabel& rhs) {
    return (lhs.label == rhs.label);
  }
  // int equality overload
  friend bool operator==(const vtxLabel& lhs, const int& rhs) {
    return (lhs.label == rhs);
  }
  // output stream overload
  friend std::ostream& operator<<(std::ostream& os, const vtxLabel& a) {
    os << a.label; return os;
  }
};

//  Unit test for the vtxLabel struct
//  Make sure vtxLabel's overloaded operators compile and work as expected.
int vtxLabelUnitTest()
{
  int ierr = 0;

  const int A = 10;
  const int B = 5;
  const int E = -15;

  vtxLabel a(A), b(B);
  vtxLabel expect(E);

  vtxLabel d(A);
  if (!(d == A)) {
    std::cout << "unitTest Error:  integer equality not as expected " << d
              << "!=" << A << std::endl;   
    ierr++;
  }
  else std::cout << "unitTest:  integer equality OK" << std::endl;

  if (!(d == a)) {
    std::cout << "unitTest Error:  vtxLabel equality not as expected " << d
              << "!=" << a << std::endl;   
    ierr++;
  }
  else std::cout << "unitTest:  vtxLabel equality OK" << std::endl;

  d += b;
  if (!(d == expect)) {
    std::cout << "unitTest Error:  += not as expected " << d << "!="
              << expect << std::endl;   
    ierr++;
  }
  else std::cout << "unitTest:  += OK" << std::endl;

  d = a + b;
  if (!(d == expect)) {
    std::cout << "unitTest Error:  addition not as expected " << a << "+" << b 
              << "=" << d << "!=" << expect << std::endl;   
    ierr++;
  }
  else std::cout << "unitTest:  addition OK" << std::endl;
  
  return ierr;
}

/////////////////////////////////////////////////////////////////////////
// ArithTraits -- arithmetic traits needed for struct vtxLabel
// Needed so that Tpetra compiles.
// Not all functions were needed; this is a subset of ArithTraits' traits.
// Modified from kokkos-kernels/src/Kokkos_ArithTraits.hpp's 
// <int> specialization

namespace Kokkos {
  namespace Details {

    template<>
    class ArithTraits<vtxLabel> {  // specialized for vtxLabel struct
    public:
      typedef vtxLabel val_type;
      typedef int mag_type;
    
      static const bool is_specialized = true;
      static const bool is_signed = true;
      static const bool is_integer = true;
      static const bool is_exact = true;
      static const bool is_complex = false;
    
      static KOKKOS_FORCEINLINE_FUNCTION bool isInf(const val_type &) {
        return false;
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isNan(const val_type &) {
        return false;
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type abs(const val_type &x) {
        return (x.label >= 0 ? x.label : -(x.label));
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type zero() { return 0; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type one() { return 1; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type min() { return INT_MIN; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type max() { return INT_MAX; }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type nan() { return -1; }
    
      // Backwards compatibility with Teuchos::ScalarTraits.
      typedef mag_type magnitudeType;
      static const bool isComplex = false;
      static const bool isOrdinal = true;
      static const bool isComparable = true;
      static const bool hasMachineParameters = false;
      static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude(
        const val_type &x) 
      {
        return abs(x);
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf(const val_type &) {
        return false;
      }
      static std::string name() { return "vtxLabel"; }
    };
  }
}

/////////////////////////////////////////////////////////////////////////////
// Teuchos::SerializationTraits are needed to copy vtxLabels into MPI buffers
// Because sizeof(vtxLabel) works for struct vtxLabel, we'll use a 
// provided serialization of vtxLabel into char*.
namespace Teuchos {
template<typename Ordinal>
struct SerializationTraits<Ordinal, vtxLabel> :
       public DirectSerializationTraits<Ordinal, vtxLabel>
{};
}

/////////////////////////////////////////////////////////////////////////////
// Class to test 
// (1) how FEMultiVector works and 
// (2) whether it works with user-defined scalar type (vtxLabel).
// Each processor owns ten vertices and copies five from another processor.
// Each processor contributes values to all fifteen vertices it stores.
// Tpetra::FEMultiVector pushes the copies to their owning processors, and 
// ADDs them to the owned values.

class FEMultiVectorTest {
public:

  typedef Tpetra::Map<> map_t;
  typedef map_t::local_ordinal_type lno_t;
  typedef map_t::global_ordinal_type gno_t;

  // Constructor assigns vertices to processors and builds maps with and 
  // without copies
  FEMultiVectorTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm_) :
    me(comm_->getRank()), np(comm_->getSize()),
    nLocalOwned(10), nLocalCopy( np > 1 ? 5 : 0), 
    nVec(2), comm(comm_)
  {
    // Each rank has 15 IDs, the last five of which overlap with the next rank.
    // (IDs and owning processors wrap-around from processor np-1 to 0.)
    const Tpetra::global_size_t nGlobal = np * nLocalOwned;
    lno_t offset = me * nLocalOwned;

    Teuchos::Array<gno_t> gids(nLocalOwned + nLocalCopy);
    for (lno_t i = 0 ; i < nLocalOwned + nLocalCopy; i++)
      gids[i] = static_cast<gno_t> (offset + i) % nGlobal;

    // Create Map of owned + copies (a.k.a. overlap map); analagous to ColumnMap
    Tpetra::global_size_t dummy =
            Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    mapWithCopies = rcp(new map_t(dummy, gids(), 0, comm));

    // Create Map of owned only (a.k.a. one-to-one map); analagous to RowMap
    mapOwned = rcp(new map_t(dummy, gids(0, nLocalOwned), 0, comm));

    // Print the entries of each map
    std::cout << me << " MAP WITH COPIES ("
                    << mapWithCopies->getGlobalNumElements() << "):  ";
    lno_t nlocal = lno_t(mapWithCopies->getNodeNumElements());
    for (lno_t i = 0; i < nlocal; i++)
      std::cout << mapWithCopies->getGlobalElement(i) << " ";
    std::cout << std::endl;

    std::cout << me << " ONE TO ONE MAP  ("
                    << mapOwned->getGlobalNumElements() << "):  ";
    nlocal = lno_t(mapOwned->getNodeNumElements());
    for (lno_t i = 0; i < nlocal; i++)
      std::cout << mapOwned->getGlobalElement(i) << " ";
    std::cout << std::endl;
  }

  // Create a Tpetra::FEMultiVector of type femv_t with two vectors
  // FEMultiVector length is number of owned vertices + number of copies of 
  //               off-processor vertices (as indicated by mapWithCopies)
  // Fill first vector with GIDs of vertices
  // Fill second vector with processor rank (me)
  // Perform communication to add copies' contributions to their owned vertices
  template <typename femv_t>
  Teuchos::RCP<femv_t> getFEMV()
  {
    // Create FEMultiVector
    typedef Tpetra::Import<lno_t, gno_t> import_t;
    Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                       mapWithCopies));
    Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, importer,
                                               nVec, true));
  
    // Store contributions to owned vertices and copies of off-processor 
    // vertices
    femv->beginFill();
    for (lno_t i = 0; i < nLocalOwned + nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      femv->replaceGlobalValue(gid, 0, gid);
      femv->replaceGlobalValue(gid, 1, me);
      // Other useful calls for filling may be:
      //   replaceLocalValue -- may be more efficient if storing edges using
      //                        LIDs rather than GIDs
      //   sumInto{Global,Local}Value -- may be able to exploit the overloaded
      //                                 addition to propagate local labels
    }
    femv->endFill(); // communication and summation of copies is done in endFill
                     // "Sums" contributions from all copies into the owned
                     // vertex on the owning processor

    printFEMV(femv, "AfterFill");

    // Now update copied vertices from their owners
    femv->doOwnedToOwnedPlusShared(Tpetra::REPLACE);

    // For the ice-sheet problem, 
    // would a Tpetra::ADD work here (consistently update the non-owned vertex's
    // values and add the non-owned vertex to the frontier if it has 
    // significant change)?  Consider changing REPLACE to ADD.
    // Sends the owner's result for each vertex to all of the copies of the
    // vertex.  After this, all copies will have the same value as the owned
    // vertex.

    printFEMV(femv, "AfterReplace");

    return femv;
  }

  template <typename femv_t>
  void printFEMV(femv_t &femv, const char *msg) 
  {
    for (int v = 0; v < nVec; v++) {
      std::cout << me << " OWNED " << msg << " FEMV[" << v << "] Owned: ";
      auto value = femv->getData(v);
      for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
      std::cout << std::endl;
    }
    femv->switchActiveMultiVector();  // Needed to print copies
    for (int v = 0; v < nVec; v++) {
      std::cout << me << " WITHCOPIES " << msg << " FEMV[" << v << "] Owned: ";
      auto value = femv->getData(v);
      for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
      std::cout << " Copies: ";
      for (lno_t i = nLocalOwned; i < nLocalOwned+nLocalCopy; i++) 
        std::cout << value[i] << " ";
      std::cout << std::endl;
    }
    femv->switchActiveMultiVector();  // Restore state
  }

  // Test using scalar_t=int to exercise push capability of FEMultiVector
  int intTest();

  // Test using scalar_t=custom data type (vtxLabel)
  int vtxLabelTest();

private:
  int me;           // my processor rank
  int np;           // number of processors
  int nLocalOwned;  // number of vertices owned by this processor
  int nLocalCopy;   // number of copies of off-processor vertices on this proc
  int nVec;         // number of vectors in multivector

  Teuchos::RCP<const Teuchos::Comm<int> > comm;  // MPI communicator

  Teuchos::RCP<const map_t> mapWithCopies;  // Tpetra::Map including owned
                                            // vertices and copies
  Teuchos::RCP<const map_t> mapOwned;       // Tpetra::Map including only owned
};


// Test using scalar_t=int to exercise push capability of FEMultiVector
// Create FEMultiVector and verify correct results with integer addition
int FEMultiVectorTest::intTest()
{
  typedef int scalar_t;
  typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;
  int ierr = 0;

  Teuchos::RCP<femv_t> femv = getFEMV<femv_t>();

  // Check results:  after ADD in endFill, 
  // -  overlapping entries of vec 0 should be 2 * gid
  //    nonoverlapping entries of vec 0 should be gid
  // -  overlapping entries of vec 1 should be me + (np + me-1) % np;
  //    nonoverlapping entries of vec 1 should be me

  // Check vector 0
  {
    auto value = femv->getData(0);
    for (lno_t i = 0; i < nLocalCopy; i++){
      gno_t gid = femv->getMap()->getGlobalElement(i);
      if (value[i] != 2*gid) {
        std::cout << me << " Error in vec 0 overlap: gid=" << gid 
                        << " value= " << value[i] << " should be " << 2*gid
                        << std::endl;
        ierr++;
      }
    }
    for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
      gno_t gid = femv->getMap()->getGlobalElement(i);
      if (value[i] != gid) {
        std::cout << me << " Error in vec 0:  gid=" << gid
                        << " value= " << value[i] << " should be " << gid
                        << std::endl;
        ierr++;
      }
    }

    femv->switchActiveMultiVector();  // Needed to view copies
    value = femv->getData(0);
    for (lno_t i = nLocalOwned; i < nLocalCopy + nLocalOwned; i++) {
      gno_t gid = femv->getMap()->getGlobalElement(i);
      if (value[i] != 2*gid) {
        std::cout << me << " Error in vec 0 copies:  gid=" << gid
                        << " value= " << value[i] << " should be " << 2*gid
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Restore state 
  }

  // Check vector 1
  {
    auto value = femv->getData(1);
    int tmp = me + (np + me - 1) % np;
    for (lno_t i = 0; i < nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (value[i] != tmp) {
        std::cout << me << " Error in vec 1 overlap:  gid=" << gid
                        << " value= " << value[i] << " should be " << tmp
                        << std::endl;
        ierr++;
      }
    }
    for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (value[i] != me) {
        std::cout << me << " Error in vec 1:  gid=" << gid
                        << " value= " << value[i] << " should be " << me
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Needed to view copies
    value = femv->getData(1);
    tmp = me + (me + 1) % np;
    for (lno_t i = nLocalOwned; i < nLocalCopy + nLocalOwned; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (value[i] != tmp)  {
        std::cout << me << " Error in vec 1:  gid=" << gid
                        << " value= " << value[i] << " should be " << tmp
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Restore State
  }
  return ierr;
}

// Test using scalar_t = custom data type vtxLabel
// Create FEMultiVector and verify correct results with vtxlabel "addition" 
int FEMultiVectorTest::vtxLabelTest()
{
  typedef vtxLabel scalar_t;
  typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;
  int ierr = 0;

  Teuchos::RCP<femv_t> femv = getFEMV<femv_t>();

  // Check results:  after ADD in endFill, 
  // -  overlapping entries of vec 0 should be -2 * gid
  //    nonoverlapping entries of vec 0 should be gid
  // -  overlapping entries of vec 1 should be -me - ((np + me-1) % np)
  //    nonoverlapping entries of vec 1 should be me
  
  // Check vector 0
  {
    auto value = femv->getData(0);
    for (lno_t i = 0; i < nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (!(value[i] == -2*gid)) {
        std::cout << me << " Error in vec 0 overlap:  gid=" << gid
                        << " value= " << value[i] << " should be " << -2*gid
                        << std::endl;
        ierr++;
      }
    }
    for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (!(value[i] == gid)) {
        std::cout << me << " Error in vec 0:  gid=" << gid
                        << " value= " << value[i] << " should be " << gid
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Needed to view copies
    value = femv->getData(0);
    for (lno_t i = nLocalOwned; i < nLocalCopy + nLocalOwned; i++) {
      gno_t gid = femv->getMap()->getGlobalElement(i);
      if (!(value[i] == -2*gid)) {
        std::cout << me << " Error in vec 0 copies:  gid=" << gid
                        << " value= " << value[i] << " should be " << -2*gid
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Restore state 
  }

  // Check vector 1
  {
    auto value = femv->getData(1);
    int tmp = -me - (np + me - 1) % np;
    for (lno_t i = 0; i < nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (!(value[i] == tmp)) {
        std::cout << me << " Error in vec 1 overlap:  gid=" << gid
                        << " value= " << value[i] << " should be " << tmp
                        << std::endl;
        ierr++;
      }
    }
    for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (!(value[i] == me)) {
        std::cout << me << " Error in vec 1:  gid=" << gid
                        << " value= " << value[i] << " should be " << me
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Needed to view copies
    value = femv->getData(1);
    tmp = -me - (me + 1) % np;
    for (lno_t i = nLocalOwned; i < nLocalCopy + nLocalOwned; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (!(value[i] == tmp))  {
        std::cout << me << " Error in vec 1:  gid=" << gid
                        << " value= " << value[i] << " should be " << tmp
                        << std::endl;
        ierr++;
      }
    }
    femv->switchActiveMultiVector();  // Restore State
  }
  return ierr;
}

/////////////////////////////////////////////////////////////////////

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();
  int ierr = 0;

  // Test that the vtxLabel struct correctly overloads operators
  if (me == 0) 
    ierr += vtxLabelUnitTest();

  // Now test Tpetra::FEMultiVector
  FEMultiVectorTest femvTest(comm);

  if (me == 0) std::cout << "Testing with int" << std::endl;
  ierr = femvTest.intTest();

  if (me == 0) std::cout << "Testing with vtxLabel" << std::endl;
  ierr = femvTest.vtxLabelTest();

  // Sum up across all processors the number of errors that occurred.
  int gerr = 0;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  // Print PASS/FAIL
  if (me == 0) {
    if (gerr == 0) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL:  " << gerr << " failures" << std::endl;
  }
  
  return 0;
}
