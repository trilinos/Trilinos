#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <string>
#include <sstream>
#include <iostream>

namespace allGood {

static int counter = 0;

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
  // int assignment
  vtxLabel& operator=(const int& other) { 
    label = other; return *this;
  }
  // += overload
  vtxLabel& operator+=(const vtxLabel& other) { 
    counter++;
    label = -label - other.label; return *this;
  }
  // addition overload
  friend vtxLabel operator+(const vtxLabel& lhs, const vtxLabel& rhs) {
    counter++;
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

}  // end namespace allGood

/////////////////////////////////////////////////////////////////////////
// ArithTraits -- arithmetic traits needed for struct vtxLabel
// Needed so that Tpetra compiles.
// Not all functions were needed; this is a subset of ArithTraits' traits.
// Modified from kokkos-kernels/src/Kokkos_ArithTraits.hpp's 
// <int> specialization

namespace Kokkos {
  namespace Details {

    template<>
    class ArithTraits<allGood::vtxLabel> {  // specialized for vtxLabel struct
    public:
      typedef allGood::vtxLabel val_type;
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
      static std::string name() { return "allGood::vtxLabel"; }
    };
  }
}


/////////////////////////////////////////////////////////////////////////////
// Teuchos::SerializationTraits are needed to copy vtxLabels into MPI buffers
// Because sizeof(vtxLabel) works for struct vtxLabel, we'll use a 
// provided serialization of vtxLabel into char*.
namespace Teuchos {
template<typename Ordinal>
struct SerializationTraits<Ordinal, allGood::vtxLabel> :
       public Teuchos::DirectSerializationTraits<Ordinal, allGood::vtxLabel>
{};
}  // end namespace Teuchos

/////////////////////////////////////////////////////////////////////////////
namespace allGood {

class OverloadTest {
public:

  typedef Tpetra::Map<> map_t;
  typedef map_t::local_ordinal_type lno_t;
  typedef map_t::global_ordinal_type gno_t;
  typedef vtxLabel scalar_t;
  typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;

  // Constructor assigns vertices to processors and builds maps with and 
  // without copies
  OverloadTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm_) :
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

    typedef Tpetra::Import<lno_t, gno_t> import_t;
    Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                       mapWithCopies));
    femv = rcp(new femv_t(mapOwned, importer, nVec, true));
  
    femv->beginFill();
    for (lno_t i = 0; i < nLocalOwned + nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      femv->replaceGlobalValue(gid, 0, gid);
      femv->replaceGlobalValue(gid, 1, me);
    }
    femv->endFill();

    printFEMV(*femv, "AfterFill");

    femv->doSourceToTarget(Tpetra::REPLACE);

    printFEMV(*femv, "AfterReplace");
  }

  void printFEMV(femv_t &femv, const char *msg) 
  {
    for (int v = 0; v < nVec; v++) {
      std::cout << me << " OWNED " << msg << " FEMV[" << v << "] Owned: ";
      auto value = femv.getData(v);
      for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
      std::cout << std::endl;
    }
    femv.switchActiveMultiVector();  // Needed to print copies
    for (int v = 0; v < nVec; v++) {
      std::cout << me << " WITHCOPIES " << msg << " FEMV[" << v << "] Owned: ";
      auto value = femv.getData(v);
      for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
      std::cout << " Copies: ";
      for (lno_t i = nLocalOwned; i < nLocalOwned+nLocalCopy; i++) 
        std::cout << value[i] << " ";
      std::cout << std::endl;
    }
    femv.switchActiveMultiVector();  // Restore state

    std::cout << me << " counter = " << counter << std::endl;
  }

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

  Teuchos::RCP<femv_t> femv;
};

}  // end namespace allGood

/////////////////////////////////////////////////////////////////////

int main(int narg, char **arg)
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  allGood::OverloadTest test(comm);

  std::cout << me << "  From Main:  counter = " 
            << allGood::counter << std::endl;

  // Print PASS/FAIL
  if (me == 0) {
    std::cout << "PASS" << std::endl;
  }
  
  return 0;
}
