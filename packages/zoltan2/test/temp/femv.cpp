#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"

#include <string>
#include <sstream>
#include <iostream>

struct vtxLabel {
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

///////////////////
// ArithTraits -- arithmetic traits needed for struct vtxLabel
// Needed so that Tpetra compiles.
// Modified from kokkos-kernels/src/Kokkos_ArithTraits.hpp's 
// <int> specialization

namespace Kokkos {
  namespace Details {

    template<>
    class ArithTraits<vtxLabel> {
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
      static KOKKOS_FORCEINLINE_FUNCTION val_type zero() { return 0; }
      static KOKKOS_FORCEINLINE_FUNCTION val_type one() { return 1; }
      static KOKKOS_FORCEINLINE_FUNCTION val_type min() { return INT_MIN; }
      static KOKKOS_FORCEINLINE_FUNCTION val_type max() { return INT_MAX; }
#ifdef NEEDED
      static KOKKOS_FORCEINLINE_FUNCTION mag_type real(const val_type &x) {
        return x.label;
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type imag(const val_type &) {
        return 0;
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type conj(const val_type &x) {
        return x;
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type
      pow(const val_type &x, const val_type &y) {
        return intPowSigned<val_type>(x, y);
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type sqrt(const val_type &x) {
        return static_cast<val_type>( ::sqrt(static_cast<double>(abs(x))));
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type cbrt(const val_type &x) {
        return static_cast<val_type>( ::cbrt(static_cast<double>(abs(x))));
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type exp(const val_type &x) {
        return static_cast<val_type>( ::exp(static_cast<double>(abs(x))));
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type log(const val_type &x) {
        return static_cast<val_type>( ::log(static_cast<double>(abs(x))));
      }
      static KOKKOS_FORCEINLINE_FUNCTION val_type log10(const val_type &x) {
        return static_cast<val_type>( ::log10(static_cast<double>(abs(x))));
      }
      static KOKKOS_FORCEINLINE_FUNCTION mag_type epsilon() { return zero(); }
#endif
      static KOKKOS_FORCEINLINE_FUNCTION val_type nan() { return -1; }
    
      // Backwards compatibility with Teuchos::ScalarTraits.
      typedef mag_type magnitudeType;
#ifdef NEEDED
      typedef val_type halfPrecision;
      typedef val_type doublePrecision;
#endif
    
      static const bool isComplex = false;
      static const bool isOrdinal = true;
      static const bool isComparable = true;
      static const bool hasMachineParameters = false;
      static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude(const val_type &x) {
        return abs(x);
      }
#ifdef NEEDED
      static KOKKOS_FORCEINLINE_FUNCTION val_type conjugate(const val_type &x) {
        return conj(x);
      }
#endif
      static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf(const val_type &) {
        return false;
      }
      static std::string name() { return "vtxLabel"; }
#ifdef NEEDED
      static KOKKOS_FORCEINLINE_FUNCTION val_type squareroot(const val_type &x){
        return sqrt(x);
      }
#endif
    };
  }
}

///////////////////
// Serialization traits needed to copy into MPI buffers
// Because sizeof(vtxLabel) is OK, we'll just use a provided serialization.
template<typename Ordinal>
struct Teuchos::SerializationTraits<Ordinal, vtxLabel> :
       public Teuchos::DirectSerializationTraits<Ordinal, vtxLabel>
{};

///////////////////
class FEMultiVectorTest {
public:
  typedef Tpetra::Map<> map_t;
  typedef map_t::local_ordinal_type lno_t;
  typedef map_t::global_ordinal_type gno_t;

  FEMultiVectorTest(Teuchos::RCP<const Teuchos::Comm<int> > &comm_) :
    me(comm_->getRank()), np(comm_->getSize()),
    nLocalOwned(10), nLocalCopy( np > 1 ? 5 : 0), 
    nVec(2), comm(comm_)
  {
    // Create a map with duplicated entries (mapWithCopies)
    // Each rank has 15 IDs, the last five of which overlap with the next rank.

    const Tpetra::global_size_t nGlobal = np * nLocalOwned;
    lno_t offset = me * nLocalOwned;

    Teuchos::Array<gno_t> gids(nLocalOwned+nLocalCopy);
    for (lno_t i = 0 ; i < nLocalOwned+nLocalCopy; i++)
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

  template <typename femv_t>
  Teuchos::RCP<femv_t> getFEMV()
  {
    // Create FEMultiVector
    typedef Tpetra::Import<lno_t, gno_t> import_t;
    Teuchos::RCP<import_t> importer = rcp(new import_t(mapOwned,
                                                       mapWithCopies));
    Teuchos::RCP<femv_t> femv = rcp(new femv_t(mapOwned, importer,
                                               nVec, true));
    std::cout << me << " FEMV " << femv->getLocalLength() << " "
                                << femv->getGlobalLength() << std::endl;
  
    femv->beginFill();
    for (lno_t i = 0; i < nLocalOwned+nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      femv->replaceGlobalValue(gid, 0, gid);
      femv->replaceGlobalValue(gid, 1, me);
    }
    femv->endFill();

    // Print
    // femv->describe(*Teuchos::getFancyOStream(Teuchos::rcpFromRef(std::cout)),
    //                Teuchos::VERB_HIGH);

    for (int v = 0; v < nVec; v++) {
      std::cout << me << " FEMV[" << v << "] Unique: ";
      auto value = femv->getData(v);
      for (lno_t i = 0; i < nLocalOwned; i++) std::cout << value[i] << " ";
      std::cout << std::endl;
    }

    return femv;
  }

  // Test using scalar_t=int to exercise push capability of FEMultiVector
  int intTest();
  // Test using scalar_t=custom data type 
  int vtxLabelTest();

private:
  int me;
  int np;
  int nLocalOwned;
  int nLocalCopy;
  int nVec;

  Teuchos::RCP<const Teuchos::Comm<int> > comm;

  Teuchos::RCP<const map_t> mapWithCopies;
  Teuchos::RCP<const map_t> mapOwned;
};

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
  }

  {
    auto value = femv->getData(1);
    int tmp = me + (np + me - 1) % np;
    for (lno_t i = 0; i < nLocalCopy; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (value[i] != tmp) 
        std::cout << me << " Error in vec 1 overlap:  gid=" << gid
                        << " value= " << value[i] << " should be " << tmp
                        << std::endl;
        ierr++;
    }
    for (lno_t i = nLocalCopy; i < nLocalOwned; i++) {
      gno_t gid = mapWithCopies->getGlobalElement(i);
      if (value[i] != me) 
        std::cout << me << " Error in vec 1:  gid=" << gid
                        << " value= " << value[i] << " should be " << me
                        << std::endl;
        ierr++;
    }
  }
  return ierr;
}

int FEMultiVectorTest::vtxLabelTest()
{
  typedef vtxLabel scalar_t;
  typedef Tpetra::FEMultiVector<scalar_t, lno_t, gno_t> femv_t;
  int ierr = 0;

  Teuchos::RCP<femv_t> femv = getFEMV<femv_t>();

  // Check results:  after ADD in endFill, 
  // -  overlapping entries of vec 0 should be MAX(local gid, received gid)
  //    nonoverlapping entries of vec 0 should be gid
  // -  overlapping entries of vec 1 should be MAX(me, (np + me-1) % np);
  //    nonoverlapping entries of vec 1 should be me
  
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
  }

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

  if (me == 0) 
    ierr += vtxLabelUnitTest();
  
  FEMultiVectorTest femvTest(comm);

  if (me == 0) std::cout << "Testing with int" << std::endl;
  ierr = femvTest.intTest();

  if (me == 0) std::cout << "Testing with vtxLabel" << std::endl;
  ierr = femvTest.vtxLabelTest();

  int gerr = 0;
  Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);

  if (me == 0) {
    if (gerr == 0) std::cout << "PASS" << std::endl;
    else std::cout << "FAIL:  " << gerr << " failures" << std::endl;
  }
  
  return 0;
}
