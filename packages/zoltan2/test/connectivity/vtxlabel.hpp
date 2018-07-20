#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"

#include <string>
#include <sstream>
#include <iostream>

#ifndef ZOLTAN2_VTXLABEL_
#define ZOLTAN2_VTXLABEL_
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
       public Teuchos::DirectSerializationTraits<Ordinal, vtxLabel>
{};
}

#endif
