#include "Tpetra_Core.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_FEMultiVector.hpp"
#include "Tpetra_MultiVector_decl.hpp"
#include "Tpetra_MultiVector_def.hpp"

#include <string>
#include <sstream>
#include <iostream>
#include <queue>
#include <vector>

#ifndef ZOLTAN2_VTXLABEL_
#define ZOLTAN2_VTXLABEL_

namespace Zoltan2{
        
  //Graph datastructure that represents the local
  //graph in ice sheet propagation problems
  template<typename lno_t,typename offset_t>
  class icePropGraph {
  public:
    //n represents the number of vertices in the local graph
    size_t nVtx;
    //m represents the number of edges in the local graph
    size_t nEdges;
    //out_array is the adjacency array for each vertex
    Teuchos::Array<lno_t> out_array;
    //out_degree_list is the offset array for indexing the out_array
    Teuchos::ArrayView<const offset_t> out_degree_list;

    //member function to return the degree of a given vertex.
    lno_t out_degree(lno_t vert){
      return (out_degree_list[vert+1] - out_degree_list[vert]);
    }
    //member function to return the neighbors for a given vertex
    lno_t* out_vertices(lno_t vert){
      lno_t* raw_out_array = out_array.getRawPtr();
      return (&raw_out_array[out_degree_list[vert]]);
    }
  }; 
        

  //Enum that denotes
  enum IcePropGrounding_Status {ICEPROPGS_FULL=2, ICEPROPGS_HALF=1, 
                                ICEPROPGS_NONE = 0};
  // Struct representing a vertex label.
  // We define our own "addition" for these labels.
  // Later, we'll create a Tpetra::FEMultiVector of these labels.
  template<typename lno_t, typename gno_t>
  class IcePropVtxLabel {
  public: 
    //The local ID for the vertex represented by this label
    lno_t id;
    
    //the global ID of a grounded vertex
    gno_t first_label;
    
    //the global ID of the vertex that sent the label
    gno_t first_sender;
    
    //this field indicates whether or not the first label is currently set
    //(this is necessary for unsigned global types)
    bool  first_used;
    
    //the global ID of a second grounded vertex
    gno_t second_label;

    //the global ID of the vertex that sent the second identifier
    gno_t second_sender;

    //this field indicates whether of not the second label is currently set
    //(this is necessary for unsigned global types)
    bool  second_used;

    //this field will be used to find biconnected components
    int bcc_name;

    //this flag indicates whether or not this vertex is a potential articulation point
    bool is_art;

    // Constructors
    IcePropVtxLabel(int idx_, int first_ = -1, int first_sender_ = -1, 
                              bool first_used_ = false, 
                              int second_ = -1, int second_sender_ = -1, 
                              bool second_used_ = false,
                              bool art_ = false, int bcc_name_ = -1) { 
      id = idx_;
      first_label = first_;
      first_sender = first_sender_;
      first_used = first_used_;
      second_label = second_;
      second_sender = second_sender_; 
      second_used = second_used_;
      is_art = art_;
      bcc_name = bcc_name_;
    }
    IcePropVtxLabel() {
      id = -1;
      first_label = -1;
      first_sender = -1;
      first_used = false;
      second_label = -1;
      second_sender = -1;
      second_used = false;
      is_art = false; 
      bcc_name = -1;
    }

    IcePropVtxLabel(volatile const IcePropVtxLabel& other){
      id = other.id;
      first_label = other.first_label;
      first_sender = other.first_sender;
      first_used = other.first_used;
      second_label = other.second_label;
      second_sender = other.second_sender;
      second_used = other.second_used;
      is_art = other.is_art;
      bcc_name = other.bcc_name;
    }
    IcePropVtxLabel(const IcePropVtxLabel& other){
      id = other.id;
      first_label = other.first_label;
      first_sender = other.first_sender;
      first_used = other.first_used;
      second_label = other.second_label;
      second_sender = other.second_sender;
      second_used = other.second_used;
      is_art = other.is_art;
      bcc_name = other.bcc_name;
    }
    // IcePropVtxLabel assignment
    volatile IcePropVtxLabel operator=(const IcePropVtxLabel& other) volatile{
      id = other.id;
      first_label = other.first_label;
      first_sender = other.first_sender;
      first_used = other.first_used;
      second_label = other.second_label;
      second_sender = other.second_sender;
      second_used = other.second_used;
      is_art = other.is_art;
      bcc_name = other.bcc_name;
      return *this; 
    }

    IcePropVtxLabel& operator=(const IcePropVtxLabel& other) {
      id = other.id;
      first_label = other.first_label;
      first_sender = other.first_sender;
      first_used = other.first_used;
      second_label = other.second_label;
      second_sender = other.second_sender;
      second_used = other.second_used;
      is_art = other.is_art;
      bcc_name = other.bcc_name;
      return *this; 
    } 

    // int assignment
    IcePropVtxLabel& operator=(const int& other) { 
      first_label = other;
      first_sender = other;
      first_used = true;
      second_label = -1;
      second_sender = -1;
      second_used = false;
      return *this;
    }

    //stub overloads for FEMultiVector ETI
    
    friend IcePropVtxLabel operator- (const IcePropVtxLabel& neg){
      IcePropVtxLabel ret;
      return ret;
    }

    friend IcePropVtxLabel operator*(const IcePropVtxLabel& one, const IcePropVtxLabel& other){
      IcePropVtxLabel ret;
      return ret;
    }

    friend IcePropVtxLabel operator+(const IcePropVtxLabel& one, const IcePropVtxLabel& other){
      IcePropVtxLabel ret;
      return ret;
    }

    friend IcePropVtxLabel operator-(const IcePropVtxLabel& one, const IcePropVtxLabel& other){
      IcePropVtxLabel ret;
      return ret;
    }

    friend IcePropVtxLabel operator/(const IcePropVtxLabel& one, const IcePropVtxLabel& other){
      IcePropVtxLabel ret;
      return ret;
    }

    friend bool operator!=(const IcePropVtxLabel& lhs, const IcePropVtxLabel& rhs){
      return !(lhs == rhs);
    }
    
    // += overload
    // for communicating copy's labels over processor boundaries.
    IcePropVtxLabel& operator+=(const IcePropVtxLabel& copy) {
      IcePropGrounding_Status owned_gs = getGroundingStatus();
      IcePropGrounding_Status copy_gs = copy.getGroundingStatus();
      //The only cases we care about are 
      //owned	copy	
      //NONE  <   HALF
      //NONE  <	FULL
      //HALF  ==	HALF
      //HALF  <	FULL
    
      //handles NONE < HALF, HALF < FULL
      if(owned_gs < copy_gs){
        first_label = copy.first_label;
        first_sender = copy.first_sender;
        first_used = copy.first_used;
        second_label = copy.second_label;
        second_sender = copy.second_sender;
        second_used = copy.second_used;
        bcc_name = copy.bcc_name;
      //handles HALF == HALF
      } else if(owned_gs == copy_gs && owned_gs == ICEPROPGS_HALF){
        if(copy.first_label != first_label){
          second_label = copy.first_label;
          second_sender = copy.first_sender;
          second_used = copy.first_used;
        }
      }
      
      return *this;
    }
    
    // IcePropVtxLabel equality overload
    friend bool operator==(const IcePropVtxLabel& lhs, const IcePropVtxLabel& rhs) {
      return ((lhs.first_label == rhs.first_label)&&(lhs.first_sender == rhs.first_sender)&&(lhs.second_label == rhs.second_label)&&(lhs.second_sender == rhs.second_sender));
    }
    // int equality overload
    friend bool operator==(const IcePropVtxLabel& lhs, const int& rhs) {
      return ((lhs.first_label == rhs)&&(lhs.first_sender == rhs));
    }
    // output stream overload
    friend std::ostream& operator<<(std::ostream& os, const IcePropVtxLabel& a) {
      os<<a.id<<": "<< a.first_label<<", "<<a.first_sender<<"; "<<a.second_label<<", "<<a.second_sender;
      return os;
    }
    IcePropGrounding_Status getGroundingStatus() const {
      return (IcePropGrounding_Status)((first_used) + (second_used));
    }
  };

}//end namespace Zoltan2
        
/////////////////////////////////////////////////////////////////////////
// ArithTraits -- arithmetic traits needed for struct IcePropVtxLabel
// Needed so that Tpetra compiles.
// Not all functions were needed; this is a subset of ArithTraits' traits.
// Modified from kokkos-kernels/src/Kokkos_ArithTraits.hpp's 
// <int> specialization
namespace Kokkos {
  namespace Details {

    template<>
    class ArithTraits<Zoltan2::IcePropVtxLabel<
                      Tpetra::Map<>::local_ordinal_type,
                      Tpetra::Map<>::global_ordinal_type> > { 
                      // specialized for IcePropVtxLabel struct
    public:
      typedef Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type,
                                 Tpetra::Map<>::global_ordinal_type> val_type;
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
	return (x.first_label >= 0 ? x.first_label : -(x.first_label));
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
      static KOKKOS_FORCEINLINE_FUNCTION val_type conj( const val_type &x){ return x; }
      static KOKKOS_FORCEINLINE_FUNCTION magnitudeType magnitude(
	const val_type &x) 
      {
	return abs(x);
      }
      static KOKKOS_FORCEINLINE_FUNCTION bool isnaninf(const val_type &) {
	return false;
      }
      static std::string name() { return "Zoltan2::IcePropVtxLabel"; }
    };
  }

}

namespace Kokkos {
  
  template<class Generator>
  struct rand<Generator, Zoltan2::IcePropVtxLabel<
                         Tpetra::Map<>::local_ordinal_type, 
                         Tpetra::Map<>::global_ordinal_type> > {
    typedef Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type, 
                                    Tpetra::Map<>::global_ordinal_type> Scalar;
    //Max value returned by draw(Generator& gen)
    KOKKOS_INLINE_FUNCTION static Scalar max() { return Scalar();}
    //Returns a value between zero and max()
    KOKKOS_INLINE_FUNCTION static Scalar draw(Generator& gen) {
      return Scalar();
    }
    //Returns a value between zero and range()
    //Note: for floating point values range can be larger than max()
    KOKKOS_INLINE_FUNCTION static Scalar draw(Generator& gen, 
                                              const Scalar& range){
      return Scalar();
    }
    //Return value between start and end
    KOKKOS_INLINE_FUNCTION static Scalar draw(Generator& gen, 
                                              const Scalar& start,
                                              const Scalar& end){
      return Scalar();
    }
  };

}
/////////////////////////////////////////////////////////////////////////////
// Teuchos::SerializationTraits are needed to copy vtxLabels into MPI buffers
// Because sizeof(vtxLabel) works for struct vtxLabel, we'll use a 
// provided serialization of vtxLabel into char*.
namespace Teuchos{
template<typename Ordinal>
struct SerializationTraits<Ordinal, Zoltan2::IcePropVtxLabel<
                                      Tpetra::Map<>::local_ordinal_type, 
                                      Tpetra::Map<>::global_ordinal_type> >:
       public Teuchos::DirectSerializationTraits<Ordinal, 
                   Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type, 
                   Tpetra::Map<>::global_ordinal_type> >
{};

template <>
struct ScalarTraits<Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type,
                                        Tpetra::Map<>::global_ordinal_type> > {
  public:
    typedef Zoltan2::IcePropVtxLabel<Tpetra::Map<>::local_ordinal_type, 
                                  Tpetra::Map<>::global_ordinal_type> val_type;
    typedef Tpetra::Map<>::global_ordinal_type magnitudeType;
    
    static const bool isComplex = false;
    static const bool isOrdinal = false;
    static const bool isComparable = false;
    static const bool hasMachineParameters = false;

    static magnitudeType eps()   { return 0;}
    static magnitudeType sfmin() { return 1;}
    static magnitudeType base()  { return 2;}
    static magnitudeType prec()  { return 0;}
    static magnitudeType t()     { return 0;}
    static magnitudeType rnd()   { return 0;}
    static magnitudeType emin()  { return 0;}
    static magnitudeType rmin()  { return 0;}
    static magnitudeType emax()  { return 0;}
    static magnitudeType rmax()  { return 0;}
    static magnitudeType magnitude(val_type a) {return 1;}
    static val_type zero() { return val_type(0);}
    static val_type one()  { return val_type(1);}
    static magnitudeType real(val_type a) {return 1;}
    static magnitudeType imag(val_type a) {return 0;}
    static val_type conjugate(val_type a) {return a;}
    static val_type nan() {return val_type();}
    static bool isnaninf(const val_type &x) {return false;}
    static void seedrandom(unsigned int s) {}
    static val_type random() {return val_type();}
    static std::string name(){return "Zoltan2::IcePropVtxLabel<lno_t, gno_t>";}
    static val_type squareroot(val_type x) {return x;}
    static val_type pow(val_type x, val_type y){ return x;}
    static val_type pi() {return val_type();}
};
}//end namespace Teuchos

#endif
