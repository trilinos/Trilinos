// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//          PyTrilinos: Python Interfaces to Trilinos Packages
//                 Copyright (2014) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia
// Corporation, the U.S. Government retains certain rights in this
// software.
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
// Questions? Contact William F. Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%{
#include "Teuchos_ScalarTraitsDecl.hpp"
#include "Teuchos_SerializationTraits.hpp"

// Specializations of Teuchos::SerializationTraits<>:
//     unsigned char, short, unsigned short, unsigned int, long,
//     long long, unsigned long long
namespace Teuchos
{
#ifdef TEUCHOS_SERIALIZATIONTRAITS_UNSIGNED_CHAR
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned char>
  : public DirectSerializationTraits<Ordinal,unsigned char>
{};
#endif
#ifdef TEUCHOS_SERIALIZATIONTRAITS_UNSIGNED_SHORT
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned short>
  : public DirectSerializationTraits<Ordinal,unsigned short>
{};
#endif
#ifdef TEUCHOS_SERIALIZATIONTRAITS_UNSIGNED_INT
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned int>
  : public DirectSerializationTraits<Ordinal,unsigned int>
{};
#endif
#ifdef TEUCHOS_SERIALIZATIONTRAITS_UNSIGNED_LONG
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned long>
  : public DirectSerializationTraits<Ordinal,unsigned long>
{};
#endif
#if defined(TEUCHOS_SERIALIZATIONTRAITS_LONG_LONG) || !defined(HAVE_TEUCHOS_LONG_LONG_INT)
template<typename Ordinal>
class SerializationTraits<Ordinal,long long>
  : public DirectSerializationTraits<Ordinal,long long>
{};
#endif
#if defined(TEUCHOS_SERIALIZATIONTRAITS_UNSIGNED_LONG_LONG) || !defined(HAVE_TEUCHOS_LONG_LONG_INT)
template<typename Ordinal>
class SerializationTraits<Ordinal,unsigned long long>
  : public DirectSerializationTraits<Ordinal,unsigned long long>
{};
#endif
}
%}

///////////////////////////////////
// Teuchos::ScalarTraits support //
///////////////////////////////////
%include "Teuchos_ScalarTraitsDecl.hpp"
// Specialization of ScalarTraits<> for type long
%inline 
%{
namespace Teuchos
{
// Type unsigned char
template<>
struct ScalarTraits<unsigned char>
{
  typedef unsigned char magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned char a)
  { return static_cast<unsigned char>(std::fabs(static_cast<double>(a))); }
  static inline unsigned char zero()  { return 0; }
  static inline unsigned char one()   { return 1; }
  static inline unsigned char conjugate(unsigned char x) { return x; }
  static inline unsigned char real(unsigned char x) { return x; }
  static inline unsigned char imag(unsigned char x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned char random() { return std::rand(); }
  static inline std::string name() { return "unsigned char"; }
  static inline unsigned char squareroot(unsigned char x)
  { return (unsigned char) std::sqrt((double) x); }
  static inline unsigned char pow(unsigned char x, unsigned char y)
  { return (unsigned char) std::pow((double)x,(double)y); }
};

// Type unsigned short (now defined in Teuchos_ScalarTraits.hpp)
#if 0
template<>
struct ScalarTraits<unsigned short>
{
  typedef unsigned short magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned short a)
  { return static_cast<unsigned short>(std::fabs(static_cast<double>(a))); }
  static inline unsigned short zero()  { return 0; }
  static inline unsigned short one()   { return 1; }
  static inline unsigned short conjugate(unsigned short x) { return x; }
  static inline unsigned short real(unsigned short x) { return x; }
  static inline unsigned short imag(unsigned short x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned short random() { return std::rand(); }
  static inline std::string name() { return "unsigned short"; }
  static inline unsigned short squareroot(unsigned short x)
  { return (unsigned short) std::sqrt((double) x); }
  static inline unsigned short pow(unsigned short x, unsigned short y)
  { return (unsigned short) std::pow((double)x,(double)y); }
};
#endif

// Type unsigned int
#ifdef TEUCHOS_SCALARTRAITS_UNSIGNED_INT
template<>
struct ScalarTraits<unsigned int>
{
  typedef unsigned int magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned int a)
  { return static_cast<unsigned int>(std::fabs(static_cast<double>(a))); }
  static inline unsigned int zero()  { return 0; }
  static inline unsigned int one()   { return 1; }
  static inline unsigned int conjugate(unsigned int x) { return x; }
  static inline unsigned int real(unsigned int x) { return x; }
  static inline unsigned int imag(unsigned int x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned int random() { return std::rand(); }
  static inline std::string name() { return "unsigned int"; }
  static inline unsigned int squareroot(unsigned int x)
  { return (unsigned int) std::sqrt((double) x); }
  static inline unsigned int pow(unsigned int x, unsigned int y)
  { return (unsigned int) std::pow((double)x,(double)y); }
};
#endif

// Type unsigned long, already taken care of in Teuchos_ScalarTraits.hpp
#ifdef TEUCHOS_SCALARTRAITS_UNSIGNED_LONG
template<>
struct ScalarTraits<unsigned long>
{
  typedef unsigned long magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned long a)
  { return static_cast<unsigned long>(std::fabs(static_cast<double>(a))); }
  static inline unsigned long zero()  { return 0; }
  static inline unsigned long one()   { return 1; }
  static inline unsigned long conjugate(unsigned long x) { return x; }
  static inline unsigned long real(unsigned long x) { return x; }
  static inline unsigned long imag(unsigned long x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned long random() { return std::rand(); }
  static inline std::string name() { return "unsigned long"; }
  static inline unsigned long squareroot(unsigned long x)
  { return (unsigned long) std::sqrt((double) x); }
  static inline unsigned long pow(unsigned long x, unsigned long y)
  { return (unsigned long) std::pow((double)x,(double)y); }
};
#endif

#if defined(TEUCHOS_SCALARTRAITS_LONG_LONG) || !defined(HAVE_TEUCHOS_LONG_LONG_INT)
// Type long long
template<>
struct ScalarTraits<long long>
{
  typedef long long magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(long long a)
  { return static_cast<long long>(std::fabs(static_cast<double>(a))); }
  static inline long long zero()  { return 0; }
  static inline long long one()   { return 1; }
  static inline long long conjugate(long long x) { return x; }
  static inline long long real(long long x) { return x; }
  static inline long long imag(long long x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline long long random() { return (-1 + 2*rand()); }
  static inline std::string name() { return "long long"; }
  static inline long long squareroot(long long x)
  { return (long long) std::sqrt((double) x); }
  static inline long long pow(long long x, long long y)
  { return (long long) std::pow((double)x,(double)y); }
};
#endif

#if defined(TEUCHOS_SCALARTRAITS_UNSIGNED_LONG_LONG) || !defined(HAVE_TEUCHOS_LONG_LONG_INT)
// Type unsigned long long
template<>
struct ScalarTraits<unsigned long long>
{
  typedef unsigned long long magnitudeType;
  static const bool isComplex = false;
  static const bool isComparable = true;
  static const bool hasMachineParameters = false;
  static inline magnitudeType magnitude(unsigned long long a)
  { return static_cast<unsigned long long>(std::fabs(static_cast<double>(a))); }
  static inline unsigned long long zero()  { return 0; }
  static inline unsigned long long one()   { return 1; }
  static inline unsigned long long conjugate(unsigned long long x) { return x; }
  static inline unsigned long long real(unsigned long long x) { return x; }
  static inline unsigned long long imag(unsigned long long x) { return 0; }
  static inline void seedrandom(unsigned int s) { std::srand(s); 
#ifdef __APPLE__
    random(); // throw away first random number to address bug 3655
#endif
  }
  static inline unsigned long long random() { return std::rand(); }
  static inline std::string name() { return "unsigned long long"; }
  static inline unsigned long long squareroot(unsigned long long x)
  { return (unsigned long long) std::sqrt((double) x); }
  static inline unsigned long long pow(unsigned long long x, unsigned long long y)
  { return (unsigned long long) std::pow((double)x,(double)y); }
};
#endif

}
%}
%ignore Teuchos::ScalarTraits< char >::eps  ;
%ignore Teuchos::ScalarTraits< char >::sfmin;
%ignore Teuchos::ScalarTraits< char >::base ;
%ignore Teuchos::ScalarTraits< char >::prec ;
%ignore Teuchos::ScalarTraits< char >::t    ;
%ignore Teuchos::ScalarTraits< char >::rnd  ;
%ignore Teuchos::ScalarTraits< char >::emin ;
%ignore Teuchos::ScalarTraits< char >::rmin ;
%ignore Teuchos::ScalarTraits< char >::emax ;
%ignore Teuchos::ScalarTraits< char >::rmax ;
%ignore Teuchos::ScalarTraits< char >::nan  ;
%ignore Teuchos::ScalarTraits< char >::pi   ;
%template(ScalarTraitsByte) Teuchos::ScalarTraits< char >;

%template(ScalarTraitsUbyte) Teuchos::ScalarTraits< unsigned char >;

%ignore Teuchos::ScalarTraits< short >::eps     ;
%ignore Teuchos::ScalarTraits< short >::sfmin   ;
%ignore Teuchos::ScalarTraits< short >::base    ;
%ignore Teuchos::ScalarTraits< short >::prec    ;
%ignore Teuchos::ScalarTraits< short >::t       ;
%ignore Teuchos::ScalarTraits< short >::rnd     ;
%ignore Teuchos::ScalarTraits< short >::emin    ;
%ignore Teuchos::ScalarTraits< short >::rmin    ;
%ignore Teuchos::ScalarTraits< short >::emax    ;
%ignore Teuchos::ScalarTraits< short >::rmax    ;
%ignore Teuchos::ScalarTraits< short >::nan     ;
%ignore Teuchos::ScalarTraits< short >::isnaninf;
%ignore Teuchos::ScalarTraits< short >::pi      ;
%template(ScalarTraitsShort) Teuchos::ScalarTraits< short >;

%ignore Teuchos::ScalarTraits< unsigned short >::eps     ;
%ignore Teuchos::ScalarTraits< unsigned short >::sfmin   ;
%ignore Teuchos::ScalarTraits< unsigned short >::base    ;
%ignore Teuchos::ScalarTraits< unsigned short >::prec    ;
%ignore Teuchos::ScalarTraits< unsigned short >::t       ;
%ignore Teuchos::ScalarTraits< unsigned short >::rnd     ;
%ignore Teuchos::ScalarTraits< unsigned short >::emin    ;
%ignore Teuchos::ScalarTraits< unsigned short >::rmin    ;
%ignore Teuchos::ScalarTraits< unsigned short >::emax    ;
%ignore Teuchos::ScalarTraits< unsigned short >::rmax    ;
%ignore Teuchos::ScalarTraits< unsigned short >::nan     ;
%ignore Teuchos::ScalarTraits< unsigned short >::isnaninf;
%ignore Teuchos::ScalarTraits< unsigned short >::pi      ;
%template(ScalarTraitsUshort) Teuchos::ScalarTraits< unsigned short >;

%ignore Teuchos::ScalarTraits< int >::eps  ;
%ignore Teuchos::ScalarTraits< int >::sfmin;
%ignore Teuchos::ScalarTraits< int >::base ;
%ignore Teuchos::ScalarTraits< int >::prec ;
%ignore Teuchos::ScalarTraits< int >::t    ;
%ignore Teuchos::ScalarTraits< int >::rnd  ;
%ignore Teuchos::ScalarTraits< int >::emin ;
%ignore Teuchos::ScalarTraits< int >::rmin ;
%ignore Teuchos::ScalarTraits< int >::emax ;
%ignore Teuchos::ScalarTraits< int >::rmax ;
%ignore Teuchos::ScalarTraits< int >::nan  ;
%ignore Teuchos::ScalarTraits< int >::pi   ;
%template(ScalarTraitsInt) Teuchos::ScalarTraits< int >;

%ignore Teuchos::ScalarTraits< unsigned int >::eps     ;
%ignore Teuchos::ScalarTraits< unsigned int >::sfmin   ;
%ignore Teuchos::ScalarTraits< unsigned int >::base    ;
%ignore Teuchos::ScalarTraits< unsigned int >::prec    ;
%ignore Teuchos::ScalarTraits< unsigned int >::t       ;
%ignore Teuchos::ScalarTraits< unsigned int >::rnd     ;
%ignore Teuchos::ScalarTraits< unsigned int >::emin    ;
%ignore Teuchos::ScalarTraits< unsigned int >::rmin    ;
%ignore Teuchos::ScalarTraits< unsigned int >::emax    ;
%ignore Teuchos::ScalarTraits< unsigned int >::rmax    ;
%ignore Teuchos::ScalarTraits< unsigned int >::nan     ;
%ignore Teuchos::ScalarTraits< unsigned int >::isnaninf;
%ignore Teuchos::ScalarTraits< unsigned int >::pi      ;
%template(ScalarTraitsUint) Teuchos::ScalarTraits< unsigned int >;

%ignore Teuchos::ScalarTraits< long >::eps     ;
%ignore Teuchos::ScalarTraits< long >::sfmin   ;
%ignore Teuchos::ScalarTraits< long >::base    ;
%ignore Teuchos::ScalarTraits< long >::prec    ;
%ignore Teuchos::ScalarTraits< long >::t       ;
%ignore Teuchos::ScalarTraits< long >::rnd     ;
%ignore Teuchos::ScalarTraits< long >::emin    ;
%ignore Teuchos::ScalarTraits< long >::rmin    ;
%ignore Teuchos::ScalarTraits< long >::emax    ;
%ignore Teuchos::ScalarTraits< long >::rmax    ;
%ignore Teuchos::ScalarTraits< long >::nan     ;
%ignore Teuchos::ScalarTraits< long >::isnaninf;
%ignore Teuchos::ScalarTraits< long >::pi      ;
%template(ScalarTraitsLong) Teuchos::ScalarTraits< long >;

%ignore Teuchos::ScalarTraits< unsigned long >::eps     ;
%ignore Teuchos::ScalarTraits< unsigned long >::sfmin   ;
%ignore Teuchos::ScalarTraits< unsigned long >::base    ;
%ignore Teuchos::ScalarTraits< unsigned long >::prec    ;
%ignore Teuchos::ScalarTraits< unsigned long >::t       ;
%ignore Teuchos::ScalarTraits< unsigned long >::rnd     ;
%ignore Teuchos::ScalarTraits< unsigned long >::emin    ;
%ignore Teuchos::ScalarTraits< unsigned long >::rmin    ;
%ignore Teuchos::ScalarTraits< unsigned long >::emax    ;
%ignore Teuchos::ScalarTraits< unsigned long >::rmax    ;
%ignore Teuchos::ScalarTraits< unsigned long >::nan     ;
%ignore Teuchos::ScalarTraits< unsigned long >::isnaninf;
%ignore Teuchos::ScalarTraits< unsigned long >::pi      ;
%template(ScalarTraitsUlong) Teuchos::ScalarTraits< unsigned long >;

%ignore Teuchos::ScalarTraits< long long >::eps     ;
%ignore Teuchos::ScalarTraits< long long >::sfmin   ;
%ignore Teuchos::ScalarTraits< long long >::base    ;
%ignore Teuchos::ScalarTraits< long long >::prec    ;
%ignore Teuchos::ScalarTraits< long long >::t       ;
%ignore Teuchos::ScalarTraits< long long >::rnd     ;
%ignore Teuchos::ScalarTraits< long long >::emin    ;
%ignore Teuchos::ScalarTraits< long long >::rmin    ;
%ignore Teuchos::ScalarTraits< long long >::emax    ;
%ignore Teuchos::ScalarTraits< long long >::rmax    ;
%ignore Teuchos::ScalarTraits< long long >::nan     ;
%ignore Teuchos::ScalarTraits< long long >::isnaninf;
%ignore Teuchos::ScalarTraits< long long >::pi      ;

%template(ScalarTraitsLongLong) Teuchos::ScalarTraits< long long >;

%ignore Teuchos::ScalarTraits< unsigned long long >::eps     ;
%ignore Teuchos::ScalarTraits< unsigned long long >::sfmin   ;
%ignore Teuchos::ScalarTraits< unsigned long long >::base    ;
%ignore Teuchos::ScalarTraits< unsigned long long >::prec    ;
%ignore Teuchos::ScalarTraits< unsigned long long >::t       ;
%ignore Teuchos::ScalarTraits< unsigned long long >::rnd     ;
%ignore Teuchos::ScalarTraits< unsigned long long >::emin    ;
%ignore Teuchos::ScalarTraits< unsigned long long >::rmin    ;
%ignore Teuchos::ScalarTraits< unsigned long long >::emax    ;
%ignore Teuchos::ScalarTraits< unsigned long long >::rmax    ;
%ignore Teuchos::ScalarTraits< unsigned long long >::nan     ;
%ignore Teuchos::ScalarTraits< unsigned long long >::isnaninf;
%ignore Teuchos::ScalarTraits< unsigned long long >::pi      ;

%template(ScalarTraitsUlongLong) Teuchos::ScalarTraits< unsigned long long >;

%template(ScalarTraitsFloat) Teuchos::ScalarTraits< float >;

%template(ScalarTraitsDouble) Teuchos::ScalarTraits< double >;

%pythoncode
%{
def ScalarTraits(scalarType):
    """
    ScalarTraits(str scalarType) -> ScalarTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which traits are requested.  Limited NumPy-style type
    specification is supported:

        'b'  byte
        'B'  unsigned byte
        'h'  short
        'H'  unsigned short
        'i'  int
        'I'  unsigned int
        'l'  long
        'L'  unsigned long
        'q'  long long
        'Q'  unsigned long long
        'f'  float
        'd'  double
    """
    if   scalarType == 'b': return ScalarTraitsByte()
    elif scalarType == 'B': return ScalarTraitsUbyte()
    elif scalarType == 'h': return ScalarTraitsShort()
    elif scalarType == 'H': return ScalarTraitsUshort()
    elif scalarType == 'i': return ScalarTraitsInt()
    elif scalarType == 'I': return ScalarTraitsUint()
    elif scalarType == 'l': return ScalarTraitsLong()
    elif scalarType == 'L': return ScalarTraitsUlong()
    elif scalarType == 'q': return ScalarTraitsLongLong()
    elif scalarType == 'Q': return ScalarTraitsUlongLong()
    elif scalarType == 'f': return ScalarTraitsFloat()
    elif scalarType == 'd': return ScalarTraitsDouble()
    raise NotImplementedError, "ScalarTraits for " + repr(scalarType) + " not supported"
%}

//////////////////////////////////////////
// Teuchos::SerializationTraits support //
//////////////////////////////////////////
%define
%teuchos_serialtraits_for_type(type,Type)
%ignore Teuchos::DirectSerializationTraits<long,type>::convertToCharPtr(const type*);
%ignore Teuchos::DirectSerializationTraits<long,type>::convertToCharPtr(const char*);
%ignore Teuchos::DirectSerializationTraits<long,type>::convertFromCharPtr(const char*);
%template(DirectSerializationTraits ## Type)
    Teuchos::DirectSerializationTraits<long, type>;
%ignore Teuchos::SerializationTraits<long,type>::convertToCharPtr(const type*);
%ignore Teuchos::SerializationTraits<long,type>::convertToCharPtr(const char*);
%ignore Teuchos::SerializationTraits<long,type>::convertFromCharPtr(const char*);
%template(SerializationTraits ## Type) Teuchos::SerializationTraits<long, type>;
%enddef
%include "Teuchos_SerializationTraits.hpp"
%teuchos_serialtraits_for_type(char              , Byte     )
%teuchos_serialtraits_for_type(unsigned char     , Ubyte    )
%teuchos_serialtraits_for_type(short             , Short    )
%teuchos_serialtraits_for_type(unsigned short    , Ushort   )
%teuchos_serialtraits_for_type(int               , Int      )
%teuchos_serialtraits_for_type(unsigned int      , Uint     )
%teuchos_serialtraits_for_type(long              , Long     )
%teuchos_serialtraits_for_type(unsigned long     , Ulong    )
%teuchos_serialtraits_for_type(long long         , LongLong )
%teuchos_serialtraits_for_type(unsigned long long, UlongLong)
%teuchos_serialtraits_for_type(float             , Float    )
%teuchos_serialtraits_for_type(double            , Double   )
%pythoncode
%{
def DirectSerializationTraits(scalarType):
    """
    DirectSerializationTraits(str scalarType) -> DirectSerializationTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which direct serialization is requested.  Limited NumPy-style type
    specification is supported:

        'b'  byte
        'B'  unsigned byte
        'h'  short
        'H'  unsigned short
        'i'  int
        'I'  unsigned int
        'l'  long
        'L'  unsigned long
        'q'  long long
        'Q'  unsigned long long
        'f'  float
        'd'  double
    """
    if   scalarType == 'b': return DirectSerializationTraitsByte()
    elif scalarType == 'B': return DirectSerializationTraitsUbyte()
    elif scalarType == 'h': return DirectSerializationTraitsShort()
    elif scalarType == 'H': return DirectSerializationTraitsUshort()
    elif scalarType == 'i': return DirectSerializationTraitsInt()
    elif scalarType == 'I': return DirectSerializationTraitsUint()
    elif scalarType == 'l': return DirectSerializationTraitsLong()
    elif scalarType == 'L': return DirectSerializationTraitsUlong()
    elif scalarType == 'q': return DirectSerializationTraitsLongLong()
    elif scalarType == 'Q': return DirectSerializationTraitsUlongLong()
    elif scalarType == 'f': return DirectSerializationTraitsFloat()
    elif scalarType == 'd': return DirectSerializationTraitsDouble()
    raise NotImplementedError, "DirectSerializationTraits for " + repr(scalarType) + \
	  " not supported"

def SerializationTraits(scalarType):
    """
    SerializationTraits(str scalarType) -> SerializationTraits<...>

    The scalarType argument is for specifying the type of scalar for
    which serialization is requested.  Limited NumPy-style type
    specification is supported:

        'b'  byte
        'B'  unsigned byte
        'h'  short
        'H'  unsigned short
        'i'  int
        'I'  unsigned int
        'l'  long
        'L'  unsigned long
        'q'  long long
        'Q'  unsigned long long
        'f'  float
        'd'  double
    """
    if   scalarType == 'b': return SerializationTraitsByte()
    elif scalarType == 'B': return SerializationTraitsUbyte()
    elif scalarType == 'h': return SerializationTraitsShort()
    elif scalarType == 'H': return SerializationTraitsUshort()
    elif scalarType == 'i': return SerializationTraitsInt()
    elif scalarType == 'I': return SerializationTraitsUint()
    elif scalarType == 'l': return SerializationTraitsLong()
    elif scalarType == 'L': return SerializationTraitsUlong()
    elif scalarType == 'q': return SerializationTraitsLongLong()
    elif scalarType == 'Q': return SerializationTraitsUlongLong()
    elif scalarType == 'f': return SerializationTraitsFloat()
    elif scalarType == 'd': return SerializationTraitsDouble()
    raise NotImplementedError, "SerializationTraits for " + repr(scalarType) + \
	  " not supported"
%}
