/*
// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#ifndef INTREPID2_KOKKOS_RANK_HPP
#define INTREPID2_KOKKOS_RANK_HPP
#include "Kokkos_Core.hpp"
#include "Kokkos_DynRankView.hpp"
#include "Sacado.hpp"

namespace Intrepid2{

template<typename field>
struct is_mdfield{
  enum {value=false};
};


/*template<typename T, typename = void>
struct conditional_eSpace : Kokkos::Serial { };

template<typename T>
struct conditional_eSpace<T, decltype(std::declval<T>().execution_space, void())> : T { };
*/

template<class T>
struct Void {
  typedef void type;
};


typedef Kokkos::DefaultExecutionSpace::size_type index_type;

template<class T, class U = void>
struct conditional_eSpace {
	typedef Kokkos::HostSpace::execution_space execution_space;
};

  // UVM confuses the space allocated in intrepid, so for now redirect to host space only
template<class T>
struct conditional_eSpace<T, typename Intrepid2::Void<typename T::execution_space>::type > {
  //typedef typename T::execution_space execution_space;
  typedef Kokkos::HostSpace::execution_space execution_space;
};

template<class A>
struct CheckType{static const bool value = false; };

template<class A>
struct Rank{static const int value = -1;};

template<class Array, class ... Params>
struct Rank<Kokkos::View<Array,Params...> >{
  static const index_type value=Kokkos::View<Array,Params...>::Rank;
};

template<class Array, class ... Params>
struct Rank<const Kokkos::View<Array,Params...> >{
  static const index_type value=Kokkos::View<Array,Params...>::Rank;
};

template<class Array, class ... Params>
struct CheckType<Kokkos::View<Array,Params...> >{
  static const bool value = true;
};

template<class Array, class ... Params>
struct CheckType<const Kokkos::View<Array,Params...> >{
  static const bool value = true;
};

template<class A,class Scalar>
struct Return_Type{
  typedef Scalar& return_type;
  typedef Scalar const_return_type;
};

template<class Array, class ... Params, class Scalar>
struct Return_Type<const Kokkos::View<Array,Params...>, Scalar>{
  typedef Kokkos::View<Array,Params...> ViewType;
  typedef typename ViewType::reference_type return_type;
  typedef typename ViewType::reference_type const_return_type;
};

template<class Array, class ... Params,class Scalar>
struct Return_Type< Kokkos::View<Array,Params...>, Scalar>{
  typedef Kokkos::View<Array,Params...> ViewType;
  typedef typename ViewType::reference_type return_type;
  typedef typename ViewType::reference_type const_return_type;
};


  // Kokkos::DynRankView specializations
  template<class Array, class ... Params>
  struct Rank<Kokkos::DynRankView<Array,Params...> >{
    static const int value = -1;
  };
  
  template<class Array, class ... Params>
  struct Rank<const Kokkos::DynRankView<Array,Params...> >{
    static const int value = -1;
  };

  template<class Array, class ... Params, class Scalar>
  struct Return_Type<const Kokkos::DynRankView<Array,Params...>, Scalar>{
    typedef Kokkos::DynRankView<Array,Params...> ViewType;
    typedef typename ViewType::reference_type return_type;
    typedef typename ViewType::reference_type const_return_type;
  };
  
  template<class Array, class ... Params, class Scalar>
  struct Return_Type< Kokkos::DynRankView<Array,Params...>, Scalar>{
    typedef Kokkos::DynRankView<Array,Params...> ViewType;
    typedef typename ViewType::reference_type return_type;
    typedef typename ViewType::reference_type const_return_type;
  };

  template<class Array, class ... Params>
  struct CheckType<Kokkos::DynRankView<Array,Params...> >{
    static const bool value = true;
  };

  template<class Array, class ... Params>
  struct CheckType<const Kokkos::DynRankView<Array,Params...> >{
    static const bool value = true;
  };

  // for array wrapper below
  template<typename field>
  struct is_mdfield;

  template<class Array, class ... Params>
  struct is_mdfield<Kokkos::DynRankView<Array,Params...>>{ 
    enum {value=true};
  };


/*

template<class arg1, class arg2, class arg3, class arg4, class Scalar>
struct Return_Type<const Kokkos::View<arg1,arg2,arg3,arg4,Kokkos::Impl::ViewSpecializeSacadoFad>, Scalar>{
      typedef Kokkos::View<arg1,arg2,arg3,arg4,Kokkos::Impl::ViewSpecializeSacadoFad> ViewType;
      typedef typename ViewType::fad_view_type return_type;
      typedef typename ViewType::fad_view_type const_return_type;
};

template<class arg1, class arg2, class arg3, class arg4, class Scalar>
struct Return_Type< Kokkos::View<arg1,arg2,arg3,arg4,Kokkos::Impl::ViewSpecializeSacadoFad>, Scalar>{
      typedef Kokkos::View<arg1,arg2,arg3,arg4,Kokkos::Impl::ViewSpecializeSacadoFad> ViewType;
      typedef typename ViewType::fad_view_type return_type;
      typedef typename ViewType::fad_view_type const_return_type;
};
*/

template<class DataT,int leftrank>
struct RankSpec{};

template<class DataT>
struct RankSpec<DataT,-1>{
public:
    static index_type returnvalue(DataT& leftvalues){
        index_type dataRank   = leftvalues.rank();
        return dataRank;
}
};
template<class DataT>
struct RankSpec<DataT, 1>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 1;
}
};
template<class DataT>
struct RankSpec<DataT, 2>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 2;
}
};
template<class DataT>
struct RankSpec<DataT, 3>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 3;
}
};
template<class DataT>
struct RankSpec<DataT, 4>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 4;
}
};
template<class DataT>
struct RankSpec<DataT, 5>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 5;
}
};
template<class DataT>
struct RankSpec<DataT, 6>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 6;
}
};
template<class DataT>
struct RankSpec<DataT, 7>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 7;
}
};
template<class DataT>
struct RankSpec<DataT, 8>{
public:
    static index_type returnvalue(DataT& leftvalues){
        return 8;
}
};
template<class DataT>
inline index_type getrank(DataT& leftvalue){
    return RankSpec<DataT,Rank<DataT>::value>::returnvalue(leftvalue);

}

// ArrayType: Kokkos::View or Intrepid2::FieldContainer or Shards
//   container (implements operator(), and has view semantics)
// ArrayRank: -1 if run-time rank, else the compile-time rank
// isconstant: whether the array is constant
template<class Scalar,class ArrayType,int ArrayRank,bool isconstant>
struct ArrayWrapper;


// -1 means rank is determined at run time.
// It's like an Intrepid FieldContainer, not like a Kokkos::View.
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,-1, false> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
   InternalArrayType view;
     typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

        ArrayWrapper<Scalar,ArrayType,1,false>runtimewrapper1;
                ArrayWrapper<Scalar,ArrayType,2,false>runtimewrapper2;
                ArrayWrapper<Scalar,ArrayType,3,false>runtimewrapper3;
                ArrayWrapper<Scalar,ArrayType,4,false>runtimewrapper4;
                ArrayWrapper<Scalar,ArrayType,5,false>runtimewrapper5;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(ArrayType& view_):
     view(view_),
     runtimewrapper1(view_),
     runtimewrapper2(view_),
     runtimewrapper3(view_),
     runtimewrapper4(view_),
     runtimewrapper5(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
    index_type rank()const{
         return view.rank();
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0) const{
        return view(i0);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1) const{
        return view(i0,i1);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1, const index_type i2) const{
        return view(i0,i1,i2);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3) const{
        return view(i0,i1,i2,i3);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4) const{
        return view(i0,i1,i2,i3,i4);
}

};

template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,1,false> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
   InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 1;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1 = 0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0)  const{
    return view(i0);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,2,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 2;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0)  const{
    return view(i0,i1);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,3,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 3;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const{
    return view(i0,i1,i2);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,4,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

   typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 4;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0)const {
    return view(i0,i1,i2,i3);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,5,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 5;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1,i2,i3,i4);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,6,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 6;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0)  const{
    return view(i0,i1,i2,i3,i4,i5);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,7,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper( ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 7;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0)  const{
    return view(i0,i1,i2,i3,i4,i5,i6);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,8,false> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, ArrayType, ArrayType& >::type InternalArrayType;
  InternalArrayType view;

   typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 8;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const{
    return view(i0,i1,i2,i3,i4,i5,i6,i7);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,-1, true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;
  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

        ArrayWrapper<Scalar,ArrayType,1,true>runtimewrapper1;
                ArrayWrapper<Scalar,ArrayType,2,true>runtimewrapper2;
                ArrayWrapper<Scalar,ArrayType,3,true>runtimewrapper3;
                ArrayWrapper<Scalar,ArrayType,4,true>runtimewrapper4;
                ArrayWrapper<Scalar,ArrayType,5,true>runtimewrapper5;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):
     view(view_),
     runtimewrapper1(view_),
     runtimewrapper2(view_),
     runtimewrapper3(view_),
     runtimewrapper4(view_),
     runtimewrapper5(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return view.rank();
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0) const{
        return view(i0);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1) const{
        return view(i0,i1);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1, const index_type i2) const{
        return view(i0,i1,i2);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3 ) const{
        return view(i0,i1,i2,i3);
}
  KOKKOS_INLINE_FUNCTION
   rtype operator() (const index_type i0, const index_type i1, const index_type i2,
                          const index_type i3, const index_type i4) const{
        return view(i0,i1,i2,i3,i4);
}

};


template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,1,true> {
    typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 1;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1 = 0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,2,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 2;
  }
  KOKKOS_INLINE_FUNCTION

  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,3,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};

  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 3;
  }
  KOKKOS_INLINE_FUNCTION

  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1,i2);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,4,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 4;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1,i2,i3);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,5,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 5;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0)  const{
    return view(i0,i1,i2,i3,i4);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,6,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 6;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1,i2,i3,i4,i5);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,7,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 7;
  }

  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1,i2,i3,i4,i5,i6);
  }
};


template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,8,true> {
  typedef typename  std::conditional<is_mdfield<ArrayType>::value, const ArrayType, const ArrayType& >::type InternalArrayType;
  InternalArrayType view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  KOKKOS_INLINE_FUNCTION
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  KOKKOS_INLINE_FUNCTION
  index_type dimension(index_type i)const{
         return view.dimension(i);
  }
  index_type rank()const{
         return 8;
  }
  KOKKOS_INLINE_FUNCTION
  rtype operator() (const index_type i0, const index_type i1=0, const index_type i2 = 0,
                          const index_type i3 = 0, const index_type i4 = 0, const index_type i5 = 0,
                          const index_type i6 = 0, const index_type i7 = 0) const {
    return view(i0,i1,i2,i3,i4,i5,i6,i7);
  }
};

}
#endif

