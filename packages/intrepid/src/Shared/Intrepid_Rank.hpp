#ifndef INTREPID_RANK_HPP
#define INTREPID_RANK_HPP

template<class A>
struct CheckType{static const bool value = false; };


template<class A>
struct Rank{static const int value = -1;};

template<class A,class Scalar>
struct Return_Type{
    typedef Scalar& return_type;
    typedef Scalar const_return_type;
    };

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

template<class DataT,int leftrank>
struct RankSpec{};

template<class DataT>
struct RankSpec<DataT,-1>{
public:
    static int returnvalue(DataT& leftvalues){
        int dataRank   = leftvalues.rank();
        return dataRank;
}
};
template<class DataT>
struct RankSpec<DataT, 1>{
public:
    static int returnvalue(DataT& leftvalues){
        return 1;
}
};
template<class DataT>
struct RankSpec<DataT, 2>{
public:
    static int returnvalue(DataT& leftvalues){
        return 2;
}
};
template<class DataT>
struct RankSpec<DataT, 3>{
public:
    static int returnvalue(DataT& leftvalues){
        return 3;
}
};
template<class DataT>
struct RankSpec<DataT, 4>{
public:
    static int returnvalue(DataT& leftvalues){
        return 4;
}
};
template<class DataT>
struct RankSpec<DataT, 5>{
public:
    static int returnvalue(DataT& leftvalues){
        return 5;
}
};
template<class DataT>
struct RankSpec<DataT, 6>{
public:
    static int returnvalue(DataT& leftvalues){
        return 6;
}
};
template<class DataT>
struct RankSpec<DataT, 7>{
public:
    static int returnvalue(DataT& leftvalues){
        return 7;
}
};
template<class DataT>
struct RankSpec<DataT, 8>{
public:
    static int returnvalue(DataT& leftvalues){
        return 8;
}
};

template<class DataT>
inline size_t getrank(DataT& leftvalue){
    return RankSpec<DataT,Rank<DataT>::value>::returnvalue(leftvalue);

}

// ArrayType: Intrepid::FieldContainer or Shards
//   container (implements operator(), and has view semantics)
// ArrayRank: -1 if run-time rank, else the compile-time rank
// isconstant: whether the array is constant
template<class Scalar,class ArrayType,int ArrayRank,bool isconstant>
struct ArrayWrapper;

// -1 means rank is determined at run time.
// It's like an Intrepid FieldContainer.
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,-1, false> {
   ArrayType& view;
     typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

        ArrayWrapper<Scalar,ArrayType,1,false>runtimewrapper1;
                ArrayWrapper<Scalar,ArrayType,2,false>runtimewrapper2;
                ArrayWrapper<Scalar,ArrayType,3,false>runtimewrapper3;
                ArrayWrapper<Scalar,ArrayType,4,false>runtimewrapper4;
                ArrayWrapper<Scalar,ArrayType,5,false>runtimewrapper5;

  ArrayWrapper( ArrayType& view_):
     view(view_),
     runtimewrapper1(view_),
     runtimewrapper2(view_),
     runtimewrapper3(view_),
     runtimewrapper4(view_),
     runtimewrapper5(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
    int rank()const{
         return view.rank();
  }

  rtype operator() (const int i0) const{
        return view(i0);
  }

   rtype operator() (const int i0, const int i1) const{
        return view(i0,i1);
  }

   rtype operator() (const int i0, const int i1, const int i2) const{
        return view(i0,i1,i2);
  }

   rtype operator() (const int i0, const int i1, const int i2,
                          const int i3) const{
        return view(i0,i1,i2,i3);
}

   rtype operator() (const int i0, const int i1, const int i2,
                          const int i3, const int i4) const{
        return view(i0,i1,i2,i3,i4);
}

};

template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,1,false> {
   ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 1;
  }

  rtype operator() (const int i0, const int i1 = 0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0)  const{
    return view(i0);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,2,false> {
   ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 2;
  }
  
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0)  const{
    return view(i0,i1);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,3,false> {
   ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 3;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const{
    return view(i0,i1,i2);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,4,false> {
   ArrayType& view;

   typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 4;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0)const {
    return view(i0,i1,i2,i3);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,5,false> {
   ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 5;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1,i2,i3,i4);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,6,false> {
   ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 6;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0)  const{
    return view(i0,i1,i2,i3,i4,i5);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,7,false> {
   ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper( ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 7;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0)  const{
    return view(i0,i1,i2,i3,i4,i5,i6);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,8,false> {
   ArrayType& view;

   typedef typename Return_Type<ArrayType, Scalar>::return_type rtype;

  ArrayWrapper(const ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 8;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const{
    return view(i0,i1,i2,i3,i4,i5,i6,i7);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,-1, true> {
  const ArrayType& view;
  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

        ArrayWrapper<Scalar,ArrayType,1,true>runtimewrapper1;
                ArrayWrapper<Scalar,ArrayType,2,true>runtimewrapper2;
                ArrayWrapper<Scalar,ArrayType,3,true>runtimewrapper3;
                ArrayWrapper<Scalar,ArrayType,4,true>runtimewrapper4;
                ArrayWrapper<Scalar,ArrayType,5,true>runtimewrapper5;

  ArrayWrapper(const ArrayType& view_):
     view(view_),
     runtimewrapper1(view_),
     runtimewrapper2(view_),
     runtimewrapper3(view_),
     runtimewrapper4(view_),
     runtimewrapper5(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return view.rank();
  }
  rtype operator() (const int i0) const{
        return view(i0);
}
   rtype operator() (const int i0, const int i1) const{
        return view(i0,i1);
}
   rtype operator() (const int i0, const int i1, const int i2) const{
        return view(i0,i1,i2);
}
   rtype operator() (const int i0, const int i1, const int i2,
                          const int i3 ) const{
        return view(i0,i1,i2,i3);
}
   rtype operator() (const int i0, const int i1, const int i2,
                          const int i3, const int i4) const{
        return view(i0,i1,i2,i3,i4);
}

};


template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,1,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  ArrayWrapper(const ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 1;
  }
  rtype operator() (const int i0, const int i1 = 0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,2,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  ArrayWrapper(const ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 2;
  }

  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1);
  }

};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar, ArrayType,3,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  ArrayWrapper(const ArrayType& view_):view(view_) {};

  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 3;
  }

  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1,i2);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,4,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;

  ArrayWrapper(const ArrayType& view_):view(view_) {};
  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 4;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1,i2,i3);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,5,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 5;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0)  const{
    return view(i0,i1,i2,i3,i4);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,6,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 6;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1,i2,i3,i4,i5);
  }
};
template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,7,true> {
  const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 7;
  }

  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1,i2,i3,i4,i5,i6);
  }
};


template<class Scalar,class ArrayType>
struct ArrayWrapper<Scalar,ArrayType,8,true> {
   const ArrayType& view;

  typedef typename Return_Type<ArrayType, Scalar>::const_return_type rtype;
  ArrayWrapper(const ArrayType& view_):view(view_) {};
  int dimension(int i)const{
         return view.dimension(i);
  }
  int rank()const{
         return 8;
  }
  rtype operator() (const int i0, const int i1=0, const int i2 = 0,
                          const int i3 = 0, const int i4 = 0, const int i5 = 0,
                          const int i6 = 0, const int i7 = 0) const {
    return view(i0,i1,i2,i3,i4,i5,i6,i7);
  }
};

#endif


#if defined(Intrepid_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The Intrepid package is deprecated"
#endif
#endif

