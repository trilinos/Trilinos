// @HEADER
// ************************************************************************
//
//                           Intrepid Package
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

/** \file   Intrepid_RealSpaceToolsDef.hpp
    \brief  Definition file for utility classes providing basic linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_REALSPACETOOLS_DEF_HPP__
#define __INTREPID2_REALSPACETOOLS_DEF_HPP__

namespace Intrepid2 {

  template<typename ExecSpaceType>
  template<class ...absArrayProperties,
           class ...inArrayProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  absval( /**/  Kokkos::DynRankView<absArrayProperties...> absArray,
          const Kokkos::DynRankView<inArrayProperties...>   inArray ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inArray.rank() > 5,
                              ">>> ERROR (RealSpaceTools::absval): Input array container has rank larger than 5.");

    INTREPID2_TEST_FOR_ABORT( inArray.rank() != absArray.rank(),
                              ">>> ERROR (RealSpaceTools::absval): Array arguments must have identical ranks!");
    for (ordinal_type i=0;i<inArray.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inArray.dimension(i) != absArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::absval): Dimensions of array arguments do not agree!");
    }
#endif

    struct Functor {
      Kokkos::DynRankView<absArrayProperties...> _absArray;
      Kokkos::DynRankView<inArrayProperties...>  _inArray;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<absArrayProperties...> &absArray_,
              Kokkos::DynRankView<inArrayProperties...>  &inArray_)
        : _absArray(absArray_), _inArray(inArray_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type i) {
        const size_type jend = inArray.dimension(1);
        const size_type kend = inArray.dimension(2);
        const size_type lend = inArray.dimension(3);
        const size_type mend = inArray.dimension(4);

        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m)
                absArray(i,j,k,l,m) = Kokkos::abs(inArray(i,j,k,l,m));
      }
    };
    const size_type iend = inArray.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(absArray, inArray) );
  }

  template<typename ExecSpaceType>
  template<class ...inoutAbsArrayProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  absval( Kokkos::DynRankView<absArrayProperties...> inoutAbsArray ) {
    RealSpaceTools<ExecSpaceType>::absval(inoutAbsArray, inoutAbsArray);
  }

  template<typename ExecSpaceType>
  template<class ...inVecProperties>
  KOKKOS_INLINE_FUNCTION
  static typename Kokkos::DynRankView<inVecProperties...>::value_type
  vectorNorm( const Kokkos::DynRankView<inVecProperties...> inVec,
              const ENorm normType ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( !(inVec.rank() >= 1 && inVec.rank() <= 5),
                              ">>> ERROR (RealSpaceTools::vectorNorm): Vector argument must have rank 1!");
#endif
    typedef typename Kokkos::DynRankView<inVecProperties...>::value_type value_type;
    
    value_type norm(0);
    const size_type iend = inVec.dimension(0);
    const size_type jend = inVec.dimension(1);
    const size_type kend = inVec.dimension(2);
    const size_type lend = inVec.dimension(3);
    const size_type mend = inVec.dimension(4);
    switch(normType) {
    case NORM_TWO:{
      for (size_type i=0;j<iend;++i)
        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m)
                norm += inVec(i,j,k,l,m)*inVec(i,j,k,l,m);
      norm = Kokkos::sqrt(norm);
      break;
    }
    case NORM_INF:{
      for (size_type i=0;j<iend;++i)
        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m) {
                const value_type current = Kokkos::abs(inVec(i,j,k,l,m));
                norm = (norm < current ? current : norm);
              }
      break;
    }
    case NORM_ONE:{
      for (size_type i=0;j<iend;++i)
        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m) 
                norm += Kokkos::abs(inVec(i,j,k,l,m));
      break;
    }
    default: {
      INTREPID2_TEST_FOR_ABORT( ( (normType != NORM_TWO) &&
                                  (normType != NORM_INF) &&
                                  (normType != NORM_ONE) ),
                                ">>> ERROR (RealSpaceTools::vectorNorm): Invalid argument normType.");
    }
    }
    return norm;
  }
  
  template<typename ExecSpaceType>
  template<class ...normArrayProperties,
           class ...inVecProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  vectorNorm( /**/  Kokkos::DynRankView<normArrayProperties...> normArray,
              const Kokkos::DynRankView<inVecProperties...>     inVecs,
              const ENorm normType ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inVecs.rank() != (normArray.rank()+1),
                              ">>> ERROR (RealSpaceTools::vectorNorm): Ranks of norm and vector array arguments are incompatible!");
    INTREPID2_TEST_FOR_ABORT( inVecs.rank() < 2 || inVecs.rank() > 3,
                              ">>> ERROR (RealSpaceTools::vectorNorm): Rank of vector array must be 2 or 3!");
    for (size_type i=0;i<inVecs.rank()-1;++i) {
      INTREPID2_TEST_FOR_ABORT( inVecs.dimension(i) != normArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::vectorNorm): Dimensions of norm and vector arguments do not agree!");
    }
#endif

    struct Functor {
      Kokkos::DynRankView<normArrayProperties...> _normArray;
      Kokkos::DynRankView<inVecProperties...>     _inVecs;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<normArrayProperties...> &normArray_,
              Kokkos::DynRankView<inVecProperties...>     &inVecs_)
        : _normArray(normArray_), _inVecs(inVecs_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        // the rank of normArray is either 1 or 2
        const size_type i = iter % _normArray.dimension(0);
        const size_type j = iter / _normArray.dimension(0);

        auto vec = ( _inVecs.rank() == 2 ? Kokkos::subdynrankview(_inVecs, i,    Kokkos::ALL()) :
                     /**/                  Kokkos::subdynrankview(_inVecs, i, j, Kokkos::ALL()) );

        _normArray(i, j) = RealSpaceTools<Kokkos::Serial>::vectorNorm(vec, _normType);
      }
    };

    // normArray rank is either 1 or 2
    const size_type loopSize = normArray.dimension(0)*normArray.dimension(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(normArray, inVecs, normType) );
  }

  template<typename ExecSpaceType>
  template<class ...transposeMatProperties,
           class ...inMatProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>
  transpose( /**/  Kokkos::DynRankView<transposeMatProperties...> transposeMats,
             const Kokkos::DynRankView<inMatProperties...>        inMats ) {

    size_t arrayRank = getrank(inMats);
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inMats.rank() != transposeMats.rank(),
                              ">>> ERROR (RealSpaceTools::transpose): Matrix array arguments do not have identical ranks!");
    INTREPID2_TEST_FOR_ABORT( inMats.rank() < 2 || inMats.rank() > 4,
                              ">>> ERROR (RealSpaceTools::transpose): Rank of matrix array must be 2, 3, or 4!");
    for (size_type i=0;i<inMats.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inMats.dimension(i) != transposeMats.dimension(i),
                                ">>> ERROR (RealSpaceTools::transpose): Dimensions of matrix arguments do not agree!");
    }
    INTREPID2_TEST_FOR_ABORT( inMats.dimension(inMats.rank()-2) != inMats.dimension(inMats.rank()-1),
                              ">>> ERROR (RealSpaceTools::transpose): Matrices are not square!");
#endif

    struct Functor {
      Kokkos::DynRankView<transposeMatProperties...> _transposeMats;
      Kokkos::DynRankView<inMatProperties...>        _inMats;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<transposeMatProperties...> &transposeMats_,
              Kokkos::DynRankView<inMatProperties...>        &inMats_)
        : _transposeMats(transposeMats_), _inMats(inMats_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        // the rank of normArray is either 1 or 2
        const size_type i = iter % _transposeMats.dimension(0);
        const size_type j = iter / _transposeMats.dimension(0);

        const size_type r = _transposeMats.rank();
        auto dst = ( r == 2 ? Kokkos::subdynrankview(_transposeMats,       Kokkos::ALL(), Kokkos::ALL()) :
                     r == 3 ? Kokkos::subdynrankview(_transposeMats, i,    Kokkos::ALL(), Kokkos::ALL()) :
                     /**/     Kokkos::subdynrankview(_transposeMats, i, j, Kokkos::ALL(), Kokkos::ALL()) );

        auto src = ( r == 2 ? Kokkos::subdynrankview(_inMats,       Kokkos::ALL(), Kokkos::ALL()) :
                     r == 3 ? Kokkos::subdynrankview(_inMats, i,    Kokkos::ALL(), Kokkos::ALL()) :
                     /**/     Kokkos::subdynrankview(_inMats, i, j, Kokkos::ALL(), Kokkos::ALL()) );

        for (size_type i=0;i<src.dimension(0);++i) {
          dst(i,i) = src(i,i);
          for (size_type j=i+1;j<src.dimension(1);++j) {
            dst(i,j) = src(j,i);
            dst(j,i) = src(i,j);
          }
        }
      }
    };
    const size_type r = transposeMats.rank();
    const size_type loopSize = ( r == 2 ? 1 :
                                 r == 3 ? transposeMats.dimension(0) :
                                 /**/     transposeMats.dimension(0)*transposeMats.dimension(1) );

    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(transposeMats, inMats) );
  }

  template<typename ExecSpaceType>
  template<class ...inverseMatProperties,
           class ...inMatProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  inverse( /**/  Kokkos::DynRankView<inverseMatProperties...> inverseMats,
           const Kokkos::DynRankView<inMatProperties...>      inMats ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inMats.rank() != inverseMats.rank(),
                              ">>> ERROR (RealSpaceTools::inverse): Matrix array arguments do not have identical ranks!");
    INTREPID2_TEST_FOR_ABORT( inMats.rank() < 2 || inMats.rank() > 4,
                              ">>> ERROR (RealSpaceTools::inverse): Rank of matrix array must be 2, 3, or 4!");
    for (size_type i=0;i<inMats.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inMats.dimension(i) != inverseMats.dimension(i),
                                ">>> ERROR (RealSpaceTools::inverse): Dimensions of matrix arguments do not agree!");
    }
    INTREPID2_TEST_FOR_ABORT( inMats.dimension(inMats.rank()-2) != inMats.dimension(inMats.rank()-1),
                              ">>> ERROR (RealSpaceTools::inverse): Matrices are not square!");
    INTREPID2_TEST_FOR_ABORT( inMats.dimension(inMats.rank()-2) < 1 || inMats.dimension(inMats.rank()-2) > 3,
                              ">>> ERROR (RealSpaceTools::inverse): Spatial dimension must be 1, 2, or 3!");
#endif

    typedef typename Kokkos::DynRankView<inMatProperties...>::value_type value_type;
    struct Functor {
      /**/  Kokkos::DynRankView<inverseMatProperties...> _inverseMats;
      const Kokkos::DynRankView<inMatProperties...>      _inMats;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<inverseMatProperties...> &inverseMats_,
              Kokkos::DynRankView<inMatProperties...>      &inMats_)
        : _inverseMats(inverseMats_), _inMats(inMats_) {}
      
      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        const size_type i = iter % _inMats.dimension(0);
        const size_type j = iter / _inMats.dimension(0);
        
        const size_type r = _inMats.rank();
        auto mat = ( r == 2 ? Kokkos::subdynrankview(_inMats,       Kokkos::ALL(), Kokkos::ALL()) :
                     r == 3 ? Kokkos::subdynrankview(_inMats, i,    Kokkos::ALL(), Kokkos::ALL()) :
                     /**/     Kokkos::subdynrankview(_inLeft, i, j, Kokkos::ALL(), Kokkos::ALL()) );
        
        // compute determinant
        const value_type val = RealSpaceTools<Kokkos::Serial>::det(mat);
#ifdef HAVE_INTREPID_DEBUG
#ifdef HAVE_INTREPID_DEBUG_INF_CHECK
        INTREPID2_TEST_FOR_ABORT( val == 0,
                                  ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#endif
#endif
        switch (dim) {
        case 1: {
          inv(0,0) = value_type(1)/mat(0.0);
          break;
        }
        case 2: {
          inv(0,0) =   mat(1,1)/val;
          inv(1,1) =   mat(0,0)/val;
          
          inv(1,0) = - mat(1,0)/val;
          inv(0,1) = - mat(0,1)/val;
          break;
        }
        case 3: {
          value_type val0, val1, val2;

          val0 =   mat(1,1)*mat(2,2) - mat(2,1)*mat(1,2);
          val1 = - mat(1,0)*mat(2,2) + mat(2,0)*mat(1,2);
          val2 =   mat(1,0)*mat(2,1) - mat(2,0)*mat(1,1);

          inv(0,0) = val0/val;
          inv(1,0) = val1/val;
          inv(2,0) = val2/val;

          val0 =   mat(2,1)*mat(0,2) - mat(0,1)*mat(2,2);
          val1 =   mat(0,0)*mat(2,2) - mat(2,0)*mat(0,2);
          val2 = - mat(0,0)*mat(2,1) + mat(2,0)*mat(0,1);

          inv(0,1) = val0/val;
          inv(1,1) = val1/val;
          inv(2,1) = val2/val;

          val0 =   mat(0,1)*mat(1,2) - mat(1,1)*mat(0,2);
          val1 = - mat(0,0)*mat(1,2) + mat(1,0)*mat(0,2);
          val2 =   mat(0,0)*mat(1,1) - mat(1,0)*mat(0,1);

          inv(0,2) = val0/val;
          inv(1,2) = val1/val;
          inv(2,2) = val2/val;
          break;
        }
        }
      }
    };
    const size_type r = inMats.rank();
    const size_type loopSize = ( r == 2 ? 1 :
                                 r == 3 ? inverseMats.dimension(0) :
                                 /**/     inverseMats.dimension(0)*inverseMats.dimension(1) );
    
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(inverseMats, inMats) );
  }

  template<typename ExecSpaceType>
  template<class ...inMatProperties>
  KOKKOS_INLINE_FUNCTION
  static typename Kokkos::DynRankView<inMatProperties...>::value_type
  RealSpaceTools<ExecSpaceType>::
  det( const Kokkos::DynRankView<inMatProperties...> inMat) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inMat.rank() != 2,
                              ">>> ERROR (RealSpaceTools::det): Rank of matrix argument must be 2!");
    INTREPID2_TEST_FOR_ABORT( inMat.dimension(0) != inMat.dimension(1),
                              ">>> ERROR (RealSpaceTools::det): Matrix is not square!");
    INTREPID2_TEST_FOR_ABORT( inMat.dimension(0) < 1 || inMat.dimension(0) > 3,
                              ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif
    const size_type dim = inMat.dimension(0);
    return ( dim == 3 ? ( inMat(0,0) * inMat(1,1) * inMat(2,2) +
                          inMat(1,0) * inMat(2,1) * inMat(0,2) +
                          inMat(2,0) * inMat(0,1) * inMat(1,2) -
                          inMat(2,0) * inMat(1,1) * inMat(0,2) -
                          inMat(0,0) * inMat(2,1) * inMat(1,2) -
                          inMat(1,0) * inMat(0,1) * inMat(2,2) ) :
             dim == 2 ? ( inMat(0,0) * inMat(1,1) - 
                          inMat(0,1) * inMat(1,0) ) :
             /**/       ( inMat(0,0) ) );
  }


  template<typename ExecSpaceType>
  template<class ...detArrayProperties,
           class ...inMatProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  det( /**/  Kokkos::DynRankView<detArrayProperties...> detArray,
       const Kokkos::DynRankView<inMatProperties...>    inMats ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inMats.rank() != detArray.rank()+2,
                              ">>> ERROR (RealSpaceTools::det): Determinant and matrix array arguments do not have compatible ranks!");
    INTREPID2_TEST_FOR_ABORT( inMats.rank() < 3 || inMats.rank() > 4,
                              ">>> ERROR (RealSpaceTools::det): Rank of matrix array must be 3 or 4!");
    for (size_type i=0;i<inMats.rank()-2;++i) {
      INTREPID2_TEST_FOR_ABORT( inMats.dimension(i) != detArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::det): Dimensions of determinant and matrix array arguments do not agree!");
    }
    INTREPID2_TEST_FOR_ABORT( inMats.dimension(inMats.rank()-2) != inMats.dimension(inMats.rank()-1),
                              ">>> ERROR (RealSpaceTools::det): Matrices are not square!");
    INTREPID2_TEST_FOR_ABORT( inMats.dimension(inMats.rank()-2) < 1 || inMats.dimension(inMats.rank()-2) > 3,
                              ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#endif

    struct Functor {
      /**/  Kokkos::DynRankView<detArrayProperties...> _detArray;
      const Kokkos::DynRankView<inMatProperties...>    _inMats;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<detArrayProperties...> &detArray_,
              Kokkos::DynRankView<inMatProperties...>    &inMats_)
        : _detArray(detArray_), _inMats(inMats_) {}
      
      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        const size_type i = iter % _inMats.dimension(0);
        const size_type j = iter / _inMats.dimension(0);
        
        const size_type r = _inMats.rank();
        auto mat = ( r == 3 ? Kokkos::subdynrankview(_inMats, i,    Kokkos::ALL(), Kokkos::ALL()) :
                     /**/     Kokkos::subdynrankview(_inLeft, i, j, Kokkos::ALL(), Kokkos::ALL()) );
        
        _detArray(i, j) = RealSpaceTools<Kokkos::Serial>::det(mat);
      }
    };
    const size_type r = inMats.rank();
    const size_type loopSize = detArray.dimension(0)*detArray.dimension(1);
    
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(detArray, inMats) );
  }

  template<typename ExecSpaceType>
  template<class ...sumArrayProperties,
           class ...inArray1Properties,
           class ...inArray2Properties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  add( /**/  Kokkos::DynRankView<sumArrayProperties...> sumArray,
       const Kokkos::DynRankView<inArray1Properties...> inArray1,
       const Kokkos::DynRankView<inArray2Properties...> inArray2 ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inArray1.rank() != inArray2.rank() ||
                              inArray1.rank() != sumArray.rank(),
                              ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
    for (size_type i=0;i<inArray1.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inArray1.dimension(i) != inArray2.dimension(i) ||
                                inArray1.dimension(i) != sumArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
    }
#endif

    struct Functor {
      /**/  Kokkos::DynRankView<sumArrayProperties...> _sumArray;
      const Kokkos::DynRankView<inArray1Properties...> _inArray1;
      const Kokkos::DynRankView<inArray2Properties...> _inArray2;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<sumArrayProperties...> &sumArray_,
              Kokkos::DynRankView<inArray1Properties...> &inArray1_,
              Kokkos::DynRankView<inArray2Properties...> &inArray2_)
        : _sumArray(sumArray_), _inArray1(inArray1_), _inArray2(inArray2_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type i) {
        const size_type jend = _sumArray.dimension(1);
        const size_type kend = _sumArray.dimension(2);
        const size_type lend = _sumArray.dimension(3);
        const size_type mend = _sumArray.dimension(4);

        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m)
                _sumArray(i,j,k,l,m) = _inArray1(i,j,k,l,m) - _inArray2(i,j,k,l,m);
      }
    };
    const size_type iend = sumArray.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(sumArray, inArray1, inArray2) );
  }

  template<typename ExecSpaceType>
  template<class ...inoutSumArrayProperties,
           class ...inArrayProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  add( /**/  Kokkos::DynRankView<inoutSumArrayProperties...> inoutSumArray,
       const Kokkos::DynRankView<inArrayProperties...>       inArray ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inArray.rank() != inoutSumArray.rank(),
                              ">>> ERROR (RealSpaceTools::sum): Array arguments must have identical ranks!");
    for (size_type i=0;i<inArray.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inArray.dimension(i) != inoutSumArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::sum): Dimensions of array arguments do not agree!");
    }
#endif
    RealSpaceTools<ExecSpaceType>::add(inoutSumArray, inoutSumArray, inArray);
  }

  template<typename ExecSpaceType>
  template<class ...diffArrayProperties,
           class ...inArray1Properties,
           class ...inArray2Properties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  subtract( /**/  Kokkos::DynRankView<diffArrayProperties...> diffArray,
            const Kokkos::DynRankView<inArray1Properties...>  inArray1,
            const Kokkos::DynRankView<inArray2Properties...>  inArray2 ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inArray1.rank() != inArray2.rank() ||
                              inArray1.rank() != diffArray.rank(),
                              ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
    for (size_type i=0;i<inArray1.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inArray1.dimension(i) != inArray2.dimension(i) ||
                                inArray1.dimension(i) != diffArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
    }
#endif

    struct Functor {
      /**/  Kokkos::DynRankView<diffArrayProperties...> _diffArray;
      const Kokkos::DynRankView<inArray1Properties...>  _inArray1;
      const Kokkos::DynRankView<inArray2Properties...>  _inArray2;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<diffArrayProperties...> &diffArray_,
              Kokkos::DynRankView<inArray1Properties...>  &inArray1_,
              Kokkos::DynRankView<inArray2Properties...>  &inArray2_)
        : _diffArray(diffArray_), _inArray1(inArray1_), _inArray2(inArray2_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type i) {
        const size_type jend = _diffArray.dimension(1);
        const size_type kend = _diffArray.dimension(2);
        const size_type lend = _diffArray.dimension(3);
        const size_type mend = _diffArray.dimension(4);

        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m)
                _diffArray(i,j,k,l,m) = _inArray1(i,j,k,l,m) - _inArray2(i,j,k,l,m);
      }
    };
    const size_type iend = diffArray.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(diffArray, inArray1, inArray2) );
  }

  template<typename ExecSpaceType>
  template<class ...inoutDiffArrayProperties,
           class ...inArrayProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  subtract( /**/  Kokkos::DynRankView<inoutDiffArrayProperties...> inoutDiffArray,
            const Kokkos::DynRankView<inArrayProperties...>        inArray ) {
#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inArray.rank() != inoutDiffArray.rank(),
                              ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
    for (size_type i=0;i<inArray.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inArray.dimension(i) != inoutDiffArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
    }
#endif
    RealSpaceTools<ExecSpaceType>::subtract(inoutDiffArray, inoutDiffArray, inArray);
  }

  template<typename ExecSpaceType>
  template<typename ValueType,
           class ...scaledArrayProperties,
           class ...inArrayProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  scale( /**/  Kokkos::DynRankView<scaledArrayProperties...> scaledArray,
         const Kokkos::DynRankView<inArrayProperties...>     inArray,
         const ValueType alpha ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inArray.rank() > 5,
                              ">>> ERROR (RealSpaceTools::scale): Input array container has rank larger than 5.");
    INTREPID2_TEST_FOR_ABORT( inArray.rank() != scaledArray.rank,
                              ">>> ERROR (RealSpaceTools::scale): Array arguments must have identical ranks!");
    for (size_type i=0;i<inArray.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inArray.dimension(i) != scaledArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::scale): Dimensions of array arguments do not agree!");
    }
#endif

    struct Functor {
      /**/  Kokkos::DynRankView<scaledArrayProperties...> _scaledArray;
      const Kokkos::DynRankView<inArrayProperties...>     _inArray;
      ValueType _alpha;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<scaledArrayProperties...> &scaledArray_,
              Kokkos::DynRankView<inArrayProperties...>     &inArray_,
              const ValueType alpha_)
        : _scaledArray(scaledArray_), _inArray(inArray_), _alpha(alpha_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type i) {
        const size_type jend = inArray.dimension(1);
        const size_type kend = inArray.dimension(2);
        const size_type lend = inArray.dimension(3);
        const size_type mend = inArray.dimension(4);

        for (size_type j=0;j<jend;++j)
          for (size_type k=0;k<kend;++k)
            for (size_type l=0;l<lend;++l)
              for (size_type m=0;m<mend;++m)
                _scaledArray(i,j,k,l,m) = _alpha*_inArray(i,j,k,l,m);
      }
    };
    const size_type iend = scaledArray.dimension(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
    Kokkos::parallel_for( policy, Functor(scaledArray, inArray, alpha) );
  }

  template<typename ExecSpaceType>
  template<typename ValueType,
           class ...inoutScaledArrayProperties,
           class ...inArrayProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  scale( /**/  Kokkos::DynRankView<inoutScaledArrayProperties...> inoutScaledArray,
         const ValueType alpha ) {
    RealSpaceTools<ExecSpaceType>::scale(inoutScaledArray, inoutScaledArray, alpha);
  }

  template<typename ExecSpaceType>
  template<class ...inVec1Properties,
           class ...inVec2Properties>
  KOKKOS_INLINE_FUNCTION
  static typename Kokkos::DynRankView<inVec1Properties...>::value_type
  RealSpaceTools<ExecSpaceType>::
  dot( const Kokkos::DynRankView<inVec1Properties...> inVec1,
       const Kokkos::DynRankView<inVec2Properties...> inVec2 ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inVec1.rank != 1 || inVec2.rank() != 1,
                              ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
    INTREPID2_TEST_FOR_ABORT( inVec1.dimension(0) != inVec2.dimension(0),
                              ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
#endif
    typedef typename Kokkos::DynRankView<inVec1Properties...>::value_type value_type;

    // designed for small problems
    value_type r_val(0);

    const size_type iend = iVec1.dimension(0);
    for (size_type i=0;i<iend;++i)
      r_val += inVec1(i)*inVec2(i);

    return r_val;
  }

  template<typename ExecSpaceType>
  template<class ...dotArrayProperties,
           class ...inVec1Properties,
           class ...inVec2Properties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  dot( /**/  Kokkos::DynRankView<dotArrayProperties...> dotArray,
       const Kokkos::DynRankView<inVec1Properties...>   inVecs1,
       const Kokkos::DynRankView<inVec2Properties...>   inVecs2 ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inVecs1.rank() != (dotArray.rank()+1),
                              ">>> ERROR (RealSpaceTools::dot): Ranks of norm and vector array arguments are incompatible!");
    INTREPID2_TEST_FOR_ABORT( inVecs1.rank() != inVecs2.rank(),
                              ">>> ERROR (RealSpaceTools::dot): Ranks of input vector arguments must be identical!");
    INTREPID2_TEST_FOR_ABORT( inVecs1.rank() < 2 || inVecs1.rank() > 3,
                              ">>> ERROR (RealSpaceTools::dot): Rank of input vector arguments must be 2 or 3!");
    for (size_type i=0;i<inVecs1.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inVecs1.dimension(i) != inVecs2.dimension(i),
                                ">>> ERROR (RealSpaceTools::dot): Dimensions of input vector arguments do not agree!");
    }
    for (size_type i=0;i<inVecs1.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inVecs1.dimension(i) != dotArray.dimension(i),
                                ">>> ERROR (RealSpaceTools::dot): Dimensions of dot-product and vector arrays do not agree!");
    }
#endif

    struct Functor {
      /**/  Kokkos::DynRankView<dotArrayProperties...> _dotArray;
      const Kokkos::DynRankView<inVec1Properties...>   _inVecs1;
      const Kokkos::DynRankView<inVec2Properties...>   _inVecs2;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<dotArrayProperties...> &dotArray_,
              Kokkos::DynRankView<inVec1Properties...>   &inVecs1_,
              Kokkos::DynRankView<inVec2Properties...>   &inVecs2_)
        : _dotArray(dotArray_), _inVecs1(inVecs1_), _inVecs2(inVecs2_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        // the rank of normArray is either 1 or 2
        const size_type i = iter % _transposeMats.dimension(0);
        const size_type j = iter / _transposeMats.dimension(0);
        
        const size_type r = _inVec1.rank();
        auto vec1 = ( r == 2 ? Kokkos::subdynrankview(_inVec1, i,    Kokkos::ALL()) :
                      /**/     Kokkos::subdynrankview(_inVec1, i, j, Kokkos::ALL()) );
        auto vec2 = ( r == 2 ? Kokkos::subdynrankview(_inVec2, i,    Kokkos::ALL()) :
                      /**/     Kokkos::subdynrankview(_inVec2, i, j, Kokkos::ALL()) );
        
        _dotArray(i,j) = RealSpaceTools<Kokkos::Serial>::dot(vec1, vec2);
      }
    };
    const size_type r = dotArray..rank();
    const size_type loopSize = dotArray.dimension(0)*dotArray.dimension(1);
    
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(dotArray, inVecs1, inVecs2) );
  }

 
  template<typename ExecSpaceType>
  template<class ...matVecProperties,
           class ...inMatProperties,
           class ...inVecProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  matvec( /**/  Kokkos::DynRankView<matVecProperties...> matVecs,
          const Kokkos::DynRankView<inMatProperties...>  inMats,
          const Kokkos::DynRankView<inVecProperties...>  inVecs ) {

#ifdef HAVE_INTREPID_DEBUG
    INTREPID2_TEST_FOR_ABORT( inMats.rank() != (inVecs.rank()+1),
                              ">>> ERROR (RealSpaceTools::matvec): Vector and matrix array arguments do not have compatible ranks!");
    INTREPID2_TEST_FOR_ABORT( inMats.rank() < 2 || inMats.rank() > 4,
                              ">>> ERROR (RealSpaceTools::matvec): Rank of matrix array must be 2, 3 or 4!");
    INTREPID2_TEST_FOR_ABORT( matVecs.rank() != inVecs.rank(),
                              ">>> ERROR (RealSpaceTools::matvec): Vector arrays must be have the same rank!");
    for (size_type i=0;i<inMats.rank()-2;++i) {
      INTREPID2_TEST_FOR_ABORT( inMats.dimension(i) != inVecs.dimension(i),
                                ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector and matrix array arguments do not agree!");
    }
    for (size_type i=0;i<inVecs.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( matVecs.dimension(i) != inVecs.dimension(i),
                                ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector array arguments do not agree!");
    }
    INTREPID2_TEST_FOR_ABORT( inMats.dimension(inMats.rank()-1) != inVecs.dimension(inVecs.rank()-1),
                              ">>> ERROR (RealSpaceTools::matvec): Matrix column dimension does not match to the length of a vector!");
#endif

    struct Functor {
      /**/  Kokkos::DynRankView<matVecProperties...> _matVecs;
      const Kokkos::DynRankView<inMatProperties...>  _inMats;
      const Kokkos::DynRankView<inVecProperties...>  _inVecs;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<matVecProperties...> &matVecs_,
              Kokkos::DynRankView<inMatProperties...>  &inMats_,
              Kokkos::DynRankView<inVecProperties...>  &inVecs_) 
        : _matVecs(matVecs_), _inMats(inMats_), _inVecs(inVecs_) {}

      KOKKOS_INLINE_FUNCTION
      ~Functor = default;

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        const size_type i = iter % _inMats.dimension(0);
        const size_type j = iter / _inMats.dimension(0);

        const size_type rm = _inMats.rank();
        auto mat    = ( rm == 2 ? Kokkos::subdynrankview(_inMats,        Kokkos::ALL(), Kokkos::ALL()) : 
                        rm == 3 ? Kokkos::subdynrankview(_inMats,  i,    Kokkos::ALL(), Kokkos::ALL()) :
                        /**/      Kokkos::subdynrankview(_inMats,  i, j, Kokkos::ALL(), Kokkos::ALL()) );

        const size_type rv = _inVecs.rank();
        auto vec    = ( rv == 1 ? Kokkos::subdynrankview(_inVecs,        Kokkos::ALL()) : 
                        rv == 2 ? Kokkos::subdynrankview(_inVecs,  i,    Kokkos::ALL()) :
                        /**/      Kokkos::subdynrankview(_inVecs,  i, j, Kokkos::ALL()) );

        auto result = ( rv == 1 ? Kokkos::subdynrankview(_matVecs,       Kokkos::ALL()) :
                        rv == 2 ? Kokkos::subdynrankview(_matVecs, i,    Kokkos::ALL()) :
                        /**/      Kokkos::subdynrankview(_matVecs, i, j, Kokkos::ALL()) );
        
        const size_type iend = result.dimension(0);
        const size_type jend = vec.dimension(0);

        for (size_type i=0;i<iend;++i) {
          auto row = Kokkos::subdynrankview(mat, i, Kokkos::ALL());
          for (size_type j=0;j<jend;++j)
            result(i) = row(j)*vec(j);
        }
      }
    };
    const size_type r = inMats.rank();
    const size_type loopSize = ( r == 2 ? 1 :
                                 r == 3 ? inMats.dimension(0) :
                                 /**/     inMats.dimension(0)*inMats.dimension(1) );
    
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(matVecs, inMats, inVecs) );
  }


  template<typename ExecSpaceType>
  template<class ...vecProdProperties,
           class ...inLeftProperties,
           class ...inRightProperties>
  KOKKOS_INLINE_FUNCTION
  static void
  RealSpaceTools<ExecSpaceType>::
  vecprod( /**/  Kokkos::DynRankView<vecProdProperties...> vecProd,
           const Kokkos::DynRankView<inLeftProperties...>  inLeft,
           const Kokkos::DynRankView<inLeftProperties...>  inRight ) {
  
#ifdef HAVE_INTREPID_DEBUG
    /*
     *   Check array rank and spatial dimension range (if applicable)
     *      (1) all array arguments are required to have matching dimensions and rank: (D), (I0,D) or (I0,I1,D)
     *      (2) spatial dimension should be 2 or 3
     */
    // (1) check rank range on inLeft and then compare the other arrays with inLeft
    INTREPID2_TEST_FOR_ABORT( inLeft.rank() < 1 || inLeft.rank() > 3,
                              ">>> ERROR (RealSpaceTools::vecprod): Rank of matrix array must be 1, 2, or 3!");
    INTREPID2_TEST_FOR_ABORT( inLeft.rank() == inRight.rank(),
                              ">>> ERROR (RealSpaceTools::vecprod): Right and left arrays must be have the same rank!");
    INTREPID2_TEST_FOR_ABORT( inLeft.rank() == vecProd.rank(),
                              ">>> ERROR (RealSpaceTools::vecprod): Left and vecProd arrays must be have the same rank!");
    for (size_type i=0;i<inLeft.rank();++i) {
      INTREPID2_TEST_FOR_ABORT( inLeft.dimension(i) != inRight.dimension(i),
                                ">>> ERROR (RealSpaceTools::vecprod): Dimensions of matrix arguments do not agree!");
      INTREPID2_TEST_FOR_ABORT( inLeft.dimension(i) != vecProd.dimension(i),
                                ">>> ERROR (RealSpaceTools::vecprod): Dimensions of matrix arguments do not agree!");
    }

    // (2) spatial dimension ordinal = array rank - 1. Suffices to check only one array because we just
    //     checked whether or not the arrays have matching dimensions.
    INTREPID2_TEST_FOR_ABORT( inLeft.dimension(inLeft.rank()-1) < 2 || 
                              inLeft.dimension(inLeft.rank()-1) > 3,
                              ">>> ERROR (RealSpaceTools::vecprod): Dimensions of arrays (rank-1) must be 2 or 3!");                              
#endif

    const bool is_vecprod_3d = true;

    struct Functor {
      /**/  Kokkos::DynRankView<vecProdProperties...> _vecProd;
      const Kokkos::DynRankView<inLeftProperties...>  _inLeft;
      const Kokkos::DynRankView<inLeftProperties...>  _inRight;
      const bool _is_vecprod_3d;

      KOKKOS_INLINE_FUNCTION
      Functor(Kokkos::DynRankView<vecProdProperties...> &vecProd_,
              Kokkos::DynRankView<inLeftProperties...>  &inLeft_,
              Kokkos::DynRankView<inLeftProperties...>  &inRight_,
              const bool is_vecprod_3d_);
      : _vecProd(vecProd_), _inLeft(inLeft_), _inRight(inRight_), _is_vecprod_3d(is_vecprod_3d_) {}
      
      KOKKOS_INLINE_FUNCTION
      ~Functor = default;
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) {
        const size_type i = iter % _inLeft.dimension(0);
        const size_type j = iter / _inLeft.dimension(0);
        
        const size_type r = _inLeft.rank();

        auto left   = ( r == 1 ? Kokkos::subdynrankview(_inLeft,         Kokkos::ALL()) : 
                        r == 2 ? Kokkos::subdynrankview(_inLeft,   i,    Kokkos::ALL()) :
                        /**/     Kokkos::subdynrankview(_inLeft,   i, j, Kokkos::ALL()) );

        auto right  = ( r == 1 ? Kokkos::subdynrankview(_inRight,        Kokkos::ALL()) : 
                        r == 2 ? Kokkos::subdynrankview(_inRight,  i,    Kokkos::ALL()) :
                        /**/     Kokkos::subdynrankview(_inRight,  i, j, Kokkos::ALL()) );

        auto result = ( r == 1 ? Kokkos::subdynrankview(_vecProd,        Kokkos::ALL()) : 
                        r == 2 ? Kokkos::subdynrankview(_vecProd,  i,    Kokkos::ALL()) :
                        /**/     Kokkos::subdynrankview(_vecProd,  i, j, Kokkos::ALL()) );
        
        // branch prediction is possible 
        if (_is_vecprod_3d) {
          result(0) = left(1)*right(2) - left(2)*right(1);
          result(1) = left(2)*right(0) - left(0)*right(2);
          result(2) = left(0)*right(1) - left(1)*right(0);
        } else {
          result(0) = left(0)*right(1) - left(1)*right(0);
        }
      }
    };
    const size_type r = inLeft.rank();
    const size_type loopSize = ( r == 1 ? 1 :
                                 r == 2 ? inMats.dimension(0) :
                                 /**/     inMats.dimension(0)*inMats.dimension(1) );
    const bool is_vecprod_3d = (inLeft.dimension(inLeft.rank() - 1) == 3);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, Functor(vecProd, inLeft, inRight, is_vecprod_3d) );
  }

} // namespace Intrepid2

#endif





// =====================================
// Too much things...
// 

//   template<class ...inVec1Properties,
//            class ...inVec2Properties>
//   KOKKOS_INLINE_FUNCTION
//   static typename Kokkos::DynRankView<inVec1Properties...>::value_type
//   RealSpaceTools<ExecSpaceType>::
//   dot( const Kokkos::DynRankView<inVec1Properties...> inVec1,
//        const Kokkos::DynRankView<inVec2Properties...> inVec2 ) {

// #ifdef HAVE_INTREPID_DEBUG
//     INTREPID2_TEST_FOR_ABORT( inVec1.rank != 1 || inVec2.rank() != 1,
//                               ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
//     INTREPID2_TEST_FOR_ABORT( inVec1.dimension(0) != inVec2.dimension(0),
//                               ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
// #endif
//     typedef typename Kokkos::DynRankView<inVec1Properties...>::value_type value_type;

//     // designed for small problems
//     value_type r_val(0);

//     // ** This is probably not necessary 
//     if (Kokkos::Impl::is_same<ExecSpace,Kokkos::Serial>::value) {
//       const size_type iend = iVec1.dimension(0);
//       for (size_type i=0;i<iend;++i)
//         r_val += inVec1(i)*inVec2(i);
//     } else {
//       struct Functor {
//         Kokkos::DynRankView<inVec1Properties...> _inVec1;
//         Kokkos::DynRankView<inVec2Properties...> _inVec2;
        
//         KOKKOS_INLINE_FUNCTION
//         Functor(const Kokkos::DynRankView<inVec1Properties...> &inVec1_,
//                 const Kokkos::DynRankView<inVec2Properties...> &inVec2_);
        
//         KOKKOS_INLINE_FUNCTION
//         ~Functor = default;
        
//         KOKKOS_INLINE_FUNCTION
//         void operator()(size_type i, value_type &dst ) const {
//           dst += _inVec1(i)*_inVec2(i);
//         }
        
//         KOKKOS_INLINE_FUNCTION
//         void join(volatile value_type &dst ,
//                   const volatile value_type &src) const {
//           dst += src;
//         }
//       };     
//       const size_type iend = inVec1.dimension(0);
//       Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
//       Kokkos::parallel_for( policy, Functor(inVec1, inVec2), r_val );
//     }
//     return r_val;
//   }


