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

/** \file   Intrepid2_RealSpaceToolsDef.hpp
    \brief  Definition file for utility classes providing basic linear algebra functionality.
    \author Created by P. Bochev, D. Ridzal, and D. Day.
    Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_REALSPACETOOLS_DEF_HPP__
#define __INTREPID2_REALSPACETOOLS_DEF_HPP__

namespace Intrepid2 {

  // ------------------------------------------------------------------------------------

  //
  // serial version code
  //

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename inVecValueType, class ...inVecProperties>
  KOKKOS_INLINE_FUNCTION
  inVecValueType
  RealSpaceTools<SpT>::Serial::
  vectorNorm( const Kokkos::DynRankView<inVecValueType,inVecProperties...> inVec,
              const ENorm normType ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( !(inVec.rank() >= 1 && inVec.rank() <= 5), dbgInfo,
                                      ">>> ERROR (RealSpaceTools::vectorNorm): Vector argument must have rank 1!" );
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return inVecValueType(0);
#endif
    }
#endif
    typedef inVecValueType value_type;

    value_type norm(0);

    const ordinal_type iend = inVec.extent(0);
    const ordinal_type jend = inVec.extent(1);
    const ordinal_type kend = inVec.extent(2);
    const ordinal_type lend = inVec.extent(3);
    const ordinal_type mend = inVec.extent(4);
    switch(normType) {
    case NORM_TWO:{
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                norm += inVec.access(i,j,k,l,m)*inVec.access(i,j,k,l,m);
      norm = sqrt(norm);
      break;
    }
    case NORM_INF:{
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m) {
                const value_type current = Util<value_type>::abs(inVec.access(i,j,k,l,m));
                norm = (norm < current ? current : norm);
              }
      break;
    }
    case NORM_ONE:{
      for (ordinal_type i=0;i<iend;++i)
        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                norm += Util<value_type>::abs(inVec.access(i,j,k,l,m));
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

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename inMatValueType, class ...inMatProperties>
  KOKKOS_INLINE_FUNCTION
  inMatValueType
  RealSpaceTools<SpT>::Serial::
  det( const Kokkos::DynRankView<inMatValueType,inMatProperties...> inMat) {

    typedef typename decltype(inMat)::non_const_value_type value_type;
#ifdef HAVE_INTREPID2_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inMat.rank() != 2, dbgInfo,
                                      ">>> ERROR (RealSpaceTools::det): Rank of matrix argument must be 2!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inMat.extent(0) != inMat.extent(1), dbgInfo,
                                      ">>> ERROR (RealSpaceTools::det): Matrix is not square!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inMat.extent(0) < 1 || inMat.extent(0) > 3, dbgInfo,
                                      ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return value_type(0);
#endif
    }
#endif
    const auto dim = inMat.extent(0);
    
    value_type r_val = 0.0;
    switch (dim) {
    case 3:
      r_val = ( inMat(0,0) * inMat(1,1) * inMat(2,2) +
                inMat(1,0) * inMat(2,1) * inMat(0,2) +
                inMat(2,0) * inMat(0,1) * inMat(1,2) -
                inMat(2,0) * inMat(1,1) * inMat(0,2) -
                inMat(0,0) * inMat(2,1) * inMat(1,2) -
                inMat(1,0) * inMat(0,1) * inMat(2,2) );
      break;
    case 2:
      r_val = ( inMat(0,0) * inMat(1,1) -
                inMat(0,1) * inMat(1,0) );
      break;
    case 1:
      r_val = ( inMat(0,0) );
      break;
    }
    return r_val;
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename inVec1ValueType, class ...inVec1Properties,
           typename inVec2ValueType, class ...inVec2Properties>
  KOKKOS_INLINE_FUNCTION
  inVec1ValueType
  RealSpaceTools<SpT>::Serial::
  dot( const Kokkos::DynRankView<inVec1ValueType,inVec1Properties...> inVec1,
       const Kokkos::DynRankView<inVec2ValueType,inVec2Properties...> inVec2 ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      bool dbgInfo = false;
      INTREPID2_TEST_FOR_DEBUG_ABORT( inVec1.rank() != 1 || inVec2.rank() != 1, dbgInfo,
                                      ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
      INTREPID2_TEST_FOR_DEBUG_ABORT( inVec1.extent(0) != inVec2.extent(0), dbgInfo,
                                      ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
#ifdef INTREPID2_TEST_FOR_EXCEPTION_OVERRIDE_TO_CONTINUE
      if (dbgInfo) return inVec1ValueType(0);
#endif
    }
#endif
    typedef inVec1ValueType value_type;

    // designed for small problems
    value_type r_val(0);

    const ordinal_type iend = inVec1.extent(0);
    for (ordinal_type i=0;i<iend;++i)
      r_val += inVec1(i)*inVec2(i);

    return r_val;
  }

  // ------------------------------------------------------------------------------------

  //
  // use parallel for
  //

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {

    /**
      \brief Functor for extractScalarValues see Intrepid2::RealSpaceTools for more
    */ 
    template<typename outputViewType,
             typename inputViewType>
    struct F_extractScalarValues {
      outputViewType _output;
      inputViewType  _input;

      KOKKOS_INLINE_FUNCTION
      F_extractScalarValues( outputViewType output_,
                       inputViewType  input_ )
        : _output(output_), _input(input_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i) const {
        const ordinal_type jend = _input.extent(1);
        const ordinal_type kend = _input.extent(2);
        const ordinal_type lend = _input.extent(3);
        const ordinal_type mend = _input.extent(4);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                _output.access(i,j,k,l,m) = get_scalar_value(_input.access(i,j,k,l,m));
      }
    };
  }

  template<typename SpT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void
  RealSpaceTools<SpT>::
  extractScalarValues(       Kokkos::DynRankView<outputValueType,outputProperties...>  output,
                       const Kokkos::DynRankView<inputValueType, inputProperties...>   input ) {
    typedef          Kokkos::DynRankView<outputValueType,outputProperties...> outputViewType;
    typedef          Kokkos::DynRankView<inputValueType,inputProperties...> inputViewType;
    typedef          FunctorRealSpaceTools::F_extractScalarValues<outputViewType,inputViewType> FunctorType;
    typedef typename ExecSpace<typename inputViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;
    
    const auto loopSize = input.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, input) );
  }

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor for clone see Intrepid2::RealSpaceTools for more
    */ 
    template<typename outputViewType,
             typename inputViewType>
    struct F_clone {
      outputViewType _output;
      inputViewType  _input;

      KOKKOS_INLINE_FUNCTION
      F_clone( outputViewType output_,
               inputViewType  input_ )
        : _output(output_), _input(input_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type iter) const {
        ordinal_type k0(0), k1(0), k2(0);
        const auto rankDiff = _output.rank() - _input.rank();
        switch (rankDiff) {
        case 3:
          unrollIndex( k0, k1, k2,
                             _output.extent(0),
                             _output.extent(1),
                             _output.extent(2),
                             iter );
          break;
        case 2:
          unrollIndex( k0, k1, 
                             _output.extent(0),
                             _output.extent(1),
                             iter );
          break;
        case 1:
          k0 = iter;
          break;
        }
        
        auto out = (rankDiff == 3 ? Kokkos::subview(_output, k0, k1, k2, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) :
                    rankDiff == 2 ? Kokkos::subview(_output, k0, k1,     Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) :
                    rankDiff == 1 ? Kokkos::subview(_output, k0,         Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()) :
                                    Kokkos::subview(_output,             Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL()));
        const ordinal_type iend = _input.extent(0);
        const ordinal_type jend = _input.extent(1);
        const ordinal_type kend = _input.extent(2);
        for (ordinal_type i=0;i<iend;++i)
          for (ordinal_type j=0;j<jend;++j)
            for (ordinal_type k=0;k<kend;++k)
              out.access(i,j,k) = _input.access(i,j,k);
      }
    };
  }

  template<typename SpT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void
  RealSpaceTools<SpT>::
  clone(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
         const Kokkos::DynRankView<inputValueType,inputProperties...> input ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      // input has at most rank 3
      INTREPID2_TEST_FOR_EXCEPTION( input.rank() > 3, std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::clone): Input container has rank larger than 3.");


      INTREPID2_TEST_FOR_EXCEPTION( input.rank() > output.rank(), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::clone): Input container rank should be smaller than ouput rank.");

      const size_type rankDiff = output.rank() - input.rank();

      // Difference of output and input rank is at most 3.
      INTREPID2_TEST_FOR_EXCEPTION( rankDiff > 3, std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::clone): Difference between output container rank and input container rank is larger than 3.");


      for (size_type i=0;i<input.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( input.extent(i) != output.extent(rankDiff + i), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::clone): Dimensions of array arguments do not agree!");
      }
    }
#endif
    typedef          Kokkos::DynRankView<outputValueType,outputProperties...>     outputViewType;
    typedef          Kokkos::DynRankView<inputValueType,inputProperties...>       inputViewType;
    typedef          FunctorRealSpaceTools::F_clone<outputViewType,inputViewType> FunctorType;
    typedef typename ExecSpace<typename inputViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    size_type loopSize = 1;
    const ordinal_type out_rank = output.rank();
    const ordinal_type in_rank = input.rank();
    const ordinal_type rankDiff = out_rank - in_rank;
    for (ordinal_type i=0;i<rankDiff;++i)
      loopSize *= output.extent(i);

    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(output, input) );
  }

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute absolute value see Intrepid2::RealSpaceTools for more
    */ 
    template<typename absArrayViewType,
             typename inArrayViewType>
    struct F_absval {
      absArrayViewType _absArray;
      inArrayViewType  _inArray;

      KOKKOS_INLINE_FUNCTION
      F_absval( absArrayViewType absArray_,
                inArrayViewType  inArray_ )
        : _absArray(absArray_), _inArray(inArray_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i) const {
        const ordinal_type jend = _inArray.extent(1);
        const ordinal_type kend = _inArray.extent(2);
        const ordinal_type lend = _inArray.extent(3);
        const ordinal_type mend = _inArray.extent(4);

        typedef typename inArrayViewType::non_const_value_type value_type;

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                _absArray.access(i,j,k,l,m) = Util<value_type>::abs(_inArray.access(i,j,k,l,m));
      }
    };
  }

  template<typename SpT>
  template<typename absArrayValueType, class ...absArrayProperties,
           typename inArrayValueType,  class ...inArrayProperties>
  void
  RealSpaceTools<SpT>::
  absval(       Kokkos::DynRankView<absArrayValueType,absArrayProperties...> absArray,
          const Kokkos::DynRankView<inArrayValueType, inArrayProperties...>   inArray ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inArray.rank() > 5, std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::absval): Input array container has rank larger than 5.");

      INTREPID2_TEST_FOR_EXCEPTION( inArray.rank() != absArray.rank(), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::absval): Array arguments must have identical ranks!");
      for (size_type i=0;i<inArray.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inArray.extent(i) != absArray.extent(i), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::absval): Dimensions of array arguments do not agree!");
      }
    }
#endif

    typedef          Kokkos::DynRankView<absArrayValueType,absArrayProperties...>            absArrayViewType;
    typedef          Kokkos::DynRankView<inArrayValueType, inArrayProperties...>             inArrayViewType;
    typedef          FunctorRealSpaceTools::F_absval<absArrayViewType,inArrayViewType>       FunctorType;
    typedef typename ExecSpace<typename inArrayViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    const auto loopSize = inArray.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(absArray, inArray) );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename inoutArrayValueType, class ...inoutAbsArrayProperties>
  void
  RealSpaceTools<SpT>::
  absval( Kokkos::DynRankView<inoutArrayValueType,inoutAbsArrayProperties...> inoutAbsArray ) {
    RealSpaceTools<SpT>::absval(inoutAbsArray, inoutAbsArray);
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute vector norm see Intrepid2::RealSpaceTools for more
    */ 
    template<typename normArrayViewType,
             typename inVecViewType>
    struct F_vectorNorm {
      normArrayViewType _normArray;
      inVecViewType     _inVecs;
      const ENorm _normType;

      KOKKOS_INLINE_FUNCTION
      F_vectorNorm( normArrayViewType normArray_,
                    inVecViewType     inVecs_,
                    const ENorm       normType_ )
        : _normArray(normArray_), _inVecs(inVecs_), _normType(normType_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type iter) const {
        ordinal_type i, j;
        unrollIndex( i, j,
                           _normArray.extent(0),
                           _normArray.extent(1),
                           iter );

        auto vec = ( _inVecs.rank() == 2 ? Kokkos::subview(_inVecs, i,    Kokkos::ALL()) :
                                           Kokkos::subview(_inVecs, i, j, Kokkos::ALL()) );

        _normArray(i, j) = RealSpaceTools<>::Serial::vectorNorm(vec, _normType);
      }
    };
  }

  template<typename SpT>
  template<typename normArrayValueType, class ...normArrayProperties,
           typename inVecValueType,     class ...inVecProperties>
  void
  RealSpaceTools<SpT>::
  vectorNorm(       Kokkos::DynRankView<normArrayValueType,normArrayProperties...> normArray,
              const Kokkos::DynRankView<inVecValueType,    inVecProperties...>     inVecs,
              const ENorm normType ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inVecs.rank() != (normArray.rank()+1), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::vectorNorm): Ranks of norm and vector array arguments are incompatible!");
      INTREPID2_TEST_FOR_EXCEPTION( inVecs.rank() < 2 || inVecs.rank() > 3, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::vectorNorm): Rank of vector array must be 2 or 3!");
      for (size_type i=0;i<inVecs.rank()-1;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inVecs.extent(i) != normArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::vectorNorm): Dimensions of norm and vector arguments do not agree!");
      }
    }
#endif

    typedef          Kokkos::DynRankView<normArrayValueType,normArrayProperties...>        normArrayViewType;
    typedef          Kokkos::DynRankView<inVecValueType,    inVecProperties...>            inVecViewType;
    typedef          FunctorRealSpaceTools::F_vectorNorm<normArrayViewType,inVecViewType>  FunctorType;
    typedef typename ExecSpace<typename inVecViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    // normArray rank is either 1 or 2
    const auto loopSize = normArray.extent(0)*normArray.extent(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(normArray, inVecs, normType) );
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute transpose see Intrepid2::RealSpaceTools for more
    */ 
    template<typename transposeMatViewType,
             typename inMatViewType>
    struct F_transpose {
      transposeMatViewType _transposeMats;
      inMatViewType        _inMats;

      KOKKOS_INLINE_FUNCTION
      F_transpose( transposeMatViewType transposeMats_,
                   inMatViewType        inMats_)
        : _transposeMats(transposeMats_), _inMats(inMats_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        // the rank of normArray is either 1 or 2
        const auto r = _transposeMats.rank();
        ordinal_type _i = iter, _j = 0;

        if ( r > 3 )
          unrollIndex( _i, _j,
                             _transposeMats.extent(0),
                             _transposeMats.extent(1),
                             iter );

        auto dst = ( r == 2 ? Kokkos::subview(_transposeMats,       Kokkos::ALL(), Kokkos::ALL()) :
                     r == 3 ? Kokkos::subview(_transposeMats, _i,    Kokkos::ALL(), Kokkos::ALL()) :
                              Kokkos::subview(_transposeMats, _i, _j, Kokkos::ALL(), Kokkos::ALL()) );

        auto src = ( r == 2 ? Kokkos::subview(_inMats,       Kokkos::ALL(), Kokkos::ALL()) :
                     r == 3 ? Kokkos::subview(_inMats, _i,    Kokkos::ALL(), Kokkos::ALL()) :
                              Kokkos::subview(_inMats, _i, _j, Kokkos::ALL(), Kokkos::ALL()) );

        for (size_type i=0;i<src.extent(0);++i) {
          dst(i, i) = src(i, i);
          for (size_type j=i+1;j<src.extent(1);++j) {
            dst(i, j) = src(j, i);
            dst(j, i) = src(i, j);
          }
        }
      }
    };
  }

  template<typename SpT>
  template<typename transposeMatValueType, class ...transposeMatProperties,
           typename inMatValueType,        class ...inMatProperties>
  void
  RealSpaceTools<SpT>::
  transpose(       Kokkos::DynRankView<transposeMatValueType,transposeMatProperties...> transposeMats,
             const Kokkos::DynRankView<inMatValueType,       inMatProperties...>        inMats ) {

#ifdef HAVE_INTREPID2_DEBUG
    {

      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() != transposeMats.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::transpose): Matrix array arguments do not have identical ranks!");
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() < 2 || inMats.rank() > 4, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::transpose): Rank of matrix array must be 2, 3, or 4!");
      for (size_type i=0;i<inMats.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(i) != transposeMats.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::transpose): Dimensions of matrix arguments do not agree!");
      }
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-2) != inMats.extent(inMats.rank()-1), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::transpose): Matrices are not square!");
    }
#endif

    typedef          Kokkos::DynRankView<transposeMatValueType,transposeMatProperties...>   transposeMatViewType;
    typedef          Kokkos::DynRankView<inMatValueType,       inMatProperties...>          inMatViewType;
    typedef          FunctorRealSpaceTools::F_transpose<transposeMatViewType,inMatViewType> FunctorType;
    typedef typename ExecSpace<typename inMatViewType::execution_space,SpT>::ExecSpaceType  ExecSpaceType;

    const auto r = transposeMats.rank();
    const auto loopSize = ( r == 2 ? 1 :
                            r == 3 ? transposeMats.extent(0) :
                                     transposeMats.extent(0)*transposeMats.extent(1) );

    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(transposeMats, inMats) );
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute inverse see Intrepid2::RealSpaceTools for more
    */ 
    template<typename inverseMatViewType,
             typename inMatViewType>
    struct F_inverse {
      typedef typename inMatViewType::non_const_value_type value_type;
      inverseMatViewType _inverseMats;
      inMatViewType      _inMats;

      KOKKOS_INLINE_FUNCTION
      F_inverse( inverseMatViewType inverseMats_,
                 inMatViewType      inMats_ )
        : _inverseMats(inverseMats_), _inMats(inMats_) {}

      template<typename matViewType,
               typename invViewType>
      KOKKOS_FORCEINLINE_FUNCTION
      static void 
      apply_inverse(       invViewType inv,
                     const matViewType mat ) {
        // compute determinant
        const value_type val = RealSpaceTools<>::Serial::det(mat);
        
#ifdef HAVE_INTREPID2_DEBUG
        {
#ifdef HAVE_INTREPID2_DEBUG_INF_CHECK
          bool dbgInfo = false;
          INTREPID2_TEST_FOR_DEBUG_ABORT( val == 0, dbgInfo,
                                          ">>> ERROR (Matrix): Inverse of a singular matrix is undefined!");
#ifdef INTREPID2_TEST_FOR_DEBUG_ABORT_OVERRIDE_TO_CONTINUE
          if (dbgInfo) return;
#endif
#endif
        }
#endif
        switch (mat.extent(0)) {
        case 1: {
          inv(0,0) = value_type(1)/mat(0,0);
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
      
      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl,
                      const ordinal_type pt) const {
        const auto mat = Kokkos::subview(_inMats,      cl, pt, Kokkos::ALL(), Kokkos::ALL());
        auto       inv = Kokkos::subview(_inverseMats, cl, pt, Kokkos::ALL(), Kokkos::ALL());

        apply_inverse( inv, mat );
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type pt) const {
        const auto mat = Kokkos::subview(_inMats,      pt, Kokkos::ALL(), Kokkos::ALL());
        auto       inv = Kokkos::subview(_inverseMats, pt, Kokkos::ALL(), Kokkos::ALL());

        apply_inverse( inv, mat );
      }
    };
  }

  template<typename SpT>
  template<typename inverseMatValueType, class ...inverseMatProperties,
           typename inMatValueType,      class ...inMatProperties>
  void
  RealSpaceTools<SpT>::
  inverse(       Kokkos::DynRankView<inverseMatValueType,inverseMatProperties...> inverseMats,
           const Kokkos::DynRankView<inMatValueType,     inMatProperties...>      inMats ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() != inverseMats.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::inverse): Matrix array arguments do not have identical ranks!");
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() < 3 || inMats.rank() > 4, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::inverse): Rank of matrix array must be 3, or 4!");
      for (size_type i=0;i<inMats.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(i) != inverseMats.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::inverse): Dimensions of matrix arguments do not agree!");
      }
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-2) != inMats.extent(inMats.rank()-1), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::inverse): Matrices are not square!");
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-2) < 1 || inMats.extent(inMats.rank()-2) > 3, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::inverse): Spatial dimension must be 1, 2, or 3!");
    }
#endif

    typedef          Kokkos::DynRankView<inverseMatValueType,inverseMatProperties...>      inverseMatViewType;
    typedef          Kokkos::DynRankView<inMatValueType,     inMatProperties...>           inMatViewType;
    typedef          FunctorRealSpaceTools::F_inverse<inverseMatViewType,inMatViewType>    FunctorType;
    typedef typename ExecSpace<typename inMatViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    switch (inMats.rank()) {
    case 3: { // output P,D,D and input P,D,D
      using range_policy_type = Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> >;
      range_policy_type policy(0, inverseMats.extent(0));
      Kokkos::parallel_for( policy, FunctorType(inverseMats, inMats) );
      break;
    }
    case 4: { // output C,P,D,D and input C,P,D,D
      using range_policy_type = Kokkos::Experimental::MDRangePolicy
        < ExecSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
      range_policy_type policy( { 0, 0 },
                                { inverseMats.extent(0), inverseMats.extent(1) } );
      Kokkos::parallel_for( policy, FunctorType(inverseMats, inMats) );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::inverse): Rank of matrix array must be 2, 3, or 4!");
    }
    }
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute determinant see Intrepid2::RealSpaceTools for more
    */ 
    template<typename detArrayViewType,
             typename inMatViewType>
    struct F_det {
      detArrayViewType _detArray;
      inMatViewType    _inMats;

      KOKKOS_INLINE_FUNCTION
      F_det( detArrayViewType detArray_,
             inMatViewType    inMats_ )
        : _detArray(detArray_), _inMats(inMats_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type pt) const {
        const auto mat = Kokkos::subview(_inMats, pt, Kokkos::ALL(), Kokkos::ALL());
        _detArray(pt) = RealSpaceTools<>::Serial::det(mat);
      }

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type cl, 
                      const ordinal_type pt) const {
        const auto mat = Kokkos::subview(_inMats, cl, pt, Kokkos::ALL(), Kokkos::ALL());
        _detArray(cl, pt) = RealSpaceTools<>::Serial::det(mat);
      }
    };
  }

  template<typename SpT>
  template<typename detArrayValueType, class ...detArrayProperties,
           typename inMatValueType,    class ...inMatProperties>
  void
  RealSpaceTools<SpT>::
  det(       Kokkos::DynRankView<detArrayValueType,detArrayProperties...> detArray,
       const Kokkos::DynRankView<inMatValueType,   inMatProperties...>    inMats ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() != detArray.rank()+2, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::det): Determinant and matrix array arguments do not have compatible ranks!");
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() < 3 || inMats.rank() > 4, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::det): Rank of matrix array must be 3 or 4!");
      for (size_type i=0;i<inMats.rank()-2;++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(i) != detArray.extent(i), std::invalid_argument,
                                       ">>> ERROR (RealSpaceTools::det): Dimensions of determinant and matrix array arguments do not agree!");
      }
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-2) != inMats.extent(inMats.rank()-1), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::det): Matrices are not square!");
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-2) < 1 || inMats.extent(inMats.rank()-2) > 3, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::det): Spatial dimension must be 1, 2, or 3!");
    }
#endif

    typedef          Kokkos::DynRankView<detArrayValueType,detArrayProperties...>          detArrayViewType;
    typedef          Kokkos::DynRankView<inMatValueType,   inMatProperties...>             inMatViewType;
    typedef          FunctorRealSpaceTools::F_det<detArrayViewType,inMatViewType>          FunctorType;
    typedef typename ExecSpace<typename inMatViewType::execution_space,SpT>::ExecSpaceType ExecSpaceType;

    switch (detArray.rank()) {
    case 1: { // output P and input P,D,D
      using range_policy_type = Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> >;
      range_policy_type policy(0, detArray.extent(0));
      Kokkos::parallel_for( policy, FunctorType(detArray, inMats) );
      break;
    }
    case 2: { // output C,P and input C,P,D,D
      using range_policy_type = Kokkos::Experimental::MDRangePolicy
        < ExecSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
      range_policy_type policy( { 0, 0 },
                                { detArray.extent(0), detArray.extent(1) } );
      Kokkos::parallel_for( policy, FunctorType(detArray, inMats) );
      break;
    }
    default: {
      INTREPID2_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::det): Rank of detArray must be 1 or 2!");
    }
    }
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to add md arrays see Intrepid2::RealSpaceTools for more
    */ 
    template<typename sumArrayViewType,
             typename inArray1Viewtype,
             typename inArray2ViewType>
    struct F_add {
      sumArrayViewType _sumArray;
      inArray1Viewtype _inArray1;
      inArray2ViewType _inArray2;

      KOKKOS_INLINE_FUNCTION
      F_add( sumArrayViewType sumArray_,
             inArray1Viewtype inArray1_,
             inArray2ViewType inArray2_ )
        : _sumArray(sumArray_), _inArray1(inArray1_), _inArray2(inArray2_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i) const {
        const ordinal_type jend = _sumArray.extent(1);
        const ordinal_type kend = _sumArray.extent(2);
        const ordinal_type lend = _sumArray.extent(3);
        const ordinal_type mend = _sumArray.extent(4);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                _sumArray.access(i,j,k,l,m) = _inArray1.access(i,j,k,l,m) + _inArray2.access(i,j,k,l,m);
      }
    };
  }

  template<typename SpT>
  template<typename sumArrayValueType, class ...sumArrayProperties,
           typename inArray1ValueType, class ...inArray1Properties,
           typename inArray2ValueType, class ...inArray2Properties>
  void
  RealSpaceTools<SpT>::
  add(       Kokkos::DynRankView<sumArrayValueType,sumArrayProperties...> sumArray,
       const Kokkos::DynRankView<inArray1ValueType,inArray1Properties...> inArray1,
       const Kokkos::DynRankView<inArray2ValueType,inArray2Properties...> inArray2 ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inArray1.rank() != inArray2.rank() ||
                                      inArray1.rank() != sumArray.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::add): Array arguments must have identical ranks!");
      for (size_type i=0;i<inArray1.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inArray1.extent(i) != inArray2.extent(i) ||
                                        inArray1.extent(i) != sumArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::add): Dimensions of array arguments do not agree!");
      }
    }
#endif

    typedef          Kokkos::DynRankView<sumArrayValueType,sumArrayProperties...>                     sumArrayViewType;
    typedef          Kokkos::DynRankView<inArray1ValueType,inArray1Properties...>                     inArray1ViewType;
    typedef          Kokkos::DynRankView<inArray2ValueType,inArray2Properties...>                     inArray2ViewType;
    typedef          FunctorRealSpaceTools::F_add<sumArrayViewType,inArray1ViewType,inArray2ViewType> FunctorType;
    typedef typename ExecSpace<typename inArray1ViewType::execution_space,SpT>::ExecSpaceType         ExecSpaceType;

    const auto loopSize = sumArray.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(sumArray, inArray1, inArray2) );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename inoutSumArrayValueType, class ...inoutSumArrayProperties,
           typename inArrayValueType,       class ...inArrayProperties>
  void
  RealSpaceTools<SpT>::
  add(       Kokkos::DynRankView<inoutSumArrayValueType,inoutSumArrayProperties...> inoutSumArray,
       const Kokkos::DynRankView<inArrayValueType,      inArrayProperties...>       inArray ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inArray.rank() != inoutSumArray.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::sum): Array arguments must have identical ranks!");
      for (size_type i=0;i<inArray.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inArray.extent(i) != inoutSumArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::sum): Dimensions of array arguments do not agree!");
      }
    }
#endif
    RealSpaceTools<SpT>::add(inoutSumArray, inoutSumArray, inArray);
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to subtract md arrays see Intrepid2::RealSpaceTools for more
    */ 
    template<typename diffArrayViewType,
             typename inArray1ViewType,
             typename inArray2ViewType>
    struct F_subtract {
            diffArrayViewType _diffArray;
      const inArray1ViewType  _inArray1;
      const inArray2ViewType  _inArray2;

      KOKKOS_INLINE_FUNCTION
      F_subtract( diffArrayViewType diffArray_,
                  inArray1ViewType  inArray1_,
                  inArray2ViewType  inArray2_ )
        : _diffArray(diffArray_), _inArray1(inArray1_), _inArray2(inArray2_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i) const {
        const ordinal_type jend = _diffArray.extent(1);
        const ordinal_type kend = _diffArray.extent(2);
        const ordinal_type lend = _diffArray.extent(3);
        const ordinal_type mend = _diffArray.extent(4);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                _diffArray.access(i,j,k,l,m) = _inArray1.access(i,j,k,l,m) - _inArray2.access(i,j,k,l,m);
      }
    };
  }

  template<typename SpT>
  template<typename diffArrayValueType, class ...diffArrayProperties,
           typename inArray1ValueType,  class ...inArray1Properties,
           typename inArray2ValueType,  class ...inArray2Properties>
  void
  RealSpaceTools<SpT>::
  subtract(       Kokkos::DynRankView<diffArrayValueType,diffArrayProperties...> diffArray,
            const Kokkos::DynRankView<inArray1ValueType, inArray1Properties...>  inArray1,
            const Kokkos::DynRankView<inArray2ValueType, inArray2Properties...>  inArray2 ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inArray1.rank() != inArray2.rank() ||
                                      inArray1.rank() != diffArray.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
      for (size_type i=0;i<inArray1.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inArray1.extent(i) != inArray2.extent(i) ||
                                        inArray1.extent(i) != diffArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
      }
    }
#endif

    typedef          Kokkos::DynRankView<diffArrayValueType,diffArrayProperties...>                         diffArrayViewType;
    typedef          Kokkos::DynRankView<inArray1ValueType, inArray1Properties...>                          inArray1ViewType;
    typedef          Kokkos::DynRankView<inArray2ValueType, inArray2Properties...>                          inArray2ViewType;
    typedef          FunctorRealSpaceTools::F_subtract<diffArrayViewType,inArray1ViewType,inArray2ViewType> FunctorType;
    typedef typename ExecSpace<typename inArray1ViewType::execution_space,SpT>::ExecSpaceType               ExecSpaceType;

    const size_type loopSize = diffArray.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(diffArray, inArray1, inArray2) );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename inoutDiffArrayValueType, class ...inoutDiffArrayProperties,
           typename inArrayValueType,        class ...inArrayProperties>
  void
  RealSpaceTools<SpT>::
  subtract(       Kokkos::DynRankView<inoutDiffArrayValueType,inoutDiffArrayProperties...> inoutDiffArray,
            const Kokkos::DynRankView<inArrayValueType,       inArrayProperties...>        inArray ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inArray.rank() != inoutDiffArray.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::subtract): Array arguments must have identical ranks!");
      for (size_type i=0;i<inArray.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inArray.extent(i) != inoutDiffArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::subtract): Dimensions of array arguments do not agree!");
      }
    }
#endif
    RealSpaceTools<SpT>::subtract(inoutDiffArray, inoutDiffArray, inArray);
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to scale md arrays see Intrepid2::RealSpaceTools for more
    */ 
    template<typename ValueType,
             typename scaledArrayViewType,
             typename inArrayViewType>
    struct F_scale {
            scaledArrayViewType _scaledArray;
      const inArrayViewType     _inArray;
      const ValueType           _alpha;

      KOKKOS_INLINE_FUNCTION
      F_scale( scaledArrayViewType scaledArray_,
               inArrayViewType     inArray_,
               const ValueType     alpha_ )
        : _scaledArray(scaledArray_), _inArray(inArray_), _alpha(alpha_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type i) const {
        const ordinal_type jend = _inArray.extent(1);
        const ordinal_type kend = _inArray.extent(2);
        const ordinal_type lend = _inArray.extent(3);
        const ordinal_type mend = _inArray.extent(4);

        for (ordinal_type j=0;j<jend;++j)
          for (ordinal_type k=0;k<kend;++k)
            for (ordinal_type l=0;l<lend;++l)
              for (ordinal_type m=0;m<mend;++m)
                _scaledArray.access(i,j,k,l,m) = _alpha*_inArray.access(i,j,k,l,m);
      }
    };
  }


  template<typename SpT>
  template<typename ValueType,
           typename scaledArrayValueType, class ...scaledArrayProperties,
           typename inArrayValueType,     class ...inArrayProperties>
  void
  RealSpaceTools<SpT>::
  scale(       Kokkos::DynRankView<scaledArrayValueType,scaledArrayProperties...> scaledArray,
         const Kokkos::DynRankView<inArrayValueType,    inArrayProperties...>     inArray,
         const ValueType alpha ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inArray.rank() > 5, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::scale): Input array container has rank larger than 5.");
      INTREPID2_TEST_FOR_EXCEPTION( inArray.rank() != scaledArray.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::scale): Array arguments must have identical ranks!");
      for (size_type i=0;i<inArray.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inArray.extent(i) != scaledArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::scale): Dimensions of array arguments do not agree!");
      }
    }
#endif

    typedef          Kokkos::DynRankView<scaledArrayValueType,scaledArrayProperties...>               scaledArrayViewtype;
    typedef          Kokkos::DynRankView<inArrayValueType,    inArrayProperties...>                   inArrayViewType;
    typedef          FunctorRealSpaceTools::F_scale<ValueType,scaledArrayViewtype,inArrayViewType> FunctorType;
    typedef typename ExecSpace<typename inArrayViewType::execution_space,SpT>::ExecSpaceType          ExecSpaceType;

    const auto loopSize = scaledArray.extent(0);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(scaledArray, inArray, alpha) );
  }

  // ------------------------------------------------------------------------------------

  template<typename SpT>
  template<typename ValueType,
           typename inoutScaledArrayValueType, class ...inoutScaledArrayProperties>
  void
  RealSpaceTools<SpT>::
  scale(       Kokkos::DynRankView<inoutScaledArrayValueType,inoutScaledArrayProperties...> inoutScaledArray,
         const ValueType alpha ) {
    RealSpaceTools<SpT>::scale(inoutScaledArray, inoutScaledArray, alpha);
  }


  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute dot product see Intrepid2::RealSpaceTools for more
    */ 
    template<typename dotArrayViewType,
             typename inVec1ViewType,
             typename inVec2ViewType>
    struct F_dot {
            dotArrayViewType _dotArray;
      const inVec1ViewType   _inVecs1;
      const inVec2ViewType   _inVecs2;

      KOKKOS_INLINE_FUNCTION
      F_dot( dotArrayViewType dotArray_,
             inVec1ViewType   inVecs1_,
             inVec2ViewType   inVecs2_ )
        : _dotArray(dotArray_), _inVecs1(inVecs1_), _inVecs2(inVecs2_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type iter) const {
        // the rank of normArray is either 1 or 2
        ordinal_type i, j;
        unrollIndex( i, j,
                           _dotArray.extent(0),
                           _dotArray.extent(1),
                           iter );

        const auto r = _inVecs1.rank();
        auto vec1 = ( r == 2 ? Kokkos::subview(_inVecs1, i,    Kokkos::ALL()) :
                               Kokkos::subview(_inVecs1, i, j, Kokkos::ALL()) );
        auto vec2 = ( r == 2 ? Kokkos::subview(_inVecs2, i,    Kokkos::ALL()) :
                               Kokkos::subview(_inVecs2, i, j, Kokkos::ALL()) );

        _dotArray(i,j) = RealSpaceTools<>::Serial::dot(vec1, vec2);
      }
    };
  }

  template<typename SpT>
  template<typename dotArrayValueType, class ...dotArrayProperties,
           typename inVec1ValueType,   class ...inVec1Properties,
           typename inVec2ValueType,   class ...inVec2Properties>
  void
  RealSpaceTools<SpT>::
  dot(       Kokkos::DynRankView<dotArrayValueType,dotArrayProperties...> dotArray,
       const Kokkos::DynRankView<inVec1ValueType,  inVec1Properties...>   inVecs1,
       const Kokkos::DynRankView<inVec2ValueType,  inVec2Properties...>   inVecs2 ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( inVecs1.rank() != (dotArray.rank()+1), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::dot): Ranks of norm and vector array arguments are incompatible!");
      INTREPID2_TEST_FOR_EXCEPTION( inVecs1.rank() != inVecs2.rank(), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::dot): Ranks of input vector arguments must be identical!");
      INTREPID2_TEST_FOR_EXCEPTION( inVecs1.rank() < 2 || inVecs1.rank() > 3, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::dot): Rank of input vector arguments must be 2 or 3!");
      for (size_type i=0;i<inVecs1.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inVecs1.extent(i) != inVecs2.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::dot): Dimensions of input vector arguments do not agree!");
      }
      for (size_type i=0;i<dotArray.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inVecs1.extent(i) != dotArray.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::dot): Dimensions of dot-product and vector arrays do not agree!");
      }
    }
#endif

    typedef          Kokkos::DynRankView<dotArrayValueType,dotArrayProperties...>                 dotArrayViewType;
    typedef          Kokkos::DynRankView<inVec1ValueType,  inVec1Properties...>                   inVec1ViewType;
    typedef          Kokkos::DynRankView<inVec2ValueType,  inVec2Properties...>                   inVec2ViewType;
    typedef          FunctorRealSpaceTools::F_dot<dotArrayViewType,inVec1ViewType,inVec2ViewType> FunctorType;
    typedef typename ExecSpace<typename inVec1ViewType::execution_space,SpT>::ExecSpaceType       ExecSpaceType;

    const auto loopSize = dotArray.extent(0)*dotArray.extent(1);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(dotArray, inVecs1, inVecs2) );
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute matvec see Intrepid2::RealSpaceTools for more
    */ 
    template<typename matVecViewType,
             typename inMatViewType,
             typename inVecViewType>
    struct F_matvec {
            matVecViewType _matVecs;
      const inMatViewType  _inMats;
      const inVecViewType  _inVecs;

      KOKKOS_INLINE_FUNCTION
      F_matvec( matVecViewType matVecs_,
                inMatViewType  inMats_,
                inVecViewType  inVecs_ )
        : _matVecs(matVecs_), _inMats(inMats_), _inVecs(inVecs_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const ordinal_type iter) const {
        const ordinal_type rm = _inMats.rank(), rv = _inVecs.rank(), rr = _matVecs.rank();
        ordinal_type _i = iter, _j = 0;

        if ( rr > 2 )
          unrollIndex( _i, _j,
                             _matVecs.extent(0),
                             _matVecs.extent(1),
                             iter );
        
        auto mat    = ( rm == 2 ? Kokkos::subview(_inMats,        Kokkos::ALL(), Kokkos::ALL()) :
                        rm == 3 ? Kokkos::subview(_inMats,  _i,    Kokkos::ALL(), Kokkos::ALL()) :
                                  Kokkos::subview(_inMats,  _i, _j, Kokkos::ALL(), Kokkos::ALL()) );

        auto vec    = ( rv == 1 ? Kokkos::subview(_inVecs,        Kokkos::ALL()) :
                        rv == 2 ? Kokkos::subview(_inVecs,  _i,    Kokkos::ALL()) :
                                  Kokkos::subview(_inVecs,  _i, _j, Kokkos::ALL()) );

        auto result = ( rr == 1 ? Kokkos::subview(_matVecs,       Kokkos::ALL()) :
                        rr == 2 ? Kokkos::subview(_matVecs, _i,    Kokkos::ALL()) :
                                  Kokkos::subview(_matVecs, _i, _j, Kokkos::ALL()) );

        const ordinal_type iend = result.extent(0);
        const ordinal_type jend = vec.extent(0);

        for (ordinal_type i=0;i<iend;++i) {
          result(i) = 0;
          for (ordinal_type j=0;j<jend;++j)
            result(i) += mat(i, j)*vec(j);
        }
      }
    };
  }

  template<typename SpT>
  template<typename matVecValueType, class ...matVecProperties,
           typename inMatValueType,  class ...inMatProperties,
           typename inVecValueType,  class ...inVecProperties>
  void
  RealSpaceTools<SpT>::
  matvec(       Kokkos::DynRankView<matVecValueType,matVecProperties...> matVecs,
          const Kokkos::DynRankView<inMatValueType, inMatProperties...>  inMats,
          const Kokkos::DynRankView<inVecValueType, inVecProperties...>  inVecs ) {

#ifdef HAVE_INTREPID2_DEBUG
    INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() < 2 || inMats.rank() > 4, std::invalid_argument,
                                  ">>> ERROR (RealSpaceTools::matvec): Rank of matrix array must be 2, 3 or 4!");
    INTREPID2_TEST_FOR_EXCEPTION( matVecs.rank() < 1 || matVecs.rank() > 3, std::invalid_argument,
                                  ">>> ERROR (RealSpaceTools::matvec): Rank of matVecs array must be 1, 2 or 3!");
    INTREPID2_TEST_FOR_EXCEPTION( inVecs.rank() < 1 || inVecs.rank() > 3, std::invalid_argument,
                                  ">>> ERROR (RealSpaceTools::matvec): Rank of inVecs array must be 1, 2 or 3!");
    if (inMats.rank() == 2) {
      // a single matrix, multiple input and output
      INTREPID2_TEST_FOR_EXCEPTION( matVecs.rank() != inVecs.rank(), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::matvec): Vector arrays must be have the same rank!");
      // output must match
      for (ordinal_type i=0;i< (static_cast<ordinal_type>(inVecs.rank())-1);++i) {
        INTREPID2_TEST_FOR_EXCEPTION( matVecs.extent(i) != inVecs.extent(i), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector array arguments do not agree!");
      }
      // matvec compatibility 
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(0) != matVecs.extent(matVecs.rank()-1) || 
                                    inMats.extent(1) != inVecs.extent(inVecs.rank()-1), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::matvec): Matvec dimensions are not compatible each other.");
    } else if (inVecs.rank() == 1) {
      // multiple matrix, single input and multiple output
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() != (matVecs.rank()+1), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::matvec): The result vector and matrix array arguments do not have compatible ranks!");
      for (ordinal_type i=0;i<(static_cast<ordinal_type>(inMats.rank())-2);++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(i) != matVecs.extent(i), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector and matrix array arguments do not agree!");
      }
      // matvec compatibility 
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-2) != matVecs.extent(matVecs.rank()-1) || 
                                    inMats.extent(inMats.rank()-1) != inVecs.extent(0), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::matvec): Matvec dimensions are not compatible each other.");
    } else {
      // multiple matrix, multiple input and multiple output
      INTREPID2_TEST_FOR_EXCEPTION( inMats.rank() != (matVecs.rank()+1), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::matvec): The result vector and matrix array arguments do not have compatible ranks!");
      INTREPID2_TEST_FOR_EXCEPTION( matVecs.rank() != inVecs.rank(), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::matvec): Vector arrays must be have the same rank!");
      for (ordinal_type i=0;i<(static_cast<ordinal_type>(inMats.rank())-2);++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(i) != matVecs.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector and matrix array arguments do not agree!");
      }
      for (size_type i=0;i<inVecs.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( matVecs.extent(i) != inVecs.extent(i), std::invalid_argument,
                                        ">>> ERROR (RealSpaceTools::matvec): Dimensions of vector array arguments do not agree!");
      }
      INTREPID2_TEST_FOR_EXCEPTION( inMats.extent(inMats.rank()-1) != inVecs.extent(inVecs.rank()-1), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::matvec): Matrix column dimension does not match to the length of a vector!");
    }
#endif

    typedef          Kokkos::DynRankView<matVecValueType,matVecProperties...>                    matVecViewType;
    typedef          Kokkos::DynRankView<inMatValueType, inMatProperties...>                     inMatViewType;
    typedef          Kokkos::DynRankView<inVecValueType, inVecProperties...>                     inVecViewType;
    typedef          FunctorRealSpaceTools::F_matvec<matVecViewType,inMatViewType,inVecViewType> FunctorType;
    typedef typename ExecSpace<typename inMatViewType::execution_space,SpT>::ExecSpaceType       ExecSpaceType;

    size_type loopSize = 1;
    const ordinal_type r = matVecs.rank() - 1;
    for (ordinal_type i=0;i<r;++i) 
      loopSize *= matVecs.extent(i);

    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(matVecs, inMats, inVecs) );
  }

  // ------------------------------------------------------------------------------------

  namespace FunctorRealSpaceTools {
    /**
      \brief Functor to compute vecprod see Intrepid2::RealSpaceTools for more
    */ 
    template<typename vecProdViewType,
             typename inLeftViewType,
             typename inRightViewType>
    struct F_vecprod {
            vecProdViewType _vecProd;
      const inLeftViewType  _inLeft;
      const inRightViewType _inRight;
      const bool _is_vecprod_3d;

      KOKKOS_INLINE_FUNCTION
      F_vecprod( vecProdViewType vecProd_,
                 inLeftViewType  inLeft_,
                 inRightViewType inRight_,
                 const bool      is_vecprod_3d_ )
        : _vecProd(vecProd_), _inLeft(inLeft_), _inRight(inRight_), _is_vecprod_3d(is_vecprod_3d_) {}

      KOKKOS_INLINE_FUNCTION
      void operator()(const size_type iter) const {
        ordinal_type i, j;
        unrollIndex( i, j,
                           _inLeft.extent(0),
                           _inLeft.extent(1),
                           iter );

        const auto r = _inLeft.rank();

        auto left   = ( r == 1 ? Kokkos::subview(_inLeft,         Kokkos::ALL()) :
                        r == 2 ? Kokkos::subview(_inLeft,   i,    Kokkos::ALL()) :
                                 Kokkos::subview(_inLeft,   i, j, Kokkos::ALL()) );

        auto right  = ( r == 1 ? Kokkos::subview(_inRight,        Kokkos::ALL()) :
                        r == 2 ? Kokkos::subview(_inRight,  i,    Kokkos::ALL()) :
                                 Kokkos::subview(_inRight,  i, j, Kokkos::ALL()) );

        auto result = ( r == 1 ? Kokkos::subview(_vecProd,        Kokkos::ALL()) :
                        r == 2 ? Kokkos::subview(_vecProd,  i,    Kokkos::ALL()) :
                                 Kokkos::subview(_vecProd,  i, j, Kokkos::ALL()) );

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
  }

  template<typename SpT>
  template<typename vecProdValueType, class ...vecProdProperties,
           typename inLeftValueType,  class ...inLeftProperties,
           typename inRightValueType, class ...inRightProperties>
  void
  RealSpaceTools<SpT>::
  vecprod(       Kokkos::DynRankView<vecProdValueType,vecProdProperties...> vecProd,
           const Kokkos::DynRankView<inLeftValueType, inLeftProperties...>  inLeft,
           const Kokkos::DynRankView<inRightValueType,inRightProperties...> inRight ) {

#ifdef HAVE_INTREPID2_DEBUG
    {
      /*
       *   Check array rank and spatial dimension range (if applicable)
       *      (1) all array arguments are required to have matching dimensions and rank: (D), (I0,D) or (I0,I1,D)
       *      (2) spatial dimension should be 2 or 3
       */
      // (1) check rank range on inLeft and then compare the other arrays with inLeft
      INTREPID2_TEST_FOR_EXCEPTION( inLeft.rank() < 1 || inLeft.rank() > 3, std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::vecprod): Rank of matrix array must be 1, 2, or 3!");
      INTREPID2_TEST_FOR_EXCEPTION( inLeft.rank() != inRight.rank(), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::vecprod): Right and left arrays must be have the same rank!");
      INTREPID2_TEST_FOR_EXCEPTION( inLeft.rank() != vecProd.rank(), std::invalid_argument,
                                    ">>> ERROR (RealSpaceTools::vecprod): Left and vecProd arrays must be have the same rank!");
      for (size_type i=0;i<inLeft.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( inLeft.extent(i) != inRight.extent(i), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::vecprod): Dimensions of matrix arguments do not agree!");
        INTREPID2_TEST_FOR_EXCEPTION( inLeft.extent(i) != vecProd.extent(i), std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::vecprod): Dimensions of matrix arguments do not agree!");
      }

      // (2) spatial dimension ordinal = array rank - 1. Suffices to check only one array because we just
      //     checked whether or not the arrays have matching dimensions.
      INTREPID2_TEST_FOR_EXCEPTION( inLeft.extent(inLeft.rank()-1) < 2 ||
                                      inLeft.extent(inLeft.rank()-1) > 3, std::invalid_argument,
                                      ">>> ERROR (RealSpaceTools::vecprod): Dimensions of arrays (rank-1) must be 2 or 3!");
    }
#endif
    typedef          Kokkos::DynRankView<vecProdValueType, vecProdProperties...>                      vecProdViewType;
    typedef          Kokkos::DynRankView<inLeftValueType, inLeftProperties...>                        inLeftViewType;
    typedef          Kokkos::DynRankView<inRightValueType, inRightProperties...>                      inRightViewType;
    typedef          FunctorRealSpaceTools::F_vecprod<vecProdViewType,inLeftViewType,inRightViewType> FunctorType;
    typedef typename ExecSpace<typename inLeftViewType::execution_space,SpT>::ExecSpaceType           ExecSpaceType;

    const auto r = inLeft.rank();
    const auto loopSize = ( r == 1 ? 1 :
                            r == 2 ? inLeft.extent(0) :
                                     inLeft.extent(0)*inLeft.extent(1) );
    const bool is_vecprod_3d = (inLeft.extent(inLeft.rank() - 1) == 3);
    Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, loopSize);
    Kokkos::parallel_for( policy, FunctorType(vecProd, inLeft, inRight, is_vecprod_3d) );
  }

  // ------------------------------------------------------------------------------------

} // namespace Intrepid2

#endif





// =====================================
// Too much things...
//

//   template<class ...inVec1Properties,
//            class ...inVec2Properties>
//   KOKKOS_INLINE_FUNCTION
//   typename Kokkos::DynRankView<inVec1Properties...>::value_type
//   RealSpaceTools<SpT>::
//   dot( const Kokkos::DynRankView<inVec1Properties...> inVec1,
//        const Kokkos::DynRankView<inVec2Properties...> inVec2 ) {

// #ifdef HAVE_INTREPID2_DEBUG
//     INTREPID2_TEST_FOR_ABORT( inVec1.rank != 1 || inVec2.rank() != 1,
//                               ">>> ERROR (RealSpaceTools::dot): Vector arguments must have rank 1!");
//     INTREPID2_TEST_FOR_ABORT( inVec1.extent(0) != inVec2.extent(0),
//                               ">>> ERROR (RealSpaceTools::dot): Dimensions of vector arguments must agree!");
// #endif
//     typedef typename Kokkos::DynRankView<inVec1Properties...>::value_type value_type;

//     // designed for small problems
//     value_type r_val(0);

//     // * This is probably not necessary
//     if (Kokkos::Impl::is_same<ExecSpace,Kokkos::Serial>::value) {
//       const ordinal_type iend = iVec1.extent(0);
//       for (ordinal_type i=0;i<iend;++i)
//         r_val += inVec1(i)*inVec2(i);
//     } else {
//       struct Functor {
//         Kokkos::DynRankView<inVec1Properties...> _inVec1;
//         Kokkos::DynRankView<inVec2Properties...> _inVec2;

//         KOKKOS_INLINE_FUNCTION
//         Functor(const Kokkos::DynRankView<inVec1Properties...> inVec1_,
//                 const Kokkos::DynRankView<inVec2Properties...> inVec2_);

//         KOKKOS_INLINE_FUNCTION
//         void operator()(ordinal_type i, value_type &dst ) const {
//           dst += _inVec1(i)*_inVec2(i);
//         }

//         KOKKOS_INLINE_FUNCTION
//         void join(volatile value_type &dst ,
//                   const volatile value_type &src) const {
//           dst += src;
//         }
//       };
//       const ordinal_type iend = inVec1.extent(0);
//       Kokkos::RangePolicy<ExecSpaceType,Kokkos::Schedule<Kokkos::Static> > policy(0, iend);
//       Kokkos::parallel_for( policy, Functor(inVec1, inVec2), r_val );
//     }
//     return r_val;
//   }


