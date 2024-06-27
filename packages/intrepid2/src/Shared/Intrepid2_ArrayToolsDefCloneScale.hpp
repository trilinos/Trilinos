// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_ArrayToolsDefCloneScale.hpp
    \brief  Definition file for clone / scale operations of the array tools interface.
    \author Created by P. Bochev and D. Ridzal.
            Kokkorized by Kyungjoo Kim
*/

#ifndef __INTREPID2_ARRAYTOOLS_DEF_CLONESCALE_HPP__
#define __INTREPID2_ARRAYTOOLS_DEF_CLONESCALE_HPP__

namespace Intrepid2 {
  
  namespace FunctorArrayTools {

    /**
      \brief Functor for clone see Intrepid2::ArrayTools for more
    */ 
    template<typename OutputViewType,
             typename inputViewType,
             ordinal_type valRank>
    struct F_clone {
            OutputViewType _output;
      const inputViewType _input;

      KOKKOS_INLINE_FUNCTION
      F_clone(OutputViewType output_,
              inputViewType input_)
        : _output(output_),
          _input(input_) {}
      
      // clone data
      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type cl,
                 const ordinal_type pt) const {
        switch (valRank) {
        case 0: {
          _output.access(cl, pt) = _input.access(pt);
          break;
        }
        case 1: {
          const ordinal_type iend = _output.extent(2);
          for (ordinal_type i=0;i<iend;++i)
            _output.access(cl, pt, i) = _input.access(pt, i);
          break;
        }
        case 2: {
          const ordinal_type
            iend = _output.extent(2),
            jend = _output.extent(3);
          for (ordinal_type i=0;i<iend;++i)
            for (ordinal_type j=0;j<jend;++j)
              _output.access(cl, pt, i, j) = _input.access(pt, i, j);
          break;
        }
        }
      }

      // clone fields
      KOKKOS_INLINE_FUNCTION
      void
      operator()(const ordinal_type cl,
                 const ordinal_type bf,
                 const ordinal_type pt) const {
        switch (valRank) {
        case 0: {
          _output.access(cl, bf, pt) = _input.access(bf, pt);
          break;
        }
        case 1: {
          const ordinal_type iend = _output.extent(3);
          for (ordinal_type i=0;i<iend;++i)
            _output.access(cl, bf, pt, i) = _input.access(bf, pt, i);
          break;
        }
        case 2: {
          const ordinal_type
            iend = _output.extent(3),
            jend = _output.extent(4);
          for (ordinal_type i=0;i<iend;++i)
            for (ordinal_type j=0;j<jend;++j)
              _output.access(cl, bf, pt, i, j) = _input.access(bf, pt, i, j);
          break;
        }
        }
      }
    };
  } 
  
  template<typename DeviceType>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void ArrayTools<DeviceType>::
  cloneFields(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
               const Kokkos::DynRankView<inputValueType, inputProperties...>  input ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( ( input.rank() < 2 || input.rank() > 4 ), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::clone): Input fields container must have rank 2, 3, or 4.");
      INTREPID2_TEST_FOR_EXCEPTION( ( output.rank() != (input.rank()+1) ), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::clone): The rank of the input fields container must be one less than the rank of the output fields container.");
      for (size_type i=0;i< input.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( (input.extent(i) != output.extent(i+1)), std::invalid_argument,
                                        ">>> ERROR (ArrayTools::clone): Dimensions of input and output fields containers do not match.");
      }
    }
#endif

    typedef Kokkos::DynRankView<outputValueType,outputProperties...> OutputViewType;
    typedef Kokkos::DynRankView<inputValueType, inputProperties...>  inputViewType; 

    using range_policy_type = Kokkos::MDRangePolicy
      < ExecSpaceType, Kokkos::Rank<3>, Kokkos::IndexType<ordinal_type> >;

    range_policy_type policy( { 0, 0, 0 },
                              { /*C*/ output.extent(0), /*F*/ output.extent(1), /*P*/ output.extent(2) } );
    const ordinal_type valRank = output.rank() - 3;
    switch (valRank) {
    case 0: {
      typedef FunctorArrayTools::F_clone<OutputViewType,inputViewType,0> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 1: {
      typedef FunctorArrayTools::F_clone<OutputViewType,inputViewType,1> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 2: {
      typedef FunctorArrayTools::F_clone<OutputViewType,inputViewType,2> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    }
  }

  template<typename DeviceType>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void ArrayTools<DeviceType>::
  cloneData(       Kokkos::DynRankView<outputValueType,outputProperties...> output,
             const Kokkos::DynRankView<inputValueType, inputProperties...>  input ) {
#ifdef HAVE_INTREPID2_DEBUG
    {
      INTREPID2_TEST_FOR_EXCEPTION( ( input.rank() < 1 || input.rank() > 3 ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::clone): Input fields container must have rank 1, 2, or 3.");
      INTREPID2_TEST_FOR_EXCEPTION( ( output.rank() != (input.rank()+1) ), std::invalid_argument,
                                    ">>> ERROR (ArrayTools::clone): The rank of the input fields container must be one less than the rank of the output fields container.");
      for (ordinal_type i=0;i<input.rank();++i) {
        INTREPID2_TEST_FOR_EXCEPTION( (input.extent(i) != output.extent(i+1)), std::invalid_argument,
                                      ">>> ERROR (ArrayTools::clone): Dimensions of input and output fields containers do not match.");
      }
    }
#endif

    typedef Kokkos::DynRankView<outputValueType,outputProperties...> OutputViewType;
    typedef Kokkos::DynRankView<inputValueType, inputProperties...>  inputViewType; 
    
    using range_policy_type = Kokkos::MDRangePolicy
      < ExecSpaceType, Kokkos::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    
    range_policy_type policy( { 0, 0 },
                              { /*C*/ output.extent(0), /*P*/ output.extent(1) } );
    const ordinal_type valRank = output.rank() - 2;
    switch (valRank) {
    case 0: {
      typedef FunctorArrayTools::F_clone<OutputViewType,inputViewType,0> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 1: {
      typedef FunctorArrayTools::F_clone<OutputViewType,inputViewType,1> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 2: {
      typedef FunctorArrayTools::F_clone<OutputViewType,inputViewType,2> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    }
  }
  
} // end namespace Intrepid2

#endif
