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
    template<typename outputViewType,
             typename inputViewType,
             ordinal_type valRank>
    struct F_clone {
            outputViewType _output;
      const inputViewType _input;

      KOKKOS_INLINE_FUNCTION
      F_clone(outputViewType output_,
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
  
  template<typename SpT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void ArrayTools<SpT>::
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

    typedef Kokkos::DynRankView<outputValueType,outputProperties...> outputViewType;
    typedef Kokkos::DynRankView<inputValueType, inputProperties...>  inputViewType; 
    typedef typename ExecSpace< typename inputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;
    
    using range_policy_type = Kokkos::Experimental::MDRangePolicy
      < ExecSpaceType, Kokkos::Experimental::Rank<3>, Kokkos::IndexType<ordinal_type> >;
    
    const ordinal_type
      C = output.extent(0),
      F = output.extent(1),
      P = output.extent(2);
    
    range_policy_type policy( { 0, 0, 0 },
                              { C, F, P } );
    const ordinal_type valRank = output.rank() - 3;
    switch (valRank) {
    case 0: {
      typedef FunctorArrayTools::F_clone<outputViewType,inputViewType,0> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 1: {
      typedef FunctorArrayTools::F_clone<outputViewType,inputViewType,1> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 2: {
      typedef FunctorArrayTools::F_clone<outputViewType,inputViewType,2> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    }
  }

  template<typename SpT>
  template<typename outputValueType, class ...outputProperties,
           typename inputValueType,  class ...inputProperties>
  void ArrayTools<SpT>::
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

    typedef Kokkos::DynRankView<outputValueType,outputProperties...> outputViewType;
    typedef Kokkos::DynRankView<inputValueType, inputProperties...>  inputViewType; 
    typedef typename ExecSpace< typename inputViewType::execution_space , SpT >::ExecSpaceType ExecSpaceType;
    
    using range_policy_type = Kokkos::Experimental::MDRangePolicy
      < ExecSpaceType, Kokkos::Experimental::Rank<2>, Kokkos::IndexType<ordinal_type> >;
    
    const ordinal_type
      C = output.extent(0),
      P = output.extent(1);
    
    range_policy_type policy( { 0, 0 },
                              { C, P } );
    const ordinal_type valRank = output.rank() - 2;
    switch (valRank) {
    case 0: {
      typedef FunctorArrayTools::F_clone<outputViewType,inputViewType,0> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 1: {
      typedef FunctorArrayTools::F_clone<outputViewType,inputViewType,1> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    case 2: {
      typedef FunctorArrayTools::F_clone<outputViewType,inputViewType,2> FunctorType;
      Kokkos::parallel_for( policy, FunctorType(output, input) );
      break;
    }
    }
  }
  
} // end namespace Intrepid2

#endif
