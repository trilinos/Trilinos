// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
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
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#ifndef PHX_MDFIELD_UTILITIES_HPP
#define PHX_MDFIELD_UTILITIES_HPP

#include "Phalanx_MDField.hpp"

namespace PHX {

  /*! \brief Iterate over entries of an MDField<DataT>, a runtime MDField.
   *
   * This iterator is one way to replace usage of the deprecated
   * operator[]. Regardless of underlying implementation details, iteration
   * proceeds so that in A(i,j,k), k is the fastest index and i is the slowest
   * one.
   *
   * Example usage:
   * \code
   *    int i = 0;
   *    for (PHAL::MDFieldIterator<ScalarT> d(array); ! d.done() ; ++d, ++i)
   *      *d = val[i];
   * \endcode
   */
  template<typename T>
  class MDFieldIterator {
  public:
    //! Reference to an entry of the MDField.
    typedef typename
      PHX::MDFieldTypeTraits<typename MDField<T>::array_type>::return_type
      reference_type;
    //! Index type.
    typedef typename PHX::MDField<T>::size_type size_type;
    //! User constructor.
    explicit MDFieldIterator(const PHX::MDField<T>& a);
    //! Increment the iterator. Efficient.
    MDFieldIterator<T>& operator++();
    //! Like all postfix ++, this one is inefficient. Convenience only.
    MDFieldIterator<T>& operator++(int);
    //! Returns whether the iterator has reached the end.
    bool done() const;
    //! Get a reference to the current value.
    reference_type ref() const;
    //! Syntactic wrapper to \c ref().
    reference_type operator*() const;
    //! Pointer type for \c operator->.
    class Ptr;
    //! Syntactic wrapper to \c (*it). .
    Ptr operator->() const;
    //! Get the index of the current value.
    size_type idx() const;
  
  private:
    // All stack-allocated variables.
    // Runtime MDField to iterate over.
    const PHX::MDField<T>& a_;
    // For index arithmetic.
    static const int max_rank_ = 7;
    size_type dimsm1_[max_rank_], idxs_[max_rank_];
    // i'th entry in a_.
    size_type i_;
    // a_'s rank.
    int rank_;
    // Whether the iterator has reached the end of the array.
    bool done_;
  };

} // namespace PHX

#include "Phalanx_MDField_Utilities_Def.hpp"

#endif // PHX_MDFIELD_UTILITIES_HPP
