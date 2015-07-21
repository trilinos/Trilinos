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


#ifndef PHX_MDFIELD_UTILITIES_DEF_HPP
#define PHX_MDFIELD_UTILITIES_DEF_HPP

namespace PHX {

  template<typename T>
  class MDFieldIterator<T>::Ptr {
    reference_type val_;
  public:
    Ptr (const reference_type& val) : val_(val) {}
    reference_type operator* () const { return val_; }
    // operator-> requires a raw pointer to end its recursive invocation. val_
    // holds an object in memory so that the pointer remains valid for the
    // duration of it->'s use of it.
    const reference_type* operator-> () const { return &val_; }
  };

  template<typename T> inline MDFieldIterator<T>::
  MDFieldIterator (const PHX::MDField<T>& a) : a_(a) {
    rank_ = a.rank();
    i_ = 0;
    done_ = a.size() == 0;
    if (done_) return;
    for (int i = 0; i < rank_; ++i) {
      dimsm1_[i] = a.dimension(i) - 1;
      idxs_[i] = 0;
    }
  }

  template<typename T> inline MDFieldIterator<T>&
  MDFieldIterator<T>::operator++ () {
#if defined(PHX_DEBUG) && ! defined(__CUDA_ARCH__)
    TEUCHOS_TEST_FOR_EXCEPTION(done_, std::logic_error,
                               "operator++ called but done()");
#endif
    // gcc 4.9.1 with -O3 and -Wall produces the warning
    //     array subscript is above array bounds [-Warray-bounds]
    // on the following line:
    //     for (int j = i+1; j < rank_; ++j) idxs_[j] = 0;
    // With -O2 this warning is not issued. I suspect something is happening in
    // the optimizer that causes the static analyzer to conclude that there
    // might be a problem.
    //   I could use
    //     #pragma GCC diagnostic
    // push, pop, and ignored "-Warray-bounds" lines, but I'm worried this will
    // cause unused-pragma warnings in other compilers.
    //   Another way to solve this problem is to provide an extra line of code
    // to the beginning of this method. This gives the static analyzer another
    // separate path to conclude this method is safe.
    //   One might think one could put this line in the ctor, but that does not
    // satisfy the analyzer.
    if (rank_ > max_rank_) throw "never happens";
    for (int i = rank_ - 1; i >= 0; --i)
      if (idxs_[i] < dimsm1_[i]) {
        ++idxs_[i];
        for (int j = i+1; j < rank_; ++j) idxs_[j] = 0;
        ++i_;
        return *this;
      }
    done_ = true;
    return *this;
  }

  template<typename T> inline MDFieldIterator<T>&
  MDFieldIterator<T>::operator++ (int) {
    MDFieldIterator<T> it(*this);
    ++(*this);
    return *this;
  }

  template<typename T> inline bool
  MDFieldIterator<T>::done () const {
    return done_;
  }

  template<typename T> inline typename MDFieldIterator<T>::reference_type
  MDFieldIterator<T>::ref () const {
    switch (rank_) {
    case 1: return a_(idxs_[0]);
    case 2: return a_(idxs_[0], idxs_[1]);
    case 3: return a_(idxs_[0], idxs_[1], idxs_[2]);
    case 4: return a_(idxs_[0], idxs_[1], idxs_[2], idxs_[3]);
    case 5: return a_(idxs_[0], idxs_[1], idxs_[2], idxs_[3], idxs_[4]);
    case 6: return a_(idxs_[0], idxs_[1], idxs_[2], idxs_[3], idxs_[4],
                      idxs_[5]);
    case 7: return a_(idxs_[0], idxs_[1], idxs_[2], idxs_[3], idxs_[4],
                      idxs_[5], idxs_[6]);
    default:
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
                                 "MDField has maximum rank 7");
    }
  }

  template<typename T> inline typename MDFieldIterator<T>::reference_type
  MDFieldIterator<T>::operator* () const {
    return ref();
  }

  template<typename T> inline typename MDFieldIterator<T>::Ptr
  MDFieldIterator<T>::operator-> () const {
    return Ptr(ref());
  }

  template<typename T> inline typename MDFieldIterator<T>::size_type
  MDFieldIterator<T>::idx () const {
    return i_;
  }

} // namespace PHX

#endif // PHX_MDFIELD_UTILITIES_DEF_HPP
