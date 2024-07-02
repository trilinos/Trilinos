// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*! \file AnasaziBasicSort.hpp
  \brief Basic implementation of the Anasazi::SortManager class
*/

#ifndef ANASAZI_BASIC_SORT_HPP
#define ANASAZI_BASIC_SORT_HPP

/*!    \class Anasazi::BasicSort
       \brief An implementation of the Anasazi::SortManager that performs a collection
       of common sorting techniques.

       \author Chris Baker, Ulrich Hetmaniuk, Rich Lehoucq, and Heidi Thornquist
*/

#include "AnasaziConfigDefs.hpp"
#include "AnasaziSortManager.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Anasazi {

  template<class MagnitudeType>
  class BasicSort : public SortManager<MagnitudeType> {

  public:

    /*! \brief Parameter list driven constructor

        This constructor accepts a paramter list with the following options:
       - \c "Sort Strategy" - a \c string specifying the desired sorting strategy. See setSortType() for valid options.
    */
    BasicSort( Teuchos::ParameterList &pl );

    /*! \brief String driven constructor

        Directly pass the string specifying sort strategy. See setSortType() for valid options.
    */
    BasicSort( const std::string &which = "LM" );

    //! Destructor
    virtual ~BasicSort();

    //! Set sort type
    /**
       @param which [in] The eigenvalues of interest for this eigenproblem.
       - \c "LM" - Largest Magnitude [ default ]
       - \c "SM" - Smallest Magnitude
       - \c "LR" - Largest Real
       - \c "SR" - Smallest Real
       - \c "LI" - Largest Imaginary
       - \c "SI" - Smallest Imaginary
    */
    void setSortType( const std::string &which );

    /*! \brief Sort real eigenvalues, optionally returning the permutation vector.

       \note This method is not valid when the sort manager is configured for "LI" or "SI" sorting
       (i.e., sorting by the imaginary components). Calling this method in that scenario will result
       in a SortManagerError exception.

       @param evals [in/out] Vector of length at least \c n containing the eigenvalues to be sorted.  <br>
                     On output, the first \c n eigenvalues will be sorted. The rest will be unchanged.

       @param perm [out] Vector of length at least \c n to store the permutation index (optional).  <br>
       If specified, on output the first \c n eigenvalues will contain the permutation indices, in the range <tt>[0,n-1]</tt>, such that <tt>evals_out[i] = evals_in[perm[i]]</tt>

       @param n [in] Number of values in evals to be sorted. If <tt>n == -1</tt>, all values will be sorted.
    */
    void sort(std::vector<MagnitudeType> &evals, Teuchos::RCP<std::vector<int> > perm = Teuchos::null, int n = -1) const;

    /*! \brief Sort complex eigenvalues, optionally returning the permutation vector.

       This routine takes two vectors, one for each part of a complex
       eigenvalue. This is helpful for solving real, non-symmetric eigenvalue
       problems.

       @param r_evals [in/out] Vector of length at least \c n containing the real part of the eigenvalues to be sorted.  <br>
                     On output, the first \c n eigenvalues will be sorted. The rest will be unchanged.

       @param i_evals [in/out] Vector of length at least \c n containing the imaginary part of the eigenvalues to be sorted.  <br>
                     On output, the first \c n eigenvalues will be sorted. The rest will be unchanged.

       @param perm [out] Vector of length at least \c n to store the permutation index (optional).  <br>
       If specified, on output the first \c n eigenvalues will contain the permutation indices, in the range <tt>[0,n-1]</tt>, such that <tt>r_evals_out[i] = r_evals_in[perm[i]]</tt>
       and similarly for \c i_evals.

       @param n [in] Number of values in \c r_evals, \c i_evals to be sorted. If <tt>n == -1</tt>, all values will be sorted, as decided by the minimum of the length of \c r_evals and the length of \c i_evals.
    */
    void sort(std::vector<MagnitudeType> &r_evals,
              std::vector<MagnitudeType> &i_evals,
              Teuchos::RCP<std::vector<int> > perm = Teuchos::null,
              int n = -1) const;

  protected:

    // enum for sort type
    enum SType {
      LM, SM,
      LR, SR,
      LI, SI
    };
    SType which_;

    // sorting methods
    template <class LTorGT>
    struct compMag {
      // for real-only LM,SM
      bool operator()(MagnitudeType, MagnitudeType);
      // for real-only LM,SM with permutation
      template <class First, class Second>
        bool operator()(std::pair<First,Second>, std::pair<First,Second>);
    };

    template <class LTorGT>
    struct compMag2 {
      // for real-imag LM,SM
      bool operator()(std::pair<MagnitudeType,MagnitudeType>, std::pair<MagnitudeType,MagnitudeType>);
      // for real-imag LM,SM with permutation
      template <class First, class Second>
        bool operator()(std::pair<First,Second>, std::pair<First,Second>);
    };

    template <class LTorGT>
    struct compAlg {
      // for real-imag LR,SR,LI,SI
      bool operator()(MagnitudeType, MagnitudeType);
      template <class First, class Second>
        bool operator()(std::pair<First,Second>, std::pair<First,Second>);
    };

    template <typename pair_type>
    struct sel1st
    {
      const typename pair_type::first_type &operator()(const pair_type &v) const;
    };

    template <typename pair_type>
    struct sel2nd
    {
      const typename pair_type::second_type &operator()(const pair_type &v) const;
    };
  };


  ////////////////////////////////////////////////////////////////////////
  //  IMPLEMENTATION
  ////////////////////////////////////////////////////////////////////////

  template<class MagnitudeType>
  BasicSort<MagnitudeType>::BasicSort(Teuchos::ParameterList &pl)
  {
    std::string which = "LM";
    which = pl.get("Sort Strategy",which);
    setSortType(which);
  }

  template<class MagnitudeType>
  BasicSort<MagnitudeType>::BasicSort(const std::string &which)
  {
    setSortType(which);
  }

  template<class MagnitudeType>
  BasicSort<MagnitudeType>::~BasicSort()
  {}

  template<class MagnitudeType>
  void BasicSort<MagnitudeType>::setSortType(const std::string &which)
  {
    // make upper case
    std::string whichlc(which);
    std::transform(which.begin(),which.end(),whichlc.begin(),(int(*)(int)) std::toupper);
    if (whichlc == "LM") {
      which_ = LM;
    }
    else if (whichlc == "SM") {
      which_ = SM;
    }
    else if (whichlc == "LR") {
      which_ = LR;
    }
    else if (whichlc == "SR") {
      which_ = SR;
    }
    else if (whichlc == "LI") {
      which_ = LI;
    }
    else if (whichlc == "SI") {
      which_ = SI;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Anasazi::BasicSort::setSortType(): sorting order is not valid");
    }
  }

  template<class MagnitudeType>
  void BasicSort<MagnitudeType>::sort(std::vector<MagnitudeType> &evals, Teuchos::RCP<std::vector<int> > perm, int n) const
  {
    TEUCHOS_TEST_FOR_EXCEPTION(n < -1, std::invalid_argument, "Anasazi::BasicSort::sort(r): n must be n >= 0 or n == -1.");
    if (n == -1) {
      n = evals.size();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(evals.size() < (unsigned int) n,
                       std::invalid_argument, "Anasazi::BasicSort::sort(r): eigenvalue vector size isn't consistent with n.");
    if (perm != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(perm->size() < (unsigned int) n,
                         std::invalid_argument, "Anasazi::BasicSort::sort(r): permutation vector size isn't consistent with n.");
    }

    typedef std::greater<MagnitudeType> greater_mt;
    typedef std::less<MagnitudeType>    less_mt;

    if (perm == Teuchos::null) {
      //
      // if permutation index is not required, just sort using the values
      //
      if (which_ == LM ) {
        std::sort(evals.begin(),evals.begin()+n,compMag<greater_mt>());
      }
      else if (which_ == SM) {
        std::sort(evals.begin(),evals.begin()+n,compMag<less_mt>());
      }
      else if (which_ == LR) {
        std::sort(evals.begin(),evals.begin()+n,compAlg<greater_mt>());
      }
      else if (which_ == SR) {
        std::sort(evals.begin(),evals.begin()+n,compAlg<less_mt>());
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, SortManagerError, "Anasazi::BasicSort::sort(r): LI or SI sorting invalid for real scalar types." );
      }
    }
    else {
      //
      // if permutation index is required, we must sort the two at once
      // in this case, we arrange a pair structure: <value,index>
      // default comparison operator for pair<t1,t2> is lexographic:
      //    compare first t1, then t2
      // this works fine for us here.
      //

      // copy the values and indices into the pair structure
      std::vector< std::pair<MagnitudeType,int> > pairs(n);
      for (int i=0; i<n; i++) {
        pairs[i] = std::make_pair(evals[i],i);
      }

      // sort the pair structure
      if (which_ == LM) {
        std::sort(pairs.begin(),pairs.begin()+n,compMag<greater_mt>());
      }
      else if (which_ == SM) {
        std::sort(pairs.begin(),pairs.begin()+n,compMag<less_mt>());
      }
      else if (which_ == LR) {
        std::sort(pairs.begin(),pairs.begin()+n,compAlg<greater_mt>());
      }
      else if (which_ == SR) {
        std::sort(pairs.begin(),pairs.begin()+n,compAlg<less_mt>());
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true, SortManagerError, "Anasazi::BasicSort::sort(r): LI or SI sorting invalid for real scalar types." );
      }

      // copy the values and indices out of the pair structure
      std::transform(pairs.begin(),pairs.end(),evals.begin(),sel1st< std::pair<MagnitudeType,int> >());
      std::transform(pairs.begin(),pairs.end(),perm->begin(),sel2nd< std::pair<MagnitudeType,int> >());
    }
  }


  template<class T1, class T2>
  class MakePairOp {
  public:
    std::pair<T1,T2> operator()( const T1 &t1, const T2 &t2 )
      { return std::make_pair(t1, t2); }
  };


  template<class MagnitudeType>
  void BasicSort<MagnitudeType>::sort(std::vector<MagnitudeType> &r_evals,
                                      std::vector<MagnitudeType> &i_evals,
                                      Teuchos::RCP< std::vector<int> > perm,
                                      int n) const
  {

    //typedef typename std::vector<MagnitudeType>::iterator r_eval_iter_t; // unused
    //typedef typename std::vector<MagnitudeType>::iterator i_eval_iter_t; // unused

    TEUCHOS_TEST_FOR_EXCEPTION(n < -1, std::invalid_argument, "Anasazi::BasicSort::sort(r,i): n must be n >= 0 or n == -1.");
    if (n == -1) {
      n = r_evals.size() < i_evals.size() ? r_evals.size() : i_evals.size();
    }
    TEUCHOS_TEST_FOR_EXCEPTION(r_evals.size() < (unsigned int) n || i_evals.size() < (unsigned int) n,
                       std::invalid_argument, "Anasazi::BasicSort::sort(r,i): eigenvalue vector size isn't consistent with n.");
    if (perm != Teuchos::null) {
      TEUCHOS_TEST_FOR_EXCEPTION(perm->size() < (unsigned int) n,
                         std::invalid_argument, "Anasazi::BasicSort::sort(r,i): permutation vector size isn't consistent with n.");
    }

    typedef std::greater<MagnitudeType> greater_mt;
    typedef std::less<MagnitudeType>    less_mt;

    //
    // put values into pairs
    //
    if (perm == Teuchos::null) {
      //
      // not permuting, so we don't need indices in the pairs
      //
      std::vector< std::pair<MagnitudeType,MagnitudeType> > pairs(n);
      // for LM,SM, the order doesn't matter
      // for LI,SI, the imaginary goes first
      // for LR,SR, the real goes in first
      if (which_ == LR || which_ == SR || which_ == LM || which_ == SM) {
        std::transform(
          r_evals.begin(), r_evals.begin()+n,
          i_evals.begin(), pairs.begin(),
          MakePairOp<MagnitudeType,MagnitudeType>());
      }
      else {
        std::transform(
          i_evals.begin(), i_evals.begin()+n,
          r_evals.begin(), pairs.begin(),
          MakePairOp<MagnitudeType,MagnitudeType>());
      }

      if (which_ == LR || which_ == LI) {
        std::sort(pairs.begin(),pairs.end(),compAlg<greater_mt>());
      }
      else if (which_ == SR || which_ == SI) {
        std::sort(pairs.begin(),pairs.end(),compAlg<less_mt>());
      }
      else if (which_ == LM) {
        std::sort(pairs.begin(),pairs.end(),compMag2<greater_mt>());
      }
      else {
        std::sort(pairs.begin(),pairs.end(),compMag2<less_mt>());
      }

      // extract the values
      // for LM,SM,LR,SR: order is (real,imag)
      // for LI,SI: order is (imag,real)
      if (which_ == LR || which_ == SR || which_ == LM || which_ == SM) {
        std::transform(pairs.begin(),pairs.end(),r_evals.begin(),sel1st< std::pair<MagnitudeType,MagnitudeType> >());
        std::transform(pairs.begin(),pairs.end(),i_evals.begin(),sel2nd< std::pair<MagnitudeType,MagnitudeType> >());
      }
      else {
        std::transform(pairs.begin(),pairs.end(),r_evals.begin(),sel2nd< std::pair<MagnitudeType,MagnitudeType> >());
        std::transform(pairs.begin(),pairs.end(),i_evals.begin(),sel1st< std::pair<MagnitudeType,MagnitudeType> >());
      }
    }
    else {
      //
      // permuting, we need indices in the pairs
      //
      std::vector< std::pair< std::pair<MagnitudeType,MagnitudeType>, int > > pairs(n);
      // for LM,SM, the order doesn't matter
      // for LI,SI, the imaginary goes first
      // for LR,SR, the real goes in first
      if (which_ == LR || which_ == SR || which_ == LM || which_ == SM) {
        for (int i=0; i<n; i++) {
          pairs[i] = std::make_pair(std::make_pair(r_evals[i],i_evals[i]),i);
        }
      }
      else {
        for (int i=0; i<n; i++) {
          pairs[i] = std::make_pair(std::make_pair(i_evals[i],r_evals[i]),i);
        }
      }

      if (which_ == LR || which_ == LI) {
        std::sort(pairs.begin(),pairs.end(),compAlg<greater_mt>());
      }
      else if (which_ == SR || which_ == SI) {
        std::sort(pairs.begin(),pairs.end(),compAlg<less_mt>());
      }
      else if (which_ == LM) {
        std::sort(pairs.begin(),pairs.end(),compMag2<greater_mt>());
      }
      else {
        std::sort(pairs.begin(),pairs.end(),compMag2<less_mt>());
      }

      // extract the values
      // for LM,SM,LR,SR: order is (real,imag)
      // for LI,SI: order is (imag,real)
      if (which_ == LR || which_ == SR || which_ == LM || which_ == SM) {
        for (int i=0; i<n; i++) {
          r_evals[i] = pairs[i].first.first;
          i_evals[i] = pairs[i].first.second;
          (*perm)[i] = pairs[i].second;
        }
      }
      else {
        for (int i=0; i<n; i++) {
          i_evals[i] = pairs[i].first.first;
          r_evals[i] = pairs[i].first.second;
          (*perm)[i] = pairs[i].second;
        }
      }
    }
  }


  template<class MagnitudeType>
  template<class LTorGT>
  bool BasicSort<MagnitudeType>::compMag<LTorGT>::operator()(MagnitudeType v1, MagnitudeType v2)
  {
    typedef Teuchos::ScalarTraits<MagnitudeType> MTT;
    LTorGT comp;
    return comp( MTT::magnitude(v1), MTT::magnitude(v2) );
  }

  template<class MagnitudeType>
  template<class LTorGT>
  bool BasicSort<MagnitudeType>::compMag2<LTorGT>::operator()(std::pair<MagnitudeType,MagnitudeType> v1, std::pair<MagnitudeType,MagnitudeType> v2)
  {
    MagnitudeType m1 = v1.first*v1.first + v1.second*v1.second;
    MagnitudeType m2 = v2.first*v2.first + v2.second*v2.second;
    LTorGT comp;
    return comp( m1, m2 );
  }

  template<class MagnitudeType>
  template<class LTorGT>
  bool BasicSort<MagnitudeType>::compAlg<LTorGT>::operator()(MagnitudeType v1, MagnitudeType v2)
  {
    LTorGT comp;
    return comp( v1, v2 );
  }

  template<class MagnitudeType>
  template<class LTorGT>
  template<class First, class Second>
  bool BasicSort<MagnitudeType>::compMag<LTorGT>::operator()(std::pair<First,Second> v1, std::pair<First,Second> v2) {
    return (*this)(v1.first,v2.first);
  }

  template<class MagnitudeType>
  template<class LTorGT>
  template<class First, class Second>
  bool BasicSort<MagnitudeType>::compMag2<LTorGT>::operator()(std::pair<First,Second> v1, std::pair<First,Second> v2) {
    return (*this)(v1.first,v2.first);
  }

  template<class MagnitudeType>
  template<class LTorGT>
  template<class First, class Second>
  bool BasicSort<MagnitudeType>::compAlg<LTorGT>::operator()(std::pair<First,Second> v1, std::pair<First,Second> v2) {
    return (*this)(v1.first,v2.first);
  }

  template <class MagnitudeType>
  template <typename pair_type>
  const typename pair_type::first_type &
  BasicSort<MagnitudeType>::sel1st<pair_type>::operator()(const pair_type &v) const
  {
    return v.first;
  }

  template <class MagnitudeType>
  template <typename pair_type>
  const typename pair_type::second_type &
  BasicSort<MagnitudeType>::sel2nd<pair_type>::operator()(const pair_type &v) const
  {
    return v.second;
  }

} // namespace Anasazi

#endif // ANASAZI_BASIC_SORT_HPP

