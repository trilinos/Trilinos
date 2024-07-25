// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RBGEN_FILTER_H
#define RBGEN_FILTER_H

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_RCP.hpp"

// finish: amend filter() to take a pointer to the IncSVDPOD object

namespace RBGen {

  //! \brief Enumerated type to list the selection criteria for singular value filters.
  enum SortType {
    LARGEST,   /*!< Select the largest singular values */
    SMALLEST   /*!< Select the smallest singular values */
  };

  //! Class for selecting desired singular values
  /*!
   *
   */
  template <class ScalarType>
  class Filter {
    public: 

    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    Filter(SortType which = LARGEST) {};

    //! Destructor.
    virtual ~Filter() {};
    //@}

    //! @name Filter Method
    //@{

    //! Selects targetted singular values.

    /*!
     * This method accepts a list of singular values (in decreasing order, as returned by GESVD) and 
     * returns a list of indices consisting of the values
     * allowed through the filter. 
     *
     * @param svals  [in] Singular values under consideration by the filter.
     * \returns A vector of ints, corresponding to the indices of the targetted singular values.
     * The indices are in ascending order and employ zero-based indexing.
     */
    virtual std::vector<int> filter(const std::vector<ScalarType> &svals) = 0;
    //@}
  };

  //! Range-based filter
  template <class ScalarType>
  class RangeFilter : public Filter<ScalarType> {
    public: 
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    RangeFilter(SortType which, int minRank = 1, int maxRank = 1) 
      : which_(which) {
      TEUCHOS_TEST_FOR_EXCEPTION(minRank < 1,std::invalid_argument,"RangeFilter: minRank must be > 1");
      TEUCHOS_TEST_FOR_EXCEPTION(maxRank < 1,std::invalid_argument,"RangeFilter: maxRank must be > 1");
      minRank_ = minRank;
      maxRank_ = maxRank;
    };

    //! Destructor.
    virtual ~RangeFilter() {};
    //@}

    //! @name Filter Method
    //@{

    //! Pass at most maxRank and at most minRank values through the filter.
    /*!
     * \note It is assumed that svals are sorted in decreasing order.
     */
    std::vector<int> filter(const std::vector<ScalarType> &svals) {
      int num = (unsigned int)svals.size();
      if      (num > maxRank_) num = maxRank_;
      else if (minRank_ < num) num = minRank_;
      std::vector<int> ret;
      ret.reserve(num);
      if (LARGEST == which_) {
        for (int i=0; i<num; ++i) {
          ret.push_back(i);
        }
      }
      else if (SMALLEST == which_) {
        for (unsigned int i=svals.size()-num; i<svals.size(); ++i) {
          ret.push_back(i);
        }
      }
      return ret;
    }
    //@}

    private: 
      SortType which_;
      int minRank_;
      int maxRank_;
  };

  //! Threshold-based filter
  template <class ScalarType>
  class ThreshFilter : public Filter<ScalarType> {
    public: 
    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    /*!
     *   The threshold definition depends on \c which:
     *     - ::LARGEST: value \c cand passes the filter if \f$cand \geq tval\f$
     *     - ::SMALLEST: value \c cand passes the filter if \f$cand \leq tval\f$
     *   where \c tval depends on \c absthresh:
     *     - \c false: 
     *        - \f$tval \doteq thresh \cdot max(svals)\f$
     *        - \f$tval \doteq thresh \cdot min(svals)\f$
     *     - \c true: \f$tval \doteq thresh\f$
     */
    ThreshFilter(SortType which, bool absthresh, 
                 typename Teuchos::ScalarTraits<ScalarType>::magnitudeType thresh)
      : which_(which), absthresh_(absthresh) {
      TEUCHOS_TEST_FOR_EXCEPTION(thresh < Teuchos::ScalarTraits<ScalarType>::zero(),
                         std::invalid_argument,"ThreshFilter: minRank must be > 1");
      thresh_ = thresh;
    };

    //! Destructor.
    virtual ~ThreshFilter() {};
    //@}

    //! @name Filter Method
    //@{

    //! Return an index for the singular values falling within the threshold.

    /*!
     *  \note It is assumed that svals are sorted in decreasing order.
     */
     std::vector<int> filter(const std::vector<ScalarType> &svals) {
      const int last = (unsigned int)svals.size() - 1;

      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType tval;
      if (absthresh_) {
        tval = thresh_;
      }
      else {
        if (LARGEST == which_) {
          tval = thresh_ * svals[0];
        }
        else if (SMALLEST == which_) {
          tval = thresh_ * svals[last];
        }
      }

      std::vector<int> ret;
      if (LARGEST == which_) {
        const int num = find(svals.begin(),svals.end(),
                             [=] (const ScalarType sval) {
                               return sval < tval;
                             }) - svals.begin();
        ret.resize(num);
        for (int i=0; i<num; ++i) {
          ret.push_back(i);
        }
      }
      else if (SMALLEST == which_) {
        const int num = svals.end() - find(svals.begin(),svals.end(),
                                           [=] (const ScalarType sval) {
                                             return sval < tval;
                                           }) + 1;
        ret.resize(num);
        for (int i=last-num+1; i<last; ++i) {
          ret.push_back(i);
        }
      }
      return ret;
    }
    //@}

    private: 
      SortType which_;
      bool absthresh_;
      typename Teuchos::ScalarTraits<ScalarType>::magnitudeType thresh_;
  };

  //! Composite filter
  template <class ScalarType>
  class CompFilter : public Filter<ScalarType> {
    public: 

    enum CompType {
      AND,
      OR
    };

    //! @name Constructor/Destructor.
    //@{

    //! Default constructor.
    CompFilter(CompType andor, 
                 const Teuchos::RCP<Filter<ScalarType> > &f1,
                 const Teuchos::RCP<Filter<ScalarType> > &f2 ) 
      : andor_(andor), f1_(f1), f2_(f2) {
      TEUCHOS_TEST_FOR_EXCEPTION(f1_ == Teuchos::null || f2_ == Teuchos::null,
                         std::invalid_argument,"CompFilter: Component filters must be non-null.");
    };

    //! Destructor.
    virtual ~CompFilter() {};
    //@}

    //! @name Filter Method
    //@{

    //! \brief 
    /*!  \note It is assumed that svals are sorted in decreasing order.
     */
      std::vector<int> filter(const std::vector<ScalarType> &svals) {
      std::vector<int> ind1 = f1_->filter(svals);
      std::vector<int> ind2 = f2_->filter(svals);
      std::vector<int> ret;
      if (AND == andor_) {
        set_intersection(ind1.begin(),ind1.end(),ind2.begin(),ind2.end(),ret.begin());
      }
      else if (OR == andor_) {
        set_union(ind1.begin(),ind1.end(),ind2.begin(),ind2.end(),ret.begin());
      }
      return ret;
    }
    //@}

    private: 
      CompType andor_;
      Teuchos::RCP<Filter<ScalarType> > f1_, f2_;
  };

}
#endif
