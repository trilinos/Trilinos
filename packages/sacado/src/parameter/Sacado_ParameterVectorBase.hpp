// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_PARAMETERVECTORBASE_HPP
#define SACADO_PARAMETERVECTORBASE_HPP

#include <vector>

#include "Teuchos_Array.hpp"

#include "Sacado_ParameterFamilyBase.hpp"

namespace Sacado {

  /*!
   * \brief A class to store the active parameters in a code in an ordered
   * fashion, along with their "base" values, i.e., the floating point
   * value upon which the templated values are based.
   */
  template <typename FamilyType, typename BaseValueType>
  class ParameterVectorBase  {

  public:

    //! Container for parameter entries
    struct Entry {

      //! Pointer to family
      Teuchos::RCP<FamilyType> family;

      //! Base value of parameter family
      BaseValueType baseValue;

      //! Constructor
      Entry(const Teuchos::RCP<FamilyType>& f, BaseValueType bv) :
        family(f), baseValue(bv) {}

    };

  protected:

    //! Vector of all parameter families
    typedef Teuchos::Array<Entry> EntryVector;

  public:

    //! Iterator typename
    typedef typename EntryVector::iterator iterator;

    //! Const iterator typename
    typedef typename EntryVector::const_iterator const_iterator;

    //! Default constructor
    ParameterVectorBase() {}

    //! Copy constructor
    ParameterVectorBase(const ParameterVectorBase& source) :
      params(source.params) {}

    //! Destructor
    virtual ~ParameterVectorBase() {}

    //! Assignment
    ParameterVectorBase& operator = (const ParameterVectorBase& source) {
      params = source.params; return *this; }

    //! Add entry
    void addParam(const Teuchos::RCP<FamilyType>& family,
                  BaseValueType baseValue) {
      params.push_back(Entry(family, baseValue));
    }

    //! Return number of parameters in vector
    unsigned int size() const { return params.size(); }

    //! Element access
    Entry& operator[] (int i) { return params[i]; }

    //! Element access
    const Entry& operator[] (int i) const { return params[i]; }

    //! Iterator pointing at beginning of vector
    iterator begin() { return params.begin(); }

    //! Iterator pointing at beginning of vector
    const_iterator begin() const { return params.begin(); }

    //! Iterator pointing at end of vector
    iterator end() { return params.end(); }

    //! Iterator pointing at end of vector
    const_iterator end() const { return params.end(); }

    //! Filter vector into types
    void
    filterParameters(ParameterVectorBase& ad,
                     ParameterVectorBase& analytic,
                     ParameterVectorBase& other,
                     std::vector<int>& index_ad,
                     std::vector<int>& index_analytic,
                     std::vector<int>& index_other) {
      index_ad.resize(0);
      index_analytic.resize(0);
      index_other.resize(0);

      typename EntryVector::iterator it;
      int i;
      for (it = params.begin(), i=0; it != params.end(); ++it, ++i) {
        if ((*it).family->supportsAD()) {
          ad.params.push_back(*it);
          index_ad.push_back(i);
        }
        else if ((*it).family->supportsAnalytic()) {
          analytic.params.push_back(*it);
          index_analytic.push_back(i);
        }
        else {
          other.params.push_back(*it);
          index_other.push_back(i);
        }
      }
    }

  protected:

    //! Parameter vector
    EntryVector params;
  };
}

#endif
