// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_PARAMETERLIBRARYBASE_HPP
#define SACADO_PARAMETERLIBRARYBASE_HPP

#include <iostream>

#include "Teuchos_Array.hpp"

#include "Sacado_ParameterFamilyBase.hpp"
#include "Sacado_ParameterVectorBase.hpp"
#include "Sacado_mpl_apply.hpp"

namespace Sacado {

  using std::string;

  /*!
   * \brief Class to provide a centralized library for setting/retrieving
   * numerical parameter values.
   */
  template <typename FamilyType, typename EntryType>
  class ParameterLibraryBase  {

  protected:

    //! Map of all parameter families
    typedef std::map<string, Teuchos::RCP<FamilyType> > FamilyMap;

  public:

    //! Iterator typename
    typedef typename FamilyMap::iterator iterator;

    //! Const iterator typename
    typedef typename FamilyMap::const_iterator const_iterator;

    //! Default constructor
    ParameterLibraryBase();

    //! Destructor
    virtual ~ParameterLibraryBase();

    //! Determine if parameter of name \em name is in the library
    bool isParameter(const std::string& name) const;

    //! Determine if parameter of name \em name has type \em type
    template <typename EvalType>
    bool isParameterForType(const std::string& name) const;

    //! Create a new parameter family
    /*!
     * Returns true if successful in adding family to library, false
     * otherwise.
     */
    bool addParameterFamily(const std::string& name, bool supports_ad,
                            bool supports_analytic);


    //! Add a new parameter using custom entry
    /*!
     * Returns true if successful in adding entry to library, false
     * otherwise.  If \c allow_overwrite is true, any existing entry will
     * be overwritten by the supplied entry.
     */
    template <typename EvalType>
    bool addEntry(const std::string& name,
                  const Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >& entry,
                  const bool allow_overwrite = false);

    //! Return parameter entry
    template <typename EvalType>
    Teuchos::RCP< typename Sacado::mpl::apply<EntryType,EvalType>::type >
    getEntry(const std::string& name);

    //! Return parameter entry
    template <typename EvalType>
    Teuchos::RCP< const typename Sacado::mpl::apply<EntryType,EvalType>::type >
    getEntry(const std::string& name) const;

    //! Return number of parameters in library
    unsigned int size() const { return library.size(); }

    //! Iterator pointing at beginning of library
    iterator begin() { return library.begin(); }

    //! Iterator pointing at beginning of library
    const_iterator begin() const { return library.begin(); }

    //! Iterator pointing at end of library
    iterator end() { return library.end(); }

    //! Iterator pointing at end of library
    const_iterator end() const { return library.end(); }

    //! Fill a vector with the supplied parameter names and values
    template <typename BaseValueType>
    void
    fillVector(const Teuchos::Array<std::string>& names,
               const Teuchos::Array<BaseValueType>& values,
               ParameterVectorBase<FamilyType,BaseValueType>& pv);

    //! Print parameter library
    /*!
     * Set print_values = true to print each parameter value for
     * each evaluation type.
     */
    void print(std::ostream& os, bool print_values = false) const;

    //! Clear the library
    void clear() { library = FamilyMap(); }

  private:

    //! Private to prohibit copying
    ParameterLibraryBase(const ParameterLibraryBase&);

    //! Private to prohibit copying
    ParameterLibraryBase& operator = (const ParameterLibraryBase&);

  protected:

    //! Scalar parameter library
    FamilyMap library;
  };

  template <typename FamilyType, typename EntryType>
  std::ostream&
  operator << (std::ostream& os,
               const ParameterLibraryBase<FamilyType, EntryType>& pl)
  {
    pl.print(os);
    return os;
  }
}

// Include template definitions
#include "Sacado_ParameterLibraryBaseImp.hpp"

#endif
