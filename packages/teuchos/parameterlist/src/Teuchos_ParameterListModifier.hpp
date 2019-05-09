// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_ParameterListModifier_hpp_
#define TEUCHOS_ParameterListModifier_hpp_

/*! \file Teuchos_ParameterListModifier.hpp
    \brief Parameter List Modifier class
*/

#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterEntryValidator.hpp"
#include "Teuchos_map.hpp"
#include "Teuchos_RCP.hpp"
namespace Teuchos {


#ifndef DOXYGEN_SHOULD_SKIP_THIS
class ParameterList;
#endif

/** \brief Abstract interface for an object that can modify both a
 *  parameter list and the parameter list being used during the
 *  validation stage.
 *
 * A parameter (sub)list modifier can be used to process optional fields and
 * dependent fields before validation.  It can also be used after validation
 * to reconcile parameters that may have dependencies on other parameters.
 */
class TEUCHOSPARAMETERLIST_LIB_DLL_EXPORT ParameterListModifier : public Describable {

public:

  //! @name Constructors/Destructor
  //@{

  //! Constructor
  ParameterListModifier() = default;

  //! Constructor that names the modifier.
  ParameterListModifier(const std::string &name);

  //! Destructor
  virtual ~ParameterListModifier();

  //@}
  //! @name Set Functions
  //@{
  //! Set the name of <tt>*this</tt> modifier.
  ParameterListModifier& setName( const std::string &name );

  //@}

  /** \brief Print documentation for this parameter list modifier.
   *
   * \param docString [in] (Multi-line) documentation std::string.
   *
   * \param out [out] The std::ostream used for the output
   *
   * The purpose of this function is to augment what is
   * in <tt>docString</tt>
   * with some description of what happens during the modification and
   * reconciliation stages of this modifier.
   */
  void printDoc(
    std::string const& docString,
    std::ostream &out
    ) const;


  /** \brief Find the parameters and/or sublists with matching base names.
   *
   * \param paramList [in] Modified parameter list to search.
   *
   * \param findParameters [out] Search through parameters
   *
   * \param findSublists [out] Search through sublists (not recursive)
   *
   * This convenience function makes it easy to search through the current level of
   * a given parameter list and find all parameters and/or sublists that begin with a
   * given name.
   */
  Array<std::string> findMatchingBaseNames(const ParameterList &paramList,
      const std::string &baseName, const bool &findParameters = true,
      const bool &findSublists = true) const;

  /** \brief Modify a parameter list and/or the valid parameter list being used to validate
   * it and throw std::exception (with a great error message) if modification fails.
   *
   * \param  paramList
   *            [in] The parameter list that needs to be validated
   * \param  validParamList
   *            [in] The parameter list being used as a template for validation.
   *
   * This function is usually called before the validation step begins in order to create optional
   * parameters and/or sublists.
   */
  virtual void modify(ParameterList &paramList, ParameterList &validParamList) const {};


  /** \brief Reconcile a parameter list and/or the valid parameter list being used to validate
   * it and throw std::exception (with a great error message) if reconciliation fails.
   *
   * \param  paramList
   *            [in,out] The parameter list that needs to be validated
   * \param  validParamList
   *            [in,out] The parameter list being used as a template for validation.
   *
   * This function is usually called after the validation step begins in order to check that
   * dependencies between multiple parameters are satisfied.
   */
  virtual void reconcile(ParameterList &paramList) const {};


  /** \brief Create parameters in the valid parameter list using a template parameter from the valid
   *  parameter list and the names of parameters in the list being validated.
   *
   *  \param parameterTemplateName [in] The name of the parameter template in \a valid_pl
   *
   *  \param paramList [in] The parameter list that needs to be validated
   *
   *  \param validParamList [in,out] The parameter list that is being used as a template for validation
   *
   *  \param excludeParameters [in] An optional list of parameter names to exclude
   *
   */
  int expandParameters(const std::string &paramTemplateName, ParameterList &paramList,
      ParameterList &validParamList,
      const Array<std::string> &excludeParameters=Array<std::string>()) const;


  /** \brief Create sublists in the valid parameter list using a template parameter from the valid
   *  parameter list and the names of sublists in the list being validated.
   *
   *  \param sublistTemplateName [in] The name of the sublist template in \a valid_pl
   *
   *  \param paramList [in] The parameter list that needs to be validated
   *
   *  \param validParamList [in,out] The parameter list that is being used as a template for validation
   *
   *  \param excludeSublists [in] An optional list of sublist names to exclude
   *
   */
  int expandSublists(const std::string &sublistTemplateName, ParameterList &paramList,
      ParameterList &validParamList,
      const Array<std::string> &excludeSublists=Array<std::string>()) const;


  /** \brief Copy a parameter into desired sublists.
   *
   * \param paramName [in] The name of the parameter to be copied.
   *
   * \param paramList [in,out] The parameter list with \a paramName.
   *
   * \param sublistNames [in] The names of any sublists in \a paramList to set the defaults in
   * using parameter \a paramName.
   *
   * \param removeParam [in] Remove \a paramName from \a paramList after defaults are set in sublists.
   */
  int setDefaultsInSublists(const std::string &paramName, ParameterList &paramList,
      const Array<std::string> &sublistNames, const bool removeParam = true) const;


  /** \brief Create sublists in the valid parameter list using a base name and the corresponding sublists
   *  in the parameter list being validated.
   *
   *  \param baseName [in] The root name of the sublists to look for and create
   *
   *  \param paramList [in] The parameter list that needs to be validated
   *
   *  \param validParamList [in,out] The parameter list that is being used as a template for validation
   *
   *  \param allowBaseName [in] Allow the parameter list \a pl to contain a parameter with the same name
   *    as base_name.
   *
   */
  int expandSublistsUsingBaseName(const std::string &baseName, ParameterList &paramList,
      ParameterList &validParamList, const bool &allowBaseName = true) const;


protected:
  std::string name_ = "ANONYMOUS";

};


// /////////////////////////////////////////////////////
// Inline and Template Function Definitions


inline
ParameterListModifier& ParameterListModifier::setName( const std::string &name_in )
{
  name_ = name_in;
  return *this;
}


} // end of Teuchos namespace

#endif
