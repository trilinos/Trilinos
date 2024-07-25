// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_InverseLibrary_hpp__
#define __Teko_InverseLibrary_hpp__

#include "Teko_InverseFactory.hpp"

// Teuchos includes
#include "Teuchos_ParameterList.hpp"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Teko includes
#include "Teko_RequestHandler.hpp"
#include "Teko_RequestHandlerContainer.hpp"

namespace Teko {

void addToStratimikosBuilder(const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& builder);

class InverseLibrary : public RequestHandlerContainer {
 public:
  InverseLibrary();

  InverseLibrary(const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& strat);

  //! add an unspecified inverse to the library
  void addInverse(const std::string& label, const Teuchos::ParameterList& pl);

  //! Add a Stratimikos solver with a label to the library
  void addStratSolver(const std::string& label, const std::string& type,
                      const Teuchos::ParameterList& pl);

  //! Add a Stratimikos preconditioner with a label to the library
  void addStratPrecond(const std::string& label, const std::string& type,
                       const Teuchos::ParameterList& pl);

  //! Add a Teko preconditioner to the library with a label
  void addBlockPrecond(const std::string& label, const std::string& type,
                       const Teuchos::ParameterList& pl);

  /** Get the fully constructed parameter list for a particular label
   *
   * \param[in] label Name used for the desired solver.
   *
   * \returns If the label is found in the library the corresponding parameter list
   *          is returned, otherwise <code>Teuchos::null</code> is returned.
   */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList(const std::string& label) const;

  //! Get the inverse factory associated with a particular label
  Teuchos::RCP<InverseFactory> getInverseFactory(const std::string& label) const;

  //! Print the inverses and parameter lists available for use
  void PrintAvailableInverses(std::ostream& os) const;

  //! Set the request handler with pointers to the appropriate callbacks
  void setRequestHandler(const Teuchos::RCP<RequestHandler>& rh) { callbackHandler_ = rh; }

  //! Get the request handler with pointers to the appropriate callbacks
  Teuchos::RCP<RequestHandler> getRequestHandler() const { return callbackHandler_; }

 protected:
  // stratimikos type Inverse objects: mapping the label to a parameter list
  std::map<std::string, Teuchos::RCP<const Teuchos::ParameterList> > stratSolver_;
  std::map<std::string, Teuchos::RCP<const Teuchos::ParameterList> > stratPrecond_;
  std::map<std::string, Teuchos::RCP<const Teuchos::ParameterList> > blockPrecond_;

  // vectors showing which string types are in Stratimikos
  std::vector<std::string> stratValidSolver_;
  std::vector<std::string> stratValidPrecond_;
  std::vector<std::string> blockValidPrecond_;

  //! For handling requests and send requests back to the user
  Teuchos::RCP<RequestHandler> callbackHandler_;

  //! This is the default builder used by stratimikos
  Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder> defaultBuilder_;

 public:
  /** \brief Build an inverse library from a parameter list.
   *
   * Build an inverse library from a parameter list. This will
   * contain all the labeled inverses specified.
   *
   * \param[in] pl Parameter list to build the library from
   * \param[in] useStratDefaults Also add the default parameters from Stratimikos
   *
   * \returns A pointer to the inverse library created.
   */
  static Teuchos::RCP<InverseLibrary> buildFromParameterList(const Teuchos::ParameterList& pl,
                                                             bool useStratDefaults = true);

  /** \brief Build an inverse library from a parameter list.
   *
   * Build an inverse library from a parameter list. This will
   * contain all the labeled inverses specified.
   *
   * \param[in] pl Parameter list to build the library from
   * \param[in] strat Stratimikos object to use
   *
   * \returns A pointer to the inverse library created.
   */
  static Teuchos::RCP<InverseLibrary> buildFromParameterList(
      const Teuchos::ParameterList& pl,
      const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& strat);

  /** \brief Build an inverse library from Stratimikos
   *
   * Build an inverse library from Stratimkos. The labels
   * will just be the names in Stratimikos. Uses the Stratimikos
   * default linear solver builder and adds extra inverse types
   *
   * \returns A pointer to the inverse library created.
   */
  static Teuchos::RCP<InverseLibrary> buildFromStratimikos();

  /** \brief Build an inverse library from Stratimikos
   *
   * Build an inverse library from Stratimkos. The labels
   * will just be the names in Stratimikos.
   *
   * \param[in] strat Stratimikos object to use
   *
   * \returns A pointer to the inverse library created.
   */
  static Teuchos::RCP<InverseLibrary> buildFromStratimikos(
      const Stratimikos::DefaultLinearSolverBuilder& strat);

  /** \brief Build an inverse library from Stratimikos
   *
   * Build an inverse library from Stratimkos. The labels
   * will just be the names in Stratimikos.
   *
   * \param[in] strat Stratimikos pointer to use
   *
   * \returns A pointer to the inverse library created.
   */
  static Teuchos::RCP<InverseLibrary> buildFromStratimikos(
      const Teuchos::RCP<Stratimikos::DefaultLinearSolverBuilder>& strat);
};

}  // end namespace Teko

#endif
