// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Belos_SolverFactory_hpp
#define __Belos_SolverFactory_hpp

#include <BelosConfigDefs.hpp>
#include <BelosOutputManager.hpp>
#include <BelosSolverManager.hpp>

#include "Belos_Details_EBelosSolverType.hpp"
#include "BelosCustomSolverFactory.hpp"

#include <Teuchos_Describable.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Belos {
namespace Impl {

//! Print the given array of strings, in YAML format, to \c out.
void
printStringArray (std::ostream& out,
                  const Teuchos::ArrayView<const std::string>& array);

//! Print the given array of strings, in YAML format, to \c out.
void
printStringArray (std::ostream& out,
                  const std::vector<std::string>& array);

//! Return the upper-case version of s.
std::string
upperCase (const std::string& s);

/// \brief Specializations of Belos::SolverFactory may inherit from
///   this class to get basic SolverFactory functionality.
///
/// This class is not for Belos users.  It's really just for Belos
/// developers who want to write a specialization of
/// Belos::SolverFactory.  Those developers may make their
/// specialization inherit from this class, in order to get basic
/// Belos::SolverFactory functionality without reimplementing
/// everything.
///
/// \tparam Scalar Same as template parameter 1 of
///   Belos::SolverFactory (which see below).
/// \tparam MV Same as template parameter 2 of
///   Belos::SolverFactory (which see below).
/// \tparam OP Same as template parameter 2 of
///   Belos::SolverFactory (which see below).
template<class Scalar, class MV, class OP>
class SolverFactoryParent :
    public Teuchos::Describable
{
protected:
  // SolverFactoryParent should never be created directly. Usually it will be
  // created through derived classes EpetraSolverFactory, TpetraSolverFactory,
  // BelosSolverFactory, or XpetraSolverFactory. If you are doing custom types
  // and include BelosSolverFactory_Generic.hpp you should explicitly use
  // GenericSolverFactory, not SolverFactory, which will avoid an error trying
  // to construct here. GenericSolverFactory is special because it registers
  // all the solver managers for any type. Note that if you are using hard coded
  // types it is possible that some the type sets will be connecting to an
  // automatic solver factory such as TpetraSolverFactory while another type
  // set could be going through the GenericSolverFactory.
  SolverFactoryParent() {}

public:
  /// \brief The type of the solver returned by create().
  ///
  /// This is a specialization of SolverManager for the same scalar,
  /// multivector, and operator types as the template parameters of
  /// this factory.
  typedef ::Belos::SolverManager<Scalar, MV, OP> solver_base_type;

  /// \brief The type of a solver factory that users may give to
  ///   addFactory() (which see below)
  typedef CustomSolverFactory<Scalar, MV, OP> custom_solver_factory_type;

protected:
  /// \brief Return an instance of the specified solver, or
  ///   Teuchos::null if this factory does not provide the requested
  ///   solver.
  ///
  /// The preferred way to customize this method is not to inherit
  /// from it, but rather to add a custom SolverFactory via
  /// addFactory() (which see below).  Nevertheless, we leave open the
  /// possibility of overriding this method, for example in order to
  /// change the order in which users' solver factories are queried.
  ///
  /// \note To implementers: DO NOT THROW if this factory does not
  ///   recognize and support the input solverName.  Instead, just
  ///   return Teuchos::null.  This will make Belos::SolverFactory's
  ///   implementation cleaner (and is also the reason why we gave
  ///   this function a different name).
  ///
  /// \param solverName [in] Name of the solver.
  ///
  /// \param solverParams [in/out] List of parameters with which to
  ///   configure the solver.  If null, we configure the solver with
  ///   default parameters.  If nonnull, the solver may modify the
  ///   list by filling in missing parameters with default values.
  ///   You can then inspect the resulting list to learn what
  ///   parameters the solver accepts.
  ///
  /// Some solvers may be accessed by multiple names ("aliases").
  /// Each solver has a canonical name, and zero or more aliases.
  /// Using some aliases (such as those that access Flexible GMRES
  /// capability in GMRES-type solvers) may make this method set
  /// certain parameters in your parameter list.
  ///
  /// The input parameter list is passed in as a Teuchos::RCP because
  /// the factory passes it to the solver, and Belos solvers want
  /// their input parameter list as a
  /// Teuchos::RCP<Teuchos::ParameterList>.  We allow a null parameter
  /// list only for convenience, and will use default parameter values
  /// in that case.
  virtual Teuchos::RCP<solver_base_type>
  getSolver (const std::string& solverName,
             const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

public:
  /// \brief Create, configure, and return the specified solver.
  ///
  /// \param solverName [in] Name of the solver.
  ///
  /// \param solverParams [in/out] List of parameters with which to
  ///   configure the solver.
  ///
  /// This method differs from getSolver() (see above) only in that it
  /// throws an exception if solverName is invalid.
  virtual Teuchos::RCP<solver_base_type>
  create (const std::string& solverName,
          const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

  /// \brief Number of supported solvers.
  ///
  /// This may differ from the number of supported solver
  /// <i>names</i>, since we may accept multiple names ("aliases") for
  /// some solvers.
  virtual int numSupportedSolvers () const;

  /// \brief List of supported solver names.
  ///
  /// The length of this list may differ from the number of supported
  /// solvers, since we may accept multiple names ("aliases") for some
  /// solvers.
  virtual Teuchos::Array<std::string> supportedSolverNames () const;

  //! Whether the given solver name names a supported solver.
  virtual bool isSupported (const std::string& solverName) const;

  /// \brief Add a custom solver factory.
  ///
  /// Any custom solver factories that you may define will override
  /// this factory.  "Override" means that if the \c solverName
  /// argument to create() or getSolver() matches any solvers that the
  /// custom factories support, then one of the custom factories will
  /// create it.
  ///
  /// \note To developers: This is an instance method, but it adds the
  ///   factory for all current and future instances.
  ///
  /// \warning This method makes no promise of reentrancy or thread
  ///   safety, with respect to other calls to this or other
  ///   factories' methods.
  void
  addFactory (const Teuchos::RCP<custom_solver_factory_type>& factory);

  /// \brief Clear all custom solver factories.
  ///
  /// Any custom solver factories that have been added to this factory.
  /// 
  /// \note To developers: Since a custom factory will be added to all 
  ///   current and future instances, there needs to be a method to clear
  ///   these factories to avoid memory leaks.
  void
  clearFactories ();

  /// \brief register a solver for Inverted Injection (DII).
  static void
  registerSolver (const std::string & solverName,
    Teuchos::RCP<SolverFactoryParent<Scalar, MV, OP>::solver_base_type > instance)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
      instance == Teuchos::null,
      std::invalid_argument, "Belos::SolverFactoryParent::registerSolver "
      "was given a null solver to register.");

    get_solverManagers()[solverName] = instance;
  }

  /// \brief is solver registered for Inverted Injection (DII).
  static bool
  isSolverRegistered (const std::string & solverName)
  {
    return (get_solverManagers().find(solverName) != get_solverManagers().end());
  }

  //! @name Implementation of Teuchos::Describable interface
  //@{

  //! A string description of this object.
  virtual std::string description() const;

  /// \brief Describe this object.
  ///
  /// At higher verbosity levels, this method will print out the list
  /// of names of supported solvers.  You can also get this list
  /// directly by using the supportedSolverNames() method.
  virtual void
  describe (Teuchos::FancyOStream& out,
            const Teuchos::EVerbosityLevel verbLevel =
            Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  //! The list of solver factories given to addFactory.
  static std::vector<Teuchos::RCP<custom_solver_factory_type> > factories_;

  static std::map<const std::string, Teuchos::RCP<typename
    SolverFactoryParent<Scalar, MV, OP>::solver_base_type> > &
    get_solverManagers() {
      static std::map<const std::string, Teuchos::RCP<typename
        SolverFactoryParent<Scalar, MV, OP>::solver_base_type> > solverManagers;
      return solverManagers;
    }
};

template<class Scalar, class MV, class OP>
std::vector<Teuchos::RCP<typename SolverFactoryParent<Scalar, MV, OP>::custom_solver_factory_type> >
SolverFactoryParent<Scalar, MV, OP>::factories_;

template<class SolverClass, class Scalar, class MV, class OP>
void registerSolverSubclassForTypes (const std::string & solverName) {
  if(!SolverFactoryParent<Scalar, MV, OP>::isSolverRegistered(solverName)) {
    Teuchos::RCP<SolverClass> solver (new SolverClass);
    SolverFactoryParent<Scalar, MV, OP>::registerSolver (solverName, solver);
  }
}

// specializations get a typedef "type"
// If this compile fails then the error is likely that BelosSolverFactory.hpp
// was included directly but the specific sub class of SolverFactoryParent was
// not included. Examples are:
//   BelosSolverFactory_Belos.hpp, BelosSolverFactory_Epetra.hpp,
//   BelosSolverFactory_Tpetra.hpp, BelosSolverFactory_Xpetra.hpp
// These were setup to be automatically included through the corresponding
// adapter includes so something may have gone wrong with that.
template<class SC, class MV, class OP>
class SolverFactorySelector {
  public:
    // TODO: This could be deleted except for the GenericSolverFactory which
    // needs to be declared for all types. So I added this but then if you
    // include GenericSolverFactory you will have SolverFactory simply point
    // to SolverFactoryParent. I changed that constructor to be protected so
    // using SolverFactory (pointing to SolverFactoryParent) will give a compile
    // error. For GenericSolverFactory you must explicity use GenericSolverFactory
    // in the code. This may be preferable because it makes it clear that
    // factory is not connecting to the standard set of types. I'm not sure how
    // to do this in a better way.
    typedef SolverFactoryParent<SC,MV,OP> type;
};

} // namespace Impl

// Derived setups such as found in BelosSolverFactory_Tpetra.hpp will define
// this specialization so that SolverFactory will be used as SolverFactoryTpetra.
template<class SC, class MV, class OP>
using SolverFactory = typename ::Belos::Impl::SolverFactorySelector<SC, MV, OP>::type;

namespace Impl {

template<class Scalar, class MV, class OP>
Teuchos::RCP<typename SolverFactoryParent<Scalar, MV, OP>::solver_base_type>
SolverFactoryParent<Scalar, MV, OP>::
create (const std::string& solverName,
        const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  using Teuchos::RCP;
  RCP<solver_base_type> solver = this->getSolver (solverName, solverParams);
  TEUCHOS_TEST_FOR_EXCEPTION
    (solver.is_null (), std::invalid_argument,
     "Invalid or unsupported Belos solver name \"" << solverName << "\".");
  return solver;
}

template<class Scalar, class MV, class OP>
Teuchos::RCP<typename SolverFactoryParent<Scalar, MV, OP>::solver_base_type>
SolverFactoryParent<Scalar, MV, OP>::
getSolver (const std::string& solverName,
           const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  using Teuchos::RCP;

  // First, check the overriding factories.
  for (std::size_t k = 0; k < factories_.size (); ++k) {
    RCP<CustomSolverFactory<Scalar, MV, OP> > factory = factories_[k];
    if (! factory.is_null ()) {
      RCP<SolverManager<Scalar, MV, OP> > solver =
        factory->getSolver (solverName, solverParams);
      if (! solver.is_null ()) {
        return solver;
      }
    }
  }

  // Upper-case version of the input solver name.
  const std::string solverNameUC = Impl::upperCase (solverName);

  // Check whether the given name is an alias.
  std::pair<std::string, bool> aliasResult =
    Details::getCanonicalNameFromAlias (solverNameUC);
  const std::string candidateCanonicalName = aliasResult.first;
  const bool isAnAlias = aliasResult.second;

  // Get the standardized name for enum reference and map
  std::string standardized_name = isAnAlias ?
                                  candidateCanonicalName :
                                  solverNameUC;

  // If the input list is null, we create a new list and use that.
  // This is OK because the effect of a null parameter list input is
  // to use default parameter values.  Thus, we can always replace a
  // null list with an empty list.
  Teuchos::RCP<Teuchos::ParameterList> pl =
    solverParams.is_null() ? Teuchos::parameterList() : solverParams;

  // Possibly modify the input parameter list as needed.
  if (isAnAlias) {
    Details::reviseParameterListForAlias (solverNameUC, *pl);
  }

  typename std::map<const std::string, Teuchos::RCP<
    typename SolverFactoryParent<Scalar, MV, OP>::solver_base_type> >::iterator
    it = get_solverManagers().find (standardized_name);

  TEUCHOS_TEST_FOR_EXCEPTION(
    it == get_solverManagers().end(),
    std::invalid_argument, "Belos solver manager " << solverNameUC <<
    " with standardized name " << standardized_name << " has not been"
    " registered.");

  TEUCHOS_TEST_FOR_EXCEPTION(
    it->second == Teuchos::null,
    std::logic_error, "Belos::SolverFactoryParent: The registered "
    "clone source for " << solverNameUC << " with standardized name "
    << standardized_name << " is null which should never happen."
    ".  Please report this bug to the Belos developers.");

  // clone the solver
  RCP<solver_base_type> solver = (it->second)->clone ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    solver == Teuchos::null,
    std::logic_error, "Belos::SolverFactoryParent: Failed "
    "to clone SolverManager with name " << solverNameUC << " with standardized"
    " name" << standardized_name << "."
    ".  Please report this bug to the Belos developers.");

  // Some solvers may not like to get a null ParameterList.  If params
  // is null, replace it with an empty parameter list.  The solver
  // will fill in default parameters for that case.  Use the name of
  // the solver's default parameters to name the new empty list.
  if (pl.is_null()) {
    pl = Teuchos::parameterList (solver->getValidParameters ()->name ());
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    pl.is_null(), std::logic_error,
    "Belos::SolverFactory: ParameterList to pass to solver is null.  This "
    "should never happen.  Please report this bug to the Belos developers.");
  solver->setParameters (pl);
  return solver;
}


template<class Scalar, class MV, class OP>
void
SolverFactoryParent<Scalar, MV, OP>::
addFactory (const Teuchos::RCP<CustomSolverFactory<Scalar, MV, OP> >& factory)
{
  factories_.push_back (factory);
}


template<class Scalar, class MV, class OP>
void
SolverFactoryParent<Scalar, MV, OP>::
clearFactories ()
{
  factories_.clear();
}


template<class Scalar, class MV, class OP>
std::string
SolverFactoryParent<Scalar, MV, OP>::
description () const
{
  using Teuchos::TypeNameTraits;

  std::ostringstream out;
  out << "\"Belos::SolverFactory\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: " << this->getObjectLabel () << ", ";
  }
  out << "Scalar: \"" << TypeNameTraits<Scalar>::name ()
      << "\", MV: \"" << TypeNameTraits<MV>::name ()
      << "\", OP: \"" << TypeNameTraits<OP>::name ()
      << "\"}";
  return out.str ();
}


template<class Scalar, class MV, class OP>
void
SolverFactoryParent<Scalar, MV, OP>::
describe (Teuchos::FancyOStream& out,
          const Teuchos::EVerbosityLevel verbLevel) const
{
  using Teuchos::TypeNameTraits;
  using std::endl;

  const Teuchos::EVerbosityLevel vl =
    (verbLevel == Teuchos::VERB_DEFAULT) ? Teuchos::VERB_LOW : verbLevel;

  if (vl == Teuchos::VERB_NONE) {
    return;
  }

  // By convention, describe() always begins with a tab.
  Teuchos::OSTab tab0 (out);
  // The description prints in YAML format.  The class name needs to
  // be protected with quotes, so that YAML doesn't get confused
  // between the colons in the class name and the colon separating
  // (key,value) pairs.
  out << "\"Belos::SolverFactory\":" << endl;
  if (this->getObjectLabel () != "") {
    out << "Label: " << this->getObjectLabel () << endl;
  }
  {
    out << "Template parameters:" << endl;
    Teuchos::OSTab tab1 (out);
    out << "Scalar: \"" << TypeNameTraits<Scalar>::name () << "\"" << endl
        << "MV: \"" << TypeNameTraits<MV>::name () << "\"" << endl
        << "OP: \"" << TypeNameTraits<OP>::name () << "\"" << endl;
  }

  // At higher verbosity levels, print out the list of supported solvers.
  if (vl > Teuchos::VERB_LOW) {
    Teuchos::OSTab tab1 (out);
    out << "Number of solvers: " << numSupportedSolvers ()
        << endl;
    out << "Canonical solver names: ";
    Impl::printStringArray (out, Details::canonicalSolverNames ());
    out << endl;

    out << "Aliases to canonical names: ";
    Impl::printStringArray (out, Details::solverNameAliases ());
    out << endl;
  }
}

template<class Scalar, class MV, class OP>
int
SolverFactoryParent<Scalar, MV, OP>::
numSupportedSolvers () const
{
  int numSupported = 0;

  // First, check the overriding factories.
  for (std::size_t k = 0; k < factories_.size (); ++k) {
    using Teuchos::RCP;
    RCP<custom_solver_factory_type> factory = factories_[k];
    if (! factory.is_null ()) {
      numSupported += factory->numSupportedSolvers ();
    }
  }

  // Now, see how many solvers this factory supports.
  return numSupported + Details::numSupportedSolvers ();
}

template<class Scalar, class MV, class OP>
Teuchos::Array<std::string>
SolverFactoryParent<Scalar, MV, OP>::
supportedSolverNames () const
{
  typedef std::vector<std::string>::const_iterator iter_type;
  Teuchos::Array<std::string> names;

  // First, check the overriding factories.
  const std::size_t numFactories = factories_.size ();
  for (std::size_t factInd = 0; factInd < numFactories; ++factInd) {
    Teuchos::RCP<custom_solver_factory_type> factory = factories_[factInd];
    if (! factory.is_null ()) {
      std::vector<std::string> supportedSolvers =
        factory->supportedSolverNames ();
      const std::size_t numSolvers = supportedSolvers.size ();
      for (std::size_t solvInd = 0; solvInd < numSolvers; ++solvInd) {
        names.push_back (supportedSolvers[solvInd]);
      }
    }
  }

  {
    std::vector<std::string> aliases = Details::solverNameAliases ();
    for (iter_type iter = aliases.begin (); iter != aliases.end (); ++iter) {
      names.push_back (*iter);
    }
  }
  {
    std::vector<std::string> canonicalNames = Details::canonicalSolverNames ();
    for (iter_type iter = canonicalNames.begin ();
         iter != canonicalNames.end (); ++iter) {
      names.push_back (*iter);
    }
  }
  return names;
}

template<class Scalar, class MV, class OP>
bool
SolverFactoryParent<Scalar, MV, OP>::
isSupported (const std::string& solverName) const
{
  // First, check the overriding factories.
  const std::size_t numFactories = factories_.size ();
  for (std::size_t factInd = 0; factInd < numFactories; ++factInd) {
    using Teuchos::RCP;
    RCP<custom_solver_factory_type> factory = factories_[factInd];
    if (! factory.is_null ()) {
      if (factory->isSupported (solverName)) {
        return true;
      }
    }
  }
  // Now, check this factory.

  // Upper-case version of the input solver name.
  const std::string solverNameUC = Impl::upperCase (solverName);

  // Check whether the given name is an alias.
  std::pair<std::string, bool> aliasResult =
    Details::getCanonicalNameFromAlias (solverNameUC);
  const std::string candidateCanonicalName = aliasResult.first;
  const bool validCanonicalName =
    (get_solverManagers().find(candidateCanonicalName) != get_solverManagers().end());
  return validCanonicalName;
}

} // namespace Impl
} // namespace Belos

// We have things like BelosSolverFactory_Tpetra.hpp which are automatically
// included through the adapter includes to maintain backwards compatibility.
// The Belos version itself doesn't have a place like that so included here
// which is awkward. It might make more sense to just copy that code here but
// it has symmetry with the other files and wanted to preserve that. To discuss.
#include "BelosSolverFactory_Belos.hpp"

#endif // __Belos_SolverFactory_hpp

