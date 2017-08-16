//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// ************************************************************************
//@HEADER

#ifndef __Belos_SolverFactory_hpp
#define __Belos_SolverFactory_hpp

#include <BelosConfigDefs.hpp>
#include <BelosOutputManager.hpp>
#include <BelosSolverManager.hpp>

#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosGCRODRSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosPseudoBlockStochasticCGSolMgr.hpp>
#include <BelosLSQRSolMgr.hpp>
#include <BelosMinresSolMgr.hpp>
#include <BelosGmresPolySolMgr.hpp>
#include <BelosPCPGSolMgr.hpp>
#include <BelosRCGSolMgr.hpp>
#include <BelosTFQMRSolMgr.hpp>
#include <BelosPseudoBlockTFQMRSolMgr.hpp>
#include <BelosFixedPointSolMgr.hpp>
#include <BelosBiCGStabSolMgr.hpp>

#include "Belos_Details_EBelosSolverType.hpp"

#include <Teuchos_Describable.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <algorithm>
#include <locale>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Belos {

/// \class SolverFactory
/// \brief Factory for all solvers which Belos supports.
/// \author Mark Hoemmen
///
/// New Belos users should start by creating an instance of this
/// class, and using it to create the solver they want.
///
/// Belos implements several different iterative solvers.  The usual
/// way in which users interact with these solvers is through
/// appropriately named subclasses of \c SolverManager.  This factory
/// class tells users which solvers are supported.  It can initialize
/// and return any supported subclass of \c SolverManager, given a
/// short name of the subclass (such as "GMRES" or "CG").
///
/// Users ask for the solver they want by a string name, and supply an
/// optional (but recommended) list of parameters
/// (Teuchos::ParameterList) for the solver.  The solver may fill in
/// the parameter list with all the valid parameters and their default
/// values, which users may later inspect and modify.  Valid solver
/// names include both "canonical names" (each maps one-to-one to a
/// specific SolverManager subclass) and "aliases."  Some aliases are
/// short nicknames for canonical names, like "GMRES" for "Pseudoblock
/// GMRES".  Other aliases refer to a canonical solver name, but also
/// modify the user's parameter list.  For example, "Flexible GMRES"
/// is an alias for "Block GMRES", and also sets the "Flexible Gmres"
/// parameter to true in the input parameter list.
///
/// <table>
/// <caption> Mapping of solver names and aliases to Belos classes </caption>
/// <tr><th> Solver name </th>               <th> Aliases </th>                                                        <th> \c SolverManager subclass </th></tr>
/// <tr><td> Pseudoblock GMRES </td>         <td> GMRES, Pseudo Block GMRES, PseudoBlockGMRES, PseudoBlockGmres </td>  <td> \c PseudoBlockGmresSolMgr </td></tr>
/// <tr><td> Block GMRES </td>               <td> Flexible GMRES </td>                                                 <td> \c BlockGmresSolMgr </td></tr>
/// <tr><td> Block CG </td>                  <td> Block CG </td>                                                       <td> \c BlockCGSolMgr </td></tr>
/// <tr><td> Pseudoblock CG </td>            <td> PseudoBlockCG, Pseudo Block CG </td>                                 <td> \c PseudoBlockCGSolMgr </td></tr>
/// <tr><td> Pseudoblock Stochastic CG </td> <td> Stochastic CG </td>                                                  <td> \c PseudoBlockStochasticCGSolMgr </td></tr>
/// <tr><td> GCRODR </td>                    <td> Recycling GMRES </td>                                                <td> \c GCRODRSolMgr </td></tr>
/// <tr><td> RCG </td>                       <td> Recycling CG </td>                                                   <td> \c RCGSolMgr </td></tr>
/// <tr><td> MINRES </td>                    <td> MINRES </td>                                                         <td> \c MinresSolMgr </td></tr>
/// <tr><td> LSQR </td>                      <td> LSQR </td>                                                           <td> \c LSQRSolMgr </td></tr>
/// <tr><td> TFQMR </td>                     <td> TFQMR, Transpose-Free QMR </td>                                      <td> \c TFQMRSolMgr </td></tr>
/// <tr><td> Pseudoblock TFQMR </td>         <td> Pseudoblock TFQMR, Pseudo Block Transpose-Free QMR </td>             <td> \c PseudoBlockTFQMRSolMgr </td></tr>
/// <tr><td> Hybrid Block GMRES </td>        <td> GmresPoly, Seed GMRES </td>                                          <td> \c GmresPolySolMgr </td></tr>
/// <tr><td> PCPG </td>                      <td> CGPoly, Seed CG </td>                                                <td> \c PCPGSolMgr </td></tr>
/// </table>
///
/// This class' template parameters are the same as those of
/// Belos::SolverManager.  Scalar is the scalar type (of entries in
/// the multivector), MV is the multivector type, and OP is the
/// operator type.  For example: Scalar=double, MV=Epetra_MultiVector,
/// and OP=Epetra_Operator will access the Epetra specialization of
/// the Belos solvers.
///
/// Here is a simple example of how to use SolverFactory to create a
/// GMRES solver for your linear system.  Your code needs to include
/// BelosSolverFactory.hpp and whatever linear algebra library header
/// files you would normally use.  Suppose that Scalar, MV, and OP
/// have been previously typedef'd to the scalar resp. multivector
/// resp. operator type in your application.
/// \code
/// using Teuchos::ParameterList;
/// using Teuchos::parameterList;
/// using Teuchos::RCP;
/// using Teuchos::rcp; // Save some typing
///
/// // The ellipses represent the code you would normally use to create
/// // the sparse matrix, preconditioner, right-hand side, and initial
/// // guess for the linear system AX=B you want to solve.
/// RCP<OP> A = ...; // The sparse matrix / operator A
/// RCP<OP> M = ...; // The (right) preconditioner M
/// RCP<MV> B = ...; // Right-hand side of AX=B
/// RCP<MV> X = ...; // Initial guess for the solution
///
/// Belos::SolverFactory<Scalar, MV, OP> factory;
/// // Make an empty new parameter list.
/// RCP<ParameterList> solverParams = parameterList();
///
/// // Set some GMRES parameters.
/// //
/// // "Num Blocks" = Maximum number of Krylov vectors to store.  This
/// // is also the restart length.  "Block" here refers to the ability
/// // of this particular solver (and many other Belos solvers) to solve
/// // multiple linear systems at a time, even though we are only solving
/// // one linear system in this example.
/// solverParams->set ("Num Blocks", 40);
/// solverParams->set ("Maximum Iterations", 400);
/// solverParams->set ("Convergence Tolerance", 1.0e-8);
///
/// // Create the GMRES solver.
/// RCP<Belos::SolverManager<Scalar, MV, OP> > solver =
///   factory.create ("GMRES", solverParams);
///
/// // Create a LinearProblem struct with the problem to solve.
/// // A, X, B, and M are passed by (smart) pointer, not copied.
/// RCP<Belos::LinearProblem<Scalar, MV, OP> > problem =
///   rcp (new Belos::LinearProblem<Scalar, MV, OP> (A, X, B));
/// problem->setRightPrec (M);
///
/// // Tell the solver what problem you want to solve.
/// solver->setProblem (problem);
///
/// // Attempt to solve the linear system.  result == Belos::Converged
/// // means that it was solved to the desired tolerance.  This call
/// // overwrites X with the computed approximate solution.
/// Belos::ReturnType result = solver->solve();
///
/// // Ask the solver how many iterations the last solve() took.
/// const int numIters = solver->getNumIters();
/// \endcode
///
/// Belos developers who have implemented a new solver (i.e., a new
/// subclass of SolverManager) and who want to make the solver
/// available through the factory should do the following:
///
/// <ol>
/// <li> Add a new symbol corresponding to their solver to the
///      Details::EBelosSolverType enum. </li>
/// <li> If necessary, specialize Details::makeSolverManagerTmpl for
///      their SolverManager subclass.  In most cases, the default
///      implementation suffices. </li>
/// <li> Add a case for their enum symbol that instantiates their
///      solver to the long switch-case statement in
///      Details::makeSolverManagerFromEnum. </li>
/// <li> In the SolverFactory constructor, define a canonical string
///      name for their solver and its mapping to the corresponding
///      enum value, following the examples and comments there.  (This
///      takes one line of code.) </li>
/// </ol>
///
template<class Scalar, class MV, class OP>
class SolverFactory : public Teuchos::Describable {
public:
  /// \brief The type of the solver returned by create().
  ///
  /// This is a specialization of SolverManager for the same scalar,
  /// multivector, and operator types as the template parameters of
  /// this factory.
  typedef SolverManager<Scalar, MV, OP> solver_base_type;

  //! Default constructor.
  SolverFactory ();

  /// \brief Create, configure, and return the specified solver.
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
  Teuchos::RCP<solver_base_type>
  create (const std::string& solverName,
          const Teuchos::RCP<Teuchos::ParameterList>& solverParams);

  /// \brief Number of supported solvers.
  ///
  /// This may differ from the number of supported solver
  /// <i>names</i>, since we may accept multiple names ("aliases") for
  /// some solvers.
  int numSupportedSolvers () const;

  /// \brief List of supported solver names.
  ///
  /// The length of this list may differ from the number of supported
  /// solvers, since we may accept multiple names ("aliases") for some
  /// solvers.
  Teuchos::Array<std::string> supportedSolverNames () const;

  //! Whether the given solver name names a supported solver.
  bool isSupported (const std::string& solverName) const;

  //! @name Implementation of Teuchos::Describable interface
  //@{

  //! A string description of this object.
  std::string description() const;

  /// \brief Describe this object.
  ///
  /// At higher verbosity levels, this method will print out the list
  /// of names of supported solvers.  You can also get this list
  /// directly by using the supportedSolverNames() method.
  void describe (Teuchos::FancyOStream& out,
                 const Teuchos::EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  //! Print the given array of strings, in YAML format, to \c out.
  static void
  printStringArray (std::ostream& out,
                    const Teuchos::ArrayView<const std::string>& array)
  {
    typedef Teuchos::ArrayView<std::string>::const_iterator iter_type;

    out << "[";
    for (iter_type iter = array.begin(); iter != array.end(); ++iter) {
      out << "\"" << *iter << "\"";
      if (iter + 1 != array.end()) {
        out << ", ";
      }
    }
    out << "]";
  }

  //! Print the given array of strings, in YAML format, to \c out.
  static void
  printStringArray (std::ostream& out,
                    const std::vector<std::string>& array)
  {
    typedef std::vector<std::string>::const_iterator iter_type;

    out << "[";
    for (iter_type iter = array.begin(); iter != array.end(); ++iter) {
      out << "\"" << *iter << "\"";
      if (iter + 1 != array.end()) {
        out << ", ";
      }
    }
    out << "]";
  }
};


namespace Details {

/// \fn makeSolverManagerTmpl
/// \brief Return a new instance of the desired SolverManager subclass.
///
/// This template function is meant to be used only by \c
/// makeSolverManagerFromEnum.  We separate it out from \c
/// makeSolverManagerFromEnum in order to avoid duplicated code for
/// instantiating different \c SolverManager subclasses with the same
/// syntax (but different template parameters).
///
/// \tparam SolverManagerBaseType A specialization of SolverManager.
///
/// \tparam SolverManagerType The specific SolverManager subclass to
///   create.  It should take the same three template parameters
///   (Scalar, MV, OP) as SolverManagerBaseType.
///
/// \param params [in/out] List of parameters with which to configure
///   the solver.  If null, we configure the solver with default
///   parameters.
template<class SolverManagerBaseType, class SolverManagerType>
Teuchos::RCP<SolverManagerBaseType>
makeSolverManagerTmpl (const Teuchos::RCP<Teuchos::ParameterList>& params);

/// \fn makeSolverManagerFromEnum
/// \brief Return a new instance of the desired SolverManager subclass.
/// \author Mark Hoemmen
///
/// The \c SolverFactory class may use this template function
/// in order to instantiate an instance of the desired subclass of \c
/// SolverManager.
///
/// \tparam Scalar The first template parameter of \c SolverManager.
/// \tparam MV The second template parameter of \c SolverManager.
/// \tparam OP The third template parameter of \c SolverManager.
///
/// \param solverType [in] Enum value representing the specific
///   SolverManager subclass to instantiate.
///
/// \param params [in/out] List of parameters with which to configure
///   the solver.  If null, we configure the solver with default
///   parameters.
template<class Scalar, class MV, class OP>
Teuchos::RCP<SolverManager<Scalar, MV, OP> >
makeSolverManagerFromEnum (const EBelosSolverType solverType,
                           const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  typedef SolverManager<Scalar, MV, OP> base_type;

  switch (solverType) {
  case SOLVER_TYPE_BLOCK_GMRES: {
    typedef BlockGmresSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_PSEUDO_BLOCK_GMRES: {
    typedef PseudoBlockGmresSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_BLOCK_CG: {
    typedef BlockCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_PSEUDO_BLOCK_CG: {
    typedef PseudoBlockCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_GCRODR: {
    typedef GCRODRSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_RCG: {
    typedef RCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_MINRES: {
    typedef MinresSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_LSQR: {
    typedef LSQRSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_STOCHASTIC_CG: {
    typedef PseudoBlockStochasticCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_TFQMR: {
    typedef TFQMRSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_PSEUDO_BLOCK_TFQMR: {
    typedef PseudoBlockTFQMRSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_GMRES_POLY: {
    typedef GmresPolySolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_PCPG: {
    typedef PCPGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_FIXED_POINT: {
    typedef FixedPointSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  case SOLVER_TYPE_BICGSTAB: {
    typedef BiCGStabSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
  }
  default: // Fall through; let the code below handle it.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Belos::SolverFactory: Invalid EBelosSolverType "
      "enum value " << solverType << ".  Please report this bug to the Belos "
      "developers.");
  }

  // Compiler guard.  This may result in a warning on some compilers
  // for an unreachable statement, but it will prevent a warning on
  // other compilers for a "missing return statement at end of
  // non-void function."
  TEUCHOS_UNREACHABLE_RETURN(Teuchos::null);
}

template<class SolverManagerBaseType, class SolverManagerType>
Teuchos::RCP<SolverManagerBaseType>
makeSolverManagerTmpl (const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;

  RCP<SolverManagerType> solver = rcp (new SolverManagerType);

  // Some solvers may not like to get a null ParameterList.  If params
  // is null, replace it with an empty parameter list.  The solver
  // will fill in default parameters for that case.  Use the name of
  // the solver's default parameters to name the new empty list.
  RCP<ParameterList> pl;
  if (params.is_null()) {
    pl = parameterList (solver->getValidParameters ()->name ());
  } else {
    pl = params;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    pl.is_null(), std::logic_error,
    "Belos::SolverFactory: ParameterList to pass to solver is null.  This "
    "should never happen.  Please report this bug to the Belos developers.");
  solver->setParameters (pl);
  return solver;
}

} // namespace Details


template<class Scalar, class MV, class OP>
SolverFactory<Scalar, MV, OP>::SolverFactory()
{}

template<class Scalar, class MV, class OP>
Teuchos::RCP<typename SolverFactory<Scalar, MV, OP>::solver_base_type>
SolverFactory<Scalar, MV, OP>::
create (const std::string& solverName,
        const Teuchos::RCP<Teuchos::ParameterList>& solverParams)
{
  const char prefix[] = "Belos::SolverFactory: ";

  // Upper-case version of the input solver name.
  std::string solverNameUC (solverName);
  {
    typedef std::string::value_type char_t;
    typedef std::ctype<char_t> facet_type;
    const facet_type& facet = std::use_facet<facet_type> (std::locale ());

    const std::string::size_type len = solverName.size ();
    for (std::string::size_type k = 0; k < len; ++k) {
      solverNameUC[k] = facet.toupper (solverName[k]);
    }
  }

  // Check whether the given name is an alias.
  std::pair<std::string, bool> aliasResult =
    Details::getCanonicalNameFromAlias (solverNameUC);
  const std::string candidateCanonicalName = aliasResult.first;
  const bool isAnAlias = aliasResult.second;

  // Get the canonical name.
  const Details::EBelosSolverType solverEnum =
    Details::getEnumFromCanonicalName (isAnAlias ?
                                       candidateCanonicalName :
                                       solverNameUC);
  const bool validCanonicalName =
    (solverEnum != Details::SOLVER_TYPE_UPPER_BOUND);

  // Check whether we found a canonical name.  If we didn't and the
  // input name is a valid alias, that's a bug.  Otherwise, the input
  // name is invalid.
  TEUCHOS_TEST_FOR_EXCEPTION
    (! validCanonicalName && isAnAlias, std::logic_error,
     prefix << "Valid alias \"" << solverName << "\" has candidate canonical "
     "name \"" << candidateCanonicalName << "\", which is not a canonical "
     "solver name.  Please report this bug to the Belos developers.");
  TEUCHOS_TEST_FOR_EXCEPTION
    (! validCanonicalName && ! isAnAlias, std::invalid_argument,
     prefix << "Invalid solver name \"" << solverName << "\".");

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

  return Details::makeSolverManagerFromEnum<Scalar, MV, OP> (solverEnum, pl);
}


template<class Scalar, class MV, class OP>
std::string
SolverFactory<Scalar, MV, OP>::description() const
{
  using Teuchos::TypeNameTraits;

  std::ostringstream out;
  out << "\"Belos::SolverFactory\": {";
  if (this->getObjectLabel () != "") {
    out << "Label: " << this->getObjectLabel () << ", ";
  }
  out << "Scalar: " << TypeNameTraits<Scalar>::name ()
      << ", MV: " << TypeNameTraits<MV>::name ()
      << ", OP: " << TypeNameTraits<OP>::name ()
      << "}";
  return out.str ();
}


template<class Scalar, class MV, class OP>
void
SolverFactory<Scalar, MV, OP>::
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
    out << "Scalar: " << TypeNameTraits<Scalar>::name () << endl
        << "MV: " << TypeNameTraits<MV>::name () << endl
        << "OP: " << TypeNameTraits<OP>::name () << endl;
  }

  // At higher verbosity levels, print out the list of supported solvers.
  if (vl > Teuchos::VERB_LOW) {
    Teuchos::OSTab tab1 (out);
    out << "Number of solvers: " << numSupportedSolvers ()
        << endl;
    out << "Canonical solver names: ";
    printStringArray (out, Details::canonicalSolverNames ());
    out << endl;

    out << "Aliases to canonical names: ";
    printStringArray (out, Details::solverNameAliases ());
    out << endl;
  }
}

template<class Scalar, class MV, class OP>
int
SolverFactory<Scalar, MV, OP>::numSupportedSolvers () const
{
  return Details::numSupportedSolvers ();
}

template<class Scalar, class MV, class OP>
Teuchos::Array<std::string>
SolverFactory<Scalar, MV, OP>::supportedSolverNames () const
{
  typedef std::vector<std::string>::const_iterator iter_type;
  Teuchos::Array<std::string> names;

  {
    std::vector<std::string> aliases = Details::solverNameAliases ();
    for (iter_type iter = aliases.begin (); iter != aliases.end (); ++iter) {
      names.push_back (*iter);
    }
  }
  {
    std::vector<std::string> canonicalNames = Details::canonicalSolverNames ();
    for (iter_type iter = canonicalNames.begin (); iter != canonicalNames.end (); ++iter) {
      names.push_back (*iter);
    }
  }
  return names;
}

} // namespace Belos

#endif // __Belos_SolverFactory_hpp

