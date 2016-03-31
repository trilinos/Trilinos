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

#include <Teuchos_Array.hpp>
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

namespace details {

/// \enum EBelosSolverType
/// \brief 1-to-1 enumeration of all supported SolverManager subclasses.
/// \author Mark Hoemmen
///
/// This enum is an implementation detail of \c SolverFactory.
/// Users of \c SolverFactory should not refer to this enum or
/// rely on the symbols or integer values therein.  We declare
/// it here for later use by \c SolverFactory.
///
/// Belos developers who have implemented a new solver (i.e., a new
/// subclass of \c SolverManager) and who want to make the solver
/// available through the \c SolverFactory should first add a new enum
/// symbol corresponding to their solver to the end of the list.  They
/// should then follow the instructions provided in the \c
/// SolverFactory documentation.
///
/// \c SolverFactory was written to be independent of the actual enum
/// values, so Belos developers are allowed to rearrange the symbols.
enum EBelosSolverType {
  SOLVER_TYPE_BLOCK_GMRES,
  SOLVER_TYPE_PSEUDO_BLOCK_GMRES,
  SOLVER_TYPE_BLOCK_CG,
  SOLVER_TYPE_PSEUDO_BLOCK_CG,
  SOLVER_TYPE_GCRODR,
  SOLVER_TYPE_RCG,
  SOLVER_TYPE_MINRES,
  SOLVER_TYPE_LSQR,
  SOLVER_TYPE_STOCHASTIC_CG,
  SOLVER_TYPE_TFQMR,
  SOLVER_TYPE_PSEUDO_BLOCK_TFQMR,
  SOLVER_TYPE_GMRES_POLY,
  SOLVER_TYPE_PCPG,
  SOLVER_TYPE_FIXED_POINT,
  SOLVER_TYPE_BICGSTAB
};

} // namespace details

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
///      Solver name | Aliases | Solver Manager Class
///  ----------- | ------- | ----------
///  Pseudoblock GMRES |  GMRES, Pseudo Block GMRES, PseudoBlockGMRES, PseudoBlockGmres | \c PseudoBlockGmresSolMgr
///  Block GMRES | Flexible GMRES | \c BlockGmresSolMgr
///  Block CG |  |  \c BlockCGSolMgr
///  Pseudoblock CG | PseudoBlockCG, Pseudo Block CG | \c PseudoBlockCGSolMgr
///  Pseudoblock Stochastic CG | Stochastic CG | \c PseudoBlockStochasticCGSolMgr
///  GCRODR | Recycling GMRES | \c GCRODRSolMgr
///  RCG | Recycling CG | \c RCGSolMgr
///  MINRES | | \c MinresSolMgr
///  LSQR | | \c LSQRSolMgr
///  TFQMR | Transpose-Free QMR | \c TFQMRSolMgr
///  Pseudoblock TFQMR | Pseudo Block Transpose-Free QMR | \c PseudoBlockTFQMRSolMgr
///  Hybrid Block GMRES | GmresPoly, Seed GMRES | \c GmresPolySolMgr
///  PCPG | CGPoly, Seed CG | \c PCPGSolMgr
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
///      details::EBelosSolverType enum. </li>
/// <li> If necessary, specialize details::makeSolverManagerTmpl for
///      their SolverManager subclass.  In most cases, the default
///      implementation suffices. </li>
/// <li> Add a case for their enum symbol that instantiates their
///      solver to the long switch-case statement in
///      details::makeSolverManagerFromEnum. </li>
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
  /// \brief Map from solver name alias to canonical solver name.
  ///
  /// The keys of this map do not necessarily include canonical solver
  /// names.  If a candidate name isn't a key in this map, then it
  /// must be a canonical name in order to be valid.  There doesn't
  /// need to be an alias for each solver.
  ///
  /// \note To Belos developers: If you want to add a new alias, first
  ///   add the mapping from alias to canonical solver name in the
  ///   SolverFactory constructor.  Then, edit
  ///   reviseParameterListForAlias() to do any modifications of the
  ///   input ParameterList associated with that alias.
  std::map<std::string, std::string> aliasToCanonicalName_;

  /// \brief Map from canonical solver name to solver enum value.
  ///
  /// Access the keys to get the list of canonical solver names.
  ///
  /// \note To Belos developers: If you add a new solver, start with
  ///   the documentation of details::EBelosSolverType for
  ///   instructions.  Each new solver needs a canonical name (a
  ///   string), which is a key into this map.  The map from canonical
  ///   name to enum value is set up in the \c SolverFactory
  ///   constructor.  The details::makeSolverManagerFromEnum()
  ///   function in turn takes the enum value and parameter list, and
  ///   returns an instance of the appropriate subclass of
  ///   SolverManager.
  std::map<std::string, details::EBelosSolverType> canonicalNameToEnum_;

  /// \brief Modify the input ParameterList appropriately for the given alias.
  ///
  /// Some aliases include modifications or special checking of the
  /// input ParameterList.  All alias-related ParameterList revision
  /// happens in this method.
  void
  reviseParameterListForAlias (const std::string& aliasName,
                               Teuchos::ParameterList& solverParams);

  //! List of canonical solver names.
  Teuchos::Array<std::string> canonicalSolverNames () const;

  //! List of supported aliases (to canonical solver names).
  Teuchos::Array<std::string> solverNameAliases () const;

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
};


namespace details {

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
  return Teuchos::null;
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

} // namespace details


template<class Scalar, class MV, class OP>
SolverFactory<Scalar, MV, OP>::SolverFactory()
{
  aliasToCanonicalName_["GMRES"] = "PSEUDOBLOCK GMRES";
  // NOTE (mfh 29 Nov 2011) Accessing the flexible capability requires
  // setting a parameter in the solver's parameter list.  This affects
  // the SolverFactory's interface, since using the "Flexible GMRES"
  // alias requires modifying the user's parameter list if necessary.
  // This is a good idea because users may not know about the
  // parameter, or may have forgotten.
  //
  // NOTE (mfh 12 Aug 2015) The keys and values need to be all uppercase.
  aliasToCanonicalName_["BLOCK GMRES"] = "BLOCK GMRES";
  aliasToCanonicalName_["FLEXIBLE GMRES"] = "BLOCK GMRES";
  aliasToCanonicalName_["CG"] = "PSEUDOBLOCK CG";
  aliasToCanonicalName_["PSEUDOBLOCKCG"] = "PSEUDOBLOCK CG";
  aliasToCanonicalName_["STOCHASTIC CG"] = "PSEUDOBLOCK STOCHASTIC CG";
  aliasToCanonicalName_["RECYCLING CG"] = "RCG";
  aliasToCanonicalName_["RECYCLING GMRES"] = "GCRODR";
  // For compatibility with Stratimikos' Belos adapter.
  aliasToCanonicalName_["PSEUDO BLOCK GMRES"] = "PSEUDOBLOCK GMRES";
  aliasToCanonicalName_["PSEUDOBLOCKGMRES"] = "PSEUDOBLOCK GMRES";
  aliasToCanonicalName_["PSEUDO BLOCK CG"] = "PSEUDOBLOCK CG";
  aliasToCanonicalName_["PSEUDOBLOCKCG"] = "PSEUDOBLOCK CG";
  aliasToCanonicalName_["TRANSPOSE-FREE QMR"] = "TFQMR";
  aliasToCanonicalName_["PSEUDO BLOCK TFQMR"] = "PSEUDOBLOCK TFQMR";
  aliasToCanonicalName_["PSEUDO BLOCK TRANSPOSE-FREE QMR"] = "PSEUDOBLOCK TFQMR";
  aliasToCanonicalName_["GMRESPOLY"] = "HYBRID BLOCK GMRES";
  aliasToCanonicalName_["SEED GMRES"] = "HYBRID BLOCK GMRES";
  aliasToCanonicalName_["CGPOLY"] = "PCPG";
  aliasToCanonicalName_["SEED CG"] = "PCPG";
  aliasToCanonicalName_["FIXED POINT"] = "FIXED POINT";
  aliasToCanonicalName_["BICGSTAB"] = "BICGSTAB";

  // Mapping from canonical solver name (a string) to its
  // corresponding enum value.  This mapping is one-to-one.
  //
  // NOTE (mfh 12 Aug 2015) The keys need to be all uppercase.
  canonicalNameToEnum_["BLOCK GMRES"] = details::SOLVER_TYPE_BLOCK_GMRES;
  canonicalNameToEnum_["PSEUDOBLOCK GMRES"] = details::SOLVER_TYPE_PSEUDO_BLOCK_GMRES;
  canonicalNameToEnum_["BLOCK CG"] = details::SOLVER_TYPE_BLOCK_CG;
  canonicalNameToEnum_["PSEUDOBLOCK CG"] = details::SOLVER_TYPE_PSEUDO_BLOCK_CG;
  canonicalNameToEnum_["PSEUDOBLOCK STOCHASTIC CG"] = details::SOLVER_TYPE_STOCHASTIC_CG;
  canonicalNameToEnum_["GCRODR"] = details::SOLVER_TYPE_GCRODR;
  canonicalNameToEnum_["RCG"] = details::SOLVER_TYPE_RCG;
  canonicalNameToEnum_["MINRES"] = details::SOLVER_TYPE_MINRES;
  canonicalNameToEnum_["LSQR"] = details::SOLVER_TYPE_LSQR;
  canonicalNameToEnum_["TFQMR"] = details::SOLVER_TYPE_TFQMR;
  canonicalNameToEnum_["PSEUDOBLOCK TFQMR"] = details::SOLVER_TYPE_PSEUDO_BLOCK_TFQMR;
  canonicalNameToEnum_["HYBRID BLOCK GMRES"] = details::SOLVER_TYPE_GMRES_POLY;
  canonicalNameToEnum_["PCPG"] = details::SOLVER_TYPE_PCPG;
  canonicalNameToEnum_["FIXED POINT"] = details::SOLVER_TYPE_FIXED_POINT;
  canonicalNameToEnum_["BICGSTAB"] = details::SOLVER_TYPE_BICGSTAB;
}


template<class Scalar, class MV, class OP>
void
SolverFactory<Scalar, MV, OP>::
reviseParameterListForAlias (const std::string& aliasName,
                             Teuchos::ParameterList& solverParams)
{
  if (aliasName == "FLEXIBLE GMRES") {
    // "Gmres" uses title case in this solver's parameter list.  For
    // our alias, we prefer the all-capitals "GMRES" that the
    // algorithm's authors (Saad and Schultz) used.
    solverParams.set ("Flexible Gmres", true);
  }
}


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
  std::map<std::string, std::string>::const_iterator aliasIter =
    aliasToCanonicalName_.find (solverNameUC);
  const bool isAnAlias = (aliasIter != aliasToCanonicalName_.end());
  const std::string candidateCanonicalName =
    isAnAlias ? aliasIter->second : solverNameUC;

  // Get the canonical name.
  std::map<std::string, details::EBelosSolverType>::const_iterator canonicalIter =
    canonicalNameToEnum_.find (candidateCanonicalName);
  const bool validCanonicalName = (canonicalIter != canonicalNameToEnum_.end());

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
    reviseParameterListForAlias (solverNameUC, *pl);
  }

  return details::makeSolverManagerFromEnum<Scalar, MV, OP> (canonicalIter->second, pl);
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
    printStringArray (out, canonicalSolverNames ());
    out << endl;

    out << "Aliases to canonical names: ";
    printStringArray (out, solverNameAliases ());
    out << endl;
  }
}

template<class Scalar, class MV, class OP>
int
SolverFactory<Scalar, MV, OP>::numSupportedSolvers () const
{
  return static_cast<int> (canonicalNameToEnum_.size());
}

template<class Scalar, class MV, class OP>
Teuchos::Array<std::string>
SolverFactory<Scalar, MV, OP>::canonicalSolverNames () const
{
  Teuchos::Array<std::string> canonicalNames;
  typedef std::map<std::string, details::EBelosSolverType>::const_iterator iter_type;
  for (iter_type iter = canonicalNameToEnum_.begin();
       iter != canonicalNameToEnum_.end(); ++iter) {
    canonicalNames.push_back (iter->first);
  }
  return canonicalNames;
}

template<class Scalar, class MV, class OP>
Teuchos::Array<std::string>
SolverFactory<Scalar, MV, OP>::solverNameAliases () const
{
  Teuchos::Array<std::string> names;
  {
    typedef std::map<std::string, std::string>::const_iterator iter_type;
    for (iter_type iter = aliasToCanonicalName_.begin();
         iter != aliasToCanonicalName_.end(); ++iter) {
      names.push_back (iter->first);
    }
  }
  return names;
}

template<class Scalar, class MV, class OP>
Teuchos::Array<std::string>
SolverFactory<Scalar, MV, OP>::supportedSolverNames () const
{
  Teuchos::Array<std::string> names;
  {
    typedef std::map<std::string, std::string>::const_iterator iter_type;
    for (iter_type iter = aliasToCanonicalName_.begin();
         iter != aliasToCanonicalName_.end(); ++iter) {
      names.push_back (iter->first);
    }
  }
  {
    typedef std::map<std::string, details::EBelosSolverType>::const_iterator iter_type;
    for (iter_type iter = canonicalNameToEnum_.begin();
         iter != canonicalNameToEnum_.end(); ++iter) {
      names.push_back (iter->first);
    }
  }
  return names;
}

} // namespace Belos

#endif // __Belos_SolverFactory_hpp

