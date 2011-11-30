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

#include <Teuchos_Array.hpp>
#include <Teuchos_Describable.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_TypeNameTraits.hpp>

#include <algorithm>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Belos {

namespace details {

/// \enum EBelosSolverType
/// \brief 1-to-1 enumeration of all supported SolverManager subclasses.
/// \author Mark Hoemmen
enum EBelosSolverType {
  SOLVER_TYPE_BLOCK_GMRES,
  SOLVER_TYPE_PSEUDO_BLOCK_GMRES,
  SOLVER_TYPE_BLOCK_CG,
  SOLVER_TYPE_PSEUDO_BLOCK_CG,
  SOLVER_TYPE_GCRODR,
  SOLVER_TYPE_RCG,
  SOLVER_TYPE_MINRES,
  SOLVER_TYPE_LSQR
};

} // namespace details

/// \class SolverFactory
/// \brief Factory for all solvers which Belos supports.
/// \author Mark Hoemmen
///
/// Belos implements several different iterative solvers.  The usual
/// way in which users interact with these solvers is through
/// appropriately named subclasses of \c SolverManager.  This
/// factory class tells users which solvers are supported.  It can
/// initialize and return any supported subclass of \c
/// SolverManager, given a short name of the subclass (such as
/// "GMRES" or "CG").
///
/// This class' template parameters are the same as those of \c
/// SolverManager: Scalar is the scalar type (of entries in the
/// multivector), MV is the multivector type, and OP is the operator
/// type.  For example: Scalar=double, MV=Epetra_MultiVector, and
/// OP=Epetra_Operator will access the Epetra specialization of the
/// Belos solvers.
///
/// This factory implements \c Teuchos::Describable.  At higher
/// verbosity levels, the describe() method will print out the list of
/// names of supported solvers.  You can also get this list directly.
template<class Scalar, class MV, class OP>
class SolverFactory : public Teuchos::Describable {
public:
  /// \typedef solver_base_type
  /// \brief The type returned by \c makeSolver().
  typedef SolverManager<Scalar, MV, OP> solver_base_type;

  //! Constructor.
  SolverFactory();

  /// \brief Create, configure, and return the specified solver.
  ///
  /// \param solverName [in] Name of the solver.
  ///
  /// \param params [in/out] List of parameters with which to configure
  ///   the solver.  If null, we configure the solver with default
  ///   parameters.
  ///
  /// It is better to provide a non-null but empty parameter list,
  /// since in that case, the solver will fill in your list with
  /// parameters and their default values.  You can then inspect the
  /// parameter names and learn how to modify their default values.
  Teuchos::RCP<solver_base_type>
  makeSolver (const std::string& solverName, 
	      const Teuchos::RCP<Teuchos::ParameterList>& params);

  /// \brief Number of supported solvers.
  /// 
  /// This may differ from the number of supported solver
  /// <i>names</i>, since we may accept multiple names ("aliases") for
  /// commonly used solvers.
  int numSupportedSolvers () const;

  /// \brief List of supported solver names.
  ///
  /// The length of this list may differ from the number of supported
  /// solvers, since we may accept multiple names ("aliases") for
  /// commonly used solvers.
  Teuchos::Array<std::string> supportedSolverNames () const;

  //! Whether the given solver name names a supported solver.
  bool isSupported (const std::string& solverName) const;

  //! @name Implementation of Teuchos::Describable interface
  //@{
  std::string description() const;

  void describe (Teuchos::FancyOStream& out, 
		 const EVerbosityLevel verbLevel = Teuchos::Describable::verbLevel_default) const;
  //@}

private:
  /// \brief Map from solver name alias to canonical solver name.
  ///
  /// A canonical solver name need not be an alias.  If a candidate
  /// name isn't a key in this map, then it must be a canonical name
  /// in order to be valid.
  std::map<std::string, std::string> aliasToCanonicalName_;

  /// \brief Map from canonical solver name to solver enum value.
  ///
  /// Access the keys to get the list of canonical solver names.
  std::map<std::string, EBelosSolverType> canonicalNameToEnum_;

  //! List of canonical solver names.
  Teuchos::Array<std::string> canonicalSolverNames () const;

  //! List of supported aliases (to canonical solver names).
  Teuchos::Array<std::string> solverNameAliases () const;
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
    break;
  } 
  case SOLVER_TYPE_PSEUDO_BLOCK_GMRES: {
    typedef PseudoBlockGmresSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  case SOLVER_TYPE_BLOCK_CG: {
    typedef BlockCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  case SOLVER_TYPE_PSEUDO_BLOCK_CG: {
    typedef PseudoBlockCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  case SOLVER_TYPE_GCRODR: {
    typedef GCRODRSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  case SOLVER_TYPE_RCG: {
    typedef RCGSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  case SOLVER_TYPE_MINRES: {
    typedef MinresSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  case SOLVER_TYPE_LSQR: {
    typedef LSQRSolMgr<Scalar, MV, OP> impl_type;
    return makeSolverManagerTmpl<base_type, impl_type> (params);
    break;
  }
  default:
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error, 
			       "Invalid EBelosSolverType enum value " << solverType 
			       << ".  Please report this bug to the Belos developers.");
    // Compiler guard.
    return Teuchos::null;
  }
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
  // is null, replace it with the solver's default parameters.
  RCP<ParameterList> pl;
  if (params.is_null()) {
    pl = parameterList (*solver->getValidParameters());
  } else {
    pl = params;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(pl.is_null(), std::logic_error, 
			     "ParameterList to pass to solver is null.  This "
			     "should never happen.  Please report this bug to "
			     "the Belos developers.");
  solver->setParameters (pl);
  return solver;
}

} // namespace details


template<class Scalar, class MV, class OP>
SolverFactory<Scalar, MV, OP>::SolverFactory()
{
  // We don't need to add the canonical names to the list of aliases.
  aliasToCanonicalName_["GMRES"] = "Pseudoblock GMRES";
  // NOTE (mfh 29 Nov 2011) Accessing the flexible capability requires
  // setting a parameter in the solver's parameter list.  This affects
  // the SolverFactory's interface, since using the "Flexible GMRES"
  // alias requires modifying the user's parameter list if necessary.
  // This is a good idea because users may not know about the
  // parameter, or may have forgotten.
  aliasToCanonicalName_["Flexible GMRES"] = "Block GMRES";
  aliasToCanonicalName_["CG"] = "Pseudoblock CG";
  aliasToCanonicalName_["Recycling CG"] = "RCG";
  aliasToCanonicalName_["Recycling GMRES"] = "GCRODR";

  canonicalNameToEnum_["Block GMRES"] = SOLVER_TYPE_BLOCK_GMRES;
  canonicalNameToEnum_["Pseudoblock GMRES"] = SOLVER_TYPE_BLOCK_GMRES;
  canonicalNameToEnum_["Block CG"] = SOLVER_TYPE_BLOCK_CG;
  canonicalNameToEnum_["Pseudoblock CG"] = SOLVER_TYPE_PSEUDO_BLOCK_CG;
  canonicalNameToEnum_["GCRODR"] = SOLVER_TYPE_GCRODR;
  canonicalNameToEnum_["RCG"] = SOLVER_TYPE_RCG;
  canonicalNameToEnum_["MINRES"] = SOLVER_TYPE_MINRES;
  canonicalNameToEnum_["LSQR"] = SOLVER_TYPE_LSQR;
}

template<class Scalar, class MV, class OP>
Teuchos::RCP<typename SolverFactory<Scalar, MV, OP>::solver_base_type>
makeSolver (const std::string& solverName, 
	    const Teuchos::RCP<Teuchos::ParameterList>& params)
{
  // Check whether the given name is an alias.
  std::map<std::string, std::string>::const_iterator aliasIter = 
    aliasToCanonicalName_.find (solverName);
  const bool isAnAlias = (aliasIter != aliasToCanonicalName_.end());
  const std::string candidateCanonicalName = 
    isAnAlias ? aliasIter->second : solverName;

  // Get the canonical name.
  std::map<std::string, EBelosSolverType>::const_iterator canonicalIter =
    canonicalNameToEnum_.find (candidateCanonicalName);
  const bool validCanonicalName = (canonicalIter == canonicalNameToEnum_.end());

  // Check whether we found a canonical name.  If we didn't and the
  // input name is a valid alias, that's a bug.  Otherwise, the input
  // name is invalid.
  TEUCHOS_TEST_FOR_EXCEPTION(! validCanonicalName && isAnAlias, std::logic_error,
    "Valid alias \"" << solverName << "\" has candidate canonical name \"" 
    << candidateCanonicalName << "\", which is not a canonical solver name.  "
    "Please report this bug to the Belos developers.");
  TEUCHOS_TEST_FOR_EXCEPTION(! validCanonicalName && ! isAnAlias, 
    std::invalid_argument, "Invalid solver name \"" << solverName << "\".");
  //
  // Special case for certain alias(es) requires possibly modifying
  // the input parameter list.  If the input list is null, we create a
  // new list and use that.  This is OK because the effect of a null
  // parameter list input is to use default parameter values.  Thus,
  // we can always replace a null list with an empty list.
  //
  Teuchos::RCP<Teuchos::ParameterList> pl = 
    params.is_null() ? Teuchos::parameterList() : params;
  if (solverName == "Flexible GMRES") {
    // "Gmres" uses title case in this solver's parameter list.
    pl->set ("Flexible Gmres", true);
  }
  return makeSolverManagerFromEnum (canonicalIter->second, pl);
}


template<class Scalar, class MV, class OP>
std::string
SolverFactory<Scalar, MV, OP>::description() const
{
  using Teuchos::TypeNameTraits;
  std::ostringstream os;
  os << "Belos::SolverFactory<" << TypeNameTraits<Scalar>::name()
     << ", " << TypeNameTraits<Scalar>::name()
     << ", " << TypeNameTraits<MV>::name()
     << ", " << TypeNameTraits<OP>::name()
     << ">";
  return os.str();
}

template<class Scalar, class MV, class OP>
void
SolverFactory<Scalar, MV, OP>::
describe (Teuchos::FancyOStream& out, 
	  const EVerbosityLevel verbLevel) const
{
  using std::endl;
  typedef Teuchos::Array<std::string>::const_iterator iter_type;

  Teuchos::OSTab tab1 (out);
  out << this->description();

  // At higher verbosity levels, print out the list of supported solvers.
  if (static_cast<int> (verbLevel) > static_cast<int> (Teuchos::VERB_LOW)) {
    out << ":" << endl;
    Teuchos::OSTab tab2 (out);
    out << "Number of supported solvers: " << numSupportedSolvers() 
	<< endl;
    out << "Supported canonical solver names:";
    {
      Teuchos::Array<std::string> names = canonicalSolverNames();
      for (iter_type iter = names.begin(); iter != names.end(); ++iter) {
	out << *iter;
	if (iter + 1 != names.end()) {
	  out << ", ";
	}
      }
    }
    out << "Supported aliases to canonical names:";
    {
      Teuchos::Array<std::string> names = solverNameAliases();
      for (iter_type iter = names.begin(); iter != names.end(); ++iter) {
	out << *iter;
	if (iter + 1 != names.end()) {
	  out << ", ";
	}
      }
    }
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
  typedef std::map<std::string, EBelosSolverType>::const_iterator iter_type;
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
    typedef std::map<std::string, EBelosSolverType>::const_iterator iter_type;
    for (iter_type iter = canonicalNameToEnum_.begin(); 
	 iter != canonicalNameToEnum_.end(); ++iter) {
      names.push_back (iter->first);
    }
  }
  return names;
}

} // namespace Belos

#endif // __Belos_SolverFactory_hpp

