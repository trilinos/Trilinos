// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef IFPACK2_DETAILS_ONELEVELFACTORY_DEF_HPP
#define IFPACK2_DETAILS_ONELEVELFACTORY_DEF_HPP

#include "Ifpack2_Factory.hpp"
#include "Ifpack2_Utilities.hpp"
#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_Details_DenseSolver.hpp"
#include "Ifpack2_Diagonal.hpp"
#include "Ifpack2_IdentitySolver.hpp"
#include "Ifpack2_ILUT.hpp"
#include "Ifpack2_MDF.hpp"
#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_Experimental_RBILUK.hpp"
#include "Ifpack2_BlockRelaxation.hpp"
#include "Ifpack2_BandedContainer.hpp"
#include "Ifpack2_DenseContainer.hpp"
#include "Ifpack2_DatabaseSchwarz.hpp"
#include "Ifpack2_SparseContainer.hpp"
#include "Ifpack2_TriDiContainer.hpp"
#include "Ifpack2_LocalSparseTriangularSolver.hpp"
#include "Ifpack2_Hiptmair.hpp"

#ifdef HAVE_IFPACK2_SHYLU_NODEFASTILU
#include "Ifpack2_Details_Fic.hpp"
#include "Ifpack2_Details_Fildl.hpp"
#include "Ifpack2_Details_Filu.hpp"
#endif // HAVE_IFPACK2_SHYLU_NODEFASTILU

#ifdef HAVE_IFPACK2_AMESOS2
#  include "Ifpack2_Details_Amesos2Wrapper.hpp"
#endif // HAVE_IFPACK2_AMESOS2

#ifdef HAVE_IFPACK2_HYPRE
#  include "Ifpack2_Hypre.hpp"
#endif // HAVE_IFPACK2_HYPRE


namespace Ifpack2 {
namespace Details {

template<class MatrixType>
Teuchos::RCP<typename OneLevelFactory<MatrixType>::prec_type>
OneLevelFactory<MatrixType>::create (const std::string& precType,
                                     const Teuchos::RCP<const row_matrix_type>& matrix) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  RCP<prec_type> prec;

  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper = canonicalize(precType);

  if (precTypeUpper == "CHEBYSHEV") {
    // We have to distinguish Ifpack2::Chebyshev from its
    // implementation class Ifpack2::Details::Chebyshev.
    prec = rcp (new ::Ifpack2::Chebyshev<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "DENSE" || precTypeUpper == "LAPACK") {
    prec = rcp (new Details::DenseSolver<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "AMESOS2") {
#ifdef HAVE_IFPACK2_AMESOS2
    prec = rcp (new Details::Amesos2Wrapper<row_matrix_type> (matrix));
#else
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Details::OneLevelFactory: "
      "You may not ask for the preconditioner \"AMESOS2\" unless "
      "you have built Trilinos with the Amesos2 package enabled.");
#endif // HAVE_IFPACK2_AMESOS2
  }
  else if (precTypeUpper == "DIAGONAL") {
    prec = rcp (new Diagonal<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "ILUT") {
    prec = rcp (new ILUT<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "RELAXATION") {
    prec = rcp (new Relaxation<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "RILUK") {
    prec = rcp (new RILUK<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "MDF") {
    prec = rcp (new MDF<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "RBILUK") {
    prec = rcp (new Experimental::RBILUK<row_matrix_type>(matrix));
  }
  else if (precTypeUpper == "FAST_IC" || precTypeUpper == "FAST_ILU" || precTypeUpper == "FAST_ILU_B" || precTypeUpper == "FAST_ILDL") {
    #ifdef HAVE_IFPACK2_SHYLU_NODEFASTILU
    {
      if(precTypeUpper == "FAST_IC")
        prec = rcp (new Details::Fic<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(matrix));
      else if(precTypeUpper == "FAST_ILU")
        prec = rcp (new Details::Filu<scalar_type, local_ordinal_type, global_ordinal_type, node_type, false>(matrix));
      else if(precTypeUpper == "FAST_ILU_B")
        prec = rcp (new Details::Filu<scalar_type, local_ordinal_type, global_ordinal_type, node_type, true>(matrix));
      else if(precTypeUpper == "FAST_ILDL")
        prec = rcp (new Details::Fildl<scalar_type, local_ordinal_type, global_ordinal_type, node_type>(matrix));
    }
    #else
    {
      throw std::invalid_argument("The Ifpack2 FastIC, FastILU and FastILDL preconditioners require the FastILU subpackage of ShyLU to be enabled\n"
                                  "To enable FastILU, set the CMake option Trilinos_ENABLE_ShyLU_NodeFastILU=ON");
    }
    #endif
  }
  else if (precTypeUpper == "KRYLOV") {
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "The \"KRYLOV\" preconditioner option has "
       "been deprecated and removed.  If you want a Krylov solver, use the "
       "Belos package.");
  }
  else if (precTypeUpper == "BLOCK_RELAXATION" ||
           precTypeUpper == "BLOCK RELAXATION" ||
           precTypeUpper == "BLOCKRELAXATION"  ||
           precTypeUpper == "DENSE_BLOCK_RELAXATION" ||
           precTypeUpper == "DENSE BLOCK RELAXATION" ||
           precTypeUpper == "DENSEBLOCKRELAXATION" ) {
    // NOTE (mfh 12 Aug 2016) Choice of "container type" is now a
    // run-time parameter.  The "ContainerType" template parameter is
    // now always Container<row_matrix_type>.
    prec = rcp (new BlockRelaxation<row_matrix_type> (matrix));
    Teuchos::ParameterList params;
    params.set ("relaxation: container", "Dense");
    prec->setParameters (params);
  }
  else if (precTypeUpper == "DATABASE SCHWARZ") {
    prec = rcp (new DatabaseSchwarz<row_matrix_type> (matrix));
  }
  else if (precTypeUpper == "SPARSE_BLOCK_RELAXATION" ||
           precTypeUpper == "SPARSE BLOCK RELAXATION" ||
           precTypeUpper == "SPARSEBLOCKRELAXATION" ) {
    // FIXME (mfh 22 May 2014) We would prefer to have the choice of
    // dense or sparse blocks (the "container type") be a run-time
    // decision.  This will require refactoring BlockRelaxation so
    // that the "container type" is not a template parameter.  For
    // now, we default to use dense blocks.
    //typedef SparseContainer<row_matrix_type, ILUT<row_matrix_type>> container_type;
#ifdef HAVE_IFPACK2_AMESOS2
    prec = rcp (new BlockRelaxation<row_matrix_type> (matrix));
    Teuchos::ParameterList params;
    params.set ("relaxation: container", "SparseAmesos2");
    prec->setParameters (params);
#else
    TEUCHOS_TEST_FOR_EXCEPTION
      (true, std::invalid_argument, "Ifpack2::Details::OneLevelFactory: "
      "\"SPARSE BLOCK RELAXATION\" requires building Trilinos with Amesos2 enabled.");
#endif
  }
  else if (precTypeUpper == "TRIDI_RELAXATION" ||
           precTypeUpper == "TRIDI RELAXATION" ||
           precTypeUpper == "TRIDIRELAXATION" ||
           precTypeUpper == "TRIDIAGONAL_RELAXATION" ||
           precTypeUpper == "TRIDIAGONAL RELAXATION" ||
           precTypeUpper == "TRIDIAGONALRELAXATION") {
    prec = rcp (new BlockRelaxation<row_matrix_type> (matrix));
    Teuchos::ParameterList params;
    params.set ("relaxation: container", "TriDi");
    prec->setParameters (params);
  }
  else if (precTypeUpper == "BANDED_RELAXATION" ||
           precTypeUpper == "BANDED RELAXATION" ||
           precTypeUpper == "BANDEDRELAXATION") {
    prec = rcp (new BlockRelaxation<row_matrix_type> (matrix));
    Teuchos::ParameterList params;
    params.set ("relaxation: container", "Banded");
    prec->setParameters (params);
  }
  else if (precTypeUpper == "IDENTITY" || precTypeUpper == "IDENTITY_SOLVER") {
    prec = rcp (new IdentitySolver<row_matrix_type> (matrix));
  }

  else if (precTypeUpper == "LOCAL SPARSE TRIANGULAR SOLVER" ||
           precTypeUpper == "LOCAL_SPARSE_TRIANGULAR_SOLVER" ||
           precTypeUpper == "LOCALSPARSETRIANGULARSOLVER" ||
           precTypeUpper == "SPARSE TRIANGULAR SOLVER" ||
           precTypeUpper == "SPARSE_TRIANGULAR_SOLVER" ||
           precTypeUpper == "SPARSETRIANGULARSOLVER") {
    prec = rcp (new LocalSparseTriangularSolver<row_matrix_type> (matrix));
  }
  else if(precTypeUpper == "HIPTMAIR") {
    prec = rcp (new Hiptmair<row_matrix_type> (matrix));
  }
#ifdef HAVE_IFPACK2_HYPRE
  else if (precTypeUpper == "HYPRE") {
    prec = rcp (new Hypre<row_matrix_type> (matrix));
  }
#endif
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Details::OneLevelFactory::create: "
      "Invalid preconditioner type \"" << precType << "\".");
  }

  TEUCHOS_TEST_FOR_EXCEPTION(
    prec.is_null (), std::logic_error, "Ifpack2::Details::OneLevelFactory::"
    "create: Return value is null right before return.  This should never "
    "happen.  Please report this bug to the Ifpack2 developers.");
  return prec;
}

template<class MatrixType>
bool
OneLevelFactory<MatrixType>::isSupported (const std::string& precType) const
{
  // precTypeUpper is the upper-case version of precType.
  std::string precTypeUpper = canonicalize(precType);
  std::vector<std::string> supportedNames = {
    "CHEBYSHEV", "DENSE", "LAPACK",
#ifdef HAVE_IFPACK2_AMESOS2
    "AMESOS2",
#endif
    "DIAGONAL", "ILUT", "RELAXATION", "RILUK", "RBILUK", "MDF",
#ifdef HAVE_IFPACK2_SHYLU_NODEFASTILU
    "FAST_IC", "FAST_ILU", "FAST_ILU_B", "FAST_ILDL",
#endif
    "BLOCK_RELAXATION", "BLOCK RELAXATION", "BLOCKRELAXATION", "DENSE_BLOCK_RELAXATION", "DENSE BLOCK RELAXATION", "DENSEBLOCKRELAXATION",
    "DATABASE SCHWARZ",
#ifdef HAVE_IFPACK2_AMESOS2
    "SPARSE_BLOCK_RELAXATION", "SPARSE BLOCK RELAXATION", "SPARSEBLOCKRELAXATION",
#endif
    #ifdef HAVE_IFPACK2_HYPRE
    "HYPRE",
#endif
    "TRIDI_RELAXATION", "TRIDI RELAXATION", "TRIDIRELAXATION", "TRIDIAGONAL_RELAXATION", "TRIDIAGONAL RELAXATION", "TRIDIAGONALRELAXATION",
    "BANDED_RELAXATION", "BANDED RELAXATION", "BANDEDRELAXATION",
    "IDENTITY", "IDENTITY_SOLVER",
    "LOCAL SPARSE TRIANGULAR SOLVER", "LOCAL_SPARSE_TRIANGULAR_SOLVER", "LOCALSPARSETRIANGULARSOLVER", "SPARSE TRIANGULAR SOLVER", "SPARSE_TRIANGULAR_SOLVER", "SPARSETRIANGULARSOLVER",
    "HIPTMAIR"
  };
  // const size_t numSupportedNames = supportedNames.size();
  // const auto end = supportedNames + numSupportedNames;
  auto it = std::find(std::begin(supportedNames), std::end(supportedNames), precTypeUpper);
  return it != std::end(supportedNames);
}

} // namespace Details
} // namespace Ifpack2

#define IFPACK2_DETAILS_ONELEVELFACTORY_INSTANT(S,LO,GO,N)              \
  template class Ifpack2::Details::OneLevelFactory< Tpetra::RowMatrix<S, LO, GO, N> >;

#endif // IFPACK2_DETAILS_ONELEVELFACTORY_DEF_HPP
