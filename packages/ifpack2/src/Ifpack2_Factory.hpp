/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
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
//@HEADER
*/

#ifndef IFPACK2_FACTORY_HPP
#define IFPACK2_FACTORY_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Relaxation.hpp"
#include "Ifpack2_Diagonal.hpp"
#include "Ifpack2_Chebyshev.hpp"
#include "Ifpack2_RILUK.hpp"
#include "Ifpack2_ILUT.hpp"
#include "Ifpack2_Krylov.hpp"
#include "Ifpack2_AdditiveSchwarz.hpp"
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
#include "Ifpack2_SupportGraph.hpp"
#endif
#include <locale>


namespace Ifpack2 {

//! \c true if the specified preconditioner type supports nonsymmetric matrices, else false.
bool supportsUnsymmetric (const std::string& prec_type);

//! A factory class to create Ifpack2 preconditioners.
/*!
Ifpack2::Factory contains just one method: create().
Using Ifpack2::Factory::create(), users can easily create a variety of
Ifpack2 preconditioners. 

The create() method has three arguments:
- a string, indicating the type of preconditioner to compute;
- a pointer to a Tpetra::RowMatrix, representing the matrix
  to be used to define the preconditioner;
- an optional integer (defaults to 0), that specifies the amount of
  overlap between processes (if the input matrix is distributed
  over multiple processes).

The first argument can assume the following values:
- "DIAGONAL": returns an instance of Ifpack2::Diagonal.
- "RELAXATION": returns an instance of Ifpack2::Relaxation.
- "CHEBYSHEV": returns an instance of Ifpack2::Chebyshev (overlap is ignored).
- "ILUT": returns an instance of Ifpack2::ILUT.
- "RILUK": returns an instance of Ifpack2::RILUK.
- "KRYLOV": returns an instance of Ifpack2::Krylov.

<P> The following fragment of code shows the
basic usage of this class.
\code
#include "Ifpack2_Factory.hpp"

...
typedef double Scalar;
typedef int    LocalOrdinal;
typedef int    GlobalOrdinal;
typedef Tpetra::DefaultPlatform::DefaultPlatformType::NodeType Node;
typedef Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> TCrsMatrix;
typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> TPrecond;
...
Ifpack2::Factory factory;

Teuchos::RCP<TCrsMatrix> A; // A is fillComplete()'d.
std::string PrecType = "ILUT"; // use incomplete LU on each process
Teuchos::RCP<TPrecond> prec = factory.create(PrecType, A);

Teuchos::ParameterList params;
params.set("fact: ilut level-of-fill", 5.0); // use ILUT(fill=5, drop=0)

prec->setParameters(List);
prec->initialize();
prec->compute();

// now prec can be used as a preconditioner

\endcode

*/
class Factory {
public:
  /** \brief Create an instance of Ifpack2_Preconditioner given the string
   * name of the preconditioner type.
   *
   * \param PrecType [in] Name of preconditioner type to be created. 
   * \param Matrix [in] Matrix used to define the preconditioner
   * \param overlap (in) Specified overlap; defaults to 0.
   *
   * Throw an exception if the preconditioner with that input name
   * does not exist.  Otherwise, return a newly created preconditioner
   * object.
   */
  template<class MatrixType>
  static
  Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                       typename MatrixType::local_ordinal_type,
                                       typename MatrixType::global_ordinal_type,
                                       typename MatrixType::node_type> >
  create(const std::string& prec_type,
         const Teuchos::RCP<const MatrixType>& matrix,
         const int overlap = 0);


  //! Clones a preconditioner for a different node type from an Ifpack2 RILUK or Chebyshev preconditioner
  template<class MatrixType, class M2>
  static
  Teuchos::RCP<Ifpack2::Preconditioner<typename M2::scalar_type,
                                       typename M2::local_ordinal_type,
                                       typename M2::global_ordinal_type,
                                       typename M2::node_type> >
  clone(const Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                       typename MatrixType::local_ordinal_type,
                                       typename MatrixType::global_ordinal_type,
                                       typename MatrixType::node_type> >& prec,
		const Teuchos::RCP<const M2>& matrix, const Teuchos::ParameterList& params = Teuchos::ParameterList());

};

/////////////////////////////////////////////
/////////////////////////////////////////////

template<class MatrixType>
Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> >
Factory::create(const std::string& prec_type,
                const Teuchos::RCP<const MatrixType>& matrix,
                const int overlap)
{
  using Teuchos::RCP;
  using Teuchos::rcp;
  typedef typename MatrixType::scalar_type Scalar;
  typedef typename MatrixType::local_ordinal_type LocalOrdinal;
  typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
  typedef typename MatrixType::node_type Node;
  RCP<Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec;

  // precTypeUpper is the upper-case version of prec_type.
  std::string precTypeUpper (prec_type);
  if (precTypeUpper.size () > 0) {
    std::locale locale;
    for (size_t k = 0; k < precTypeUpper.size (); ++k) {
      precTypeUpper[k] = std::toupper<char> (precTypeUpper[k], locale);
    }
  }

  const bool one_mpi_rank = (matrix->getComm ()->getSize () == 1);

  if (precTypeUpper == "ILUT") {
    // Note: ILUT doesn't work for multiple MPI ranks... you have to use AdditiveSchwarz.
    if (one_mpi_rank) {
      prec = rcp (new Ifpack2::ILUT<MatrixType> (matrix));
    }
    else {
      prec = rcp (new Ifpack2::AdditiveSchwarz<MatrixType, Ifpack2::ILUT<MatrixType> > (matrix, overlap));
    }
  }
  else if (precTypeUpper == "RILUK") {
    prec = rcp (new Ifpack2::RILUK<MatrixType> (matrix));
  }
  else if (precTypeUpper == "RELAXATION") {
    prec = rcp (new Ifpack2::Relaxation<MatrixType> (matrix));
  }
  else if (precTypeUpper == "CHEBYSHEV") {
    prec = rcp (new Ifpack2::Chebyshev<MatrixType> (matrix));
  }
  else if (precTypeUpper == "DIAGONAL") {
    prec = rcp (new Ifpack2::Diagonal<MatrixType> (matrix));
  }
  else if (precTypeUpper == "SCHWARZ") {
    prec = rcp (new Ifpack2::AdditiveSchwarz<MatrixType,Ifpack2::ILUT<MatrixType> > (matrix, overlap));
  }
  else if (precTypeUpper == "KRYLOV") {
    prec = rcp (new Ifpack2::Krylov<MatrixType, Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node> > (matrix));
  }
#if defined(HAVE_IFPACK2_EXPERIMENTAL) && defined(HAVE_IFPACK2_SUPPORTGRAPH)
  else if (precTypeUpper == "SUPPORTGRAPH") {
    if (one_mpi_rank) {
      prec = rcp (new Ifpack2::SupportGraph<MatrixType> (matrix));
    }
    else {
      prec = rcp (new Ifpack2::AdditiveSchwarz<MatrixType, Ifpack2::SupportGraph<MatrixType> > (matrix, overlap));
    }
  }
#endif
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "Ifpack2::Factory::create: "
      "Invalid preconditioner type \"" << prec_type << "\".");
  }
  return prec;
}

template<class MatrixType, class M2>
Teuchos::RCP<Ifpack2::Preconditioner<typename M2::scalar_type, typename M2::local_ordinal_type,typename M2::global_ordinal_type,typename M2::node_type> >
Factory::clone(const Teuchos::RCP<Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                typename MatrixType::local_ordinal_type,
                typename MatrixType::global_ordinal_type,
                typename MatrixType::node_type> >& prec,
                const Teuchos::RCP<const M2>& matrix, const Teuchos::ParameterList& params) {
  typedef typename M2::scalar_type scalar_type;
  typedef typename M2::local_ordinal_type local_ordinal_type;
  typedef typename M2::global_ordinal_type global_ordinal_type;
  typedef typename M2::node_type new_node_type;

  Teuchos::RCP<Ifpack2::Preconditioner<scalar_type, local_ordinal_type,global_ordinal_type, new_node_type> > new_prec;
  Teuchos::RCP<Ifpack2::Chebyshev<MatrixType> > chebyPrec;
  chebyPrec = Teuchos::rcp_dynamic_cast<Ifpack2::Chebyshev<MatrixType> >(prec);
  if (chebyPrec != Teuchos::null){
	new_prec = chebyPrec->clone(matrix, params);
	return new_prec;
  }
  Teuchos::RCP<Ifpack2::RILUK<MatrixType> > luPrec;
  luPrec = Teuchos::rcp_dynamic_cast<Ifpack2::RILUK<MatrixType> >(prec);
  if (luPrec != Teuchos::null){	
	new_prec = luPrec->clone(matrix);
	return new_prec;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    true, std::invalid_argument, "Ifpack2::Factory::clone: "
    "Invalid preconditioner type to clone \"");	
}

} //namespace Ifpack2

#endif
