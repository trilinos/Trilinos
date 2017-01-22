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

#ifndef IFPACK2_DETAILS_ONELEVELFACTORY_DECL_HPP
#define IFPACK2_DETAILS_ONELEVELFACTORY_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Tpetra_RowMatrix.hpp"
#include <type_traits> // std::is_same

namespace Ifpack2 {
namespace Details {

/*!
\class OneLevelFactory
\brief "Factory" for creating single-level preconditioners.

\warning This class is an implementation detail of Ifpack2.  Users
  must not rely on it continuing to exist, on its name, or on its
  contents.

This class exists so that Ifpack2's "multi-level" (nested)
preconditioners -- AdditiveSchwarz and SupportGraph (if enabled) --
can create a default inner preconditioner.  These preconditioners
can't use Ifpack2::Factory, because that would introduce a circular
dependency between them and Ifpack2::Factory.  (Ifpack2::Factory has
to be able to create AdditiveSchwarz, for example.)  Ifpack2 indicates
that a preconditioner is multi-level if it implements
Details::NestedPreconditioner.

We resolve the circular dependency using this class, OneLevelFactory.
OneLevelFactory knows how to create all of Ifpack2's "single-level"
(not multi-level) preconditioners.  This includes Diagonal,
Relaxation, and the incomplete factorizations (like ILUT and RILUK).

The preferred way to create Ifpack2 preconditioners is to use Factory.
Users who create a multi-level preconditioner using Factory will be
able to specify any arbitrary Ifpack2 preconditioner as an inner
preconditioner for the "outer" multi-level preconditioner.  However,
if users create a multi-level preconditioner directly (by calling its
constructor, not by using Factory), the preconditioner will use
OneLevelFactory to create its inner solver.  This means that if users
create a multi-level preconditioner by calling its constructor, then
they are limited to single-level inner preconditioners.  Full
generality is only possible using Factory, or by creating an arbitrary
inner preconditioner themselves and giving it to the "outer"
preconditioner (by calling its setInnerPreconditioner() method).

This class' create() method lets users create an instance of any
single-level Ifpack2 preconditioner.

The create() method has two arguments:
<ol>
<li> a string, indicating the type of preconditioner to compute; and </li>
<li> a Tpetra::RowMatrix, representing the input matrix (to be used to
     define the preconditioner) </li>
</ol>

The first argument is not case sensitive.  It can assume the following
values:
<ul>
<li> "AMESOS2": returns an instance of Details::Amesos2Wrapper (if the
     Amesos2 package is enabled) </li>
<li> "BANDED RELAXATION": returns an instance of BlockRelaxation with banded matrix blocks </li>
<li> "CHEBYSHEV": returns an instance of Chebyshev </li>
<li> "DENSE" or "LAPACK": returns an instance of Details::DenseSolver </li>
<li> "DENSE BLOCK RELAXATION": returns an instance of BlockRelaxation with dense blocks </li>
<li> "DIAGONAL": returns an instance of Diagonal </li>
<li> "ILUT": returns an instance of ILUT </li>
<li> "LOCAL SPARSE TRIDIAGONAL SOLVER": returns an instance of LocalSparseTridiagonalSolver </li>
<li> "RBILUK": returns an instance of RBILUK (ILK(k) preconditioner for BlockCrsMatrix) </li>
<li> "RELAXATION": returns an instance of Relaxation </li>
<li> "RILUK": returns an instance of RILUK (ILU(k) preconditioner) </li>
<li> "SPARSE BLOCK RELAXATION": returns an instance of BlockRelaxation with sparse blocks </li>
<li> "TRIDIAGONAL RELAXATION": returns an instance of BlockRelaxation with tridiagonal matrix blocks </li>
</ul>

*/
template<class MatrixType>
class OneLevelFactory {
public:
  typedef typename MatrixType::scalar_type scalar_type;
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  typedef typename MatrixType::node_type node_type;
  typedef ::Ifpack2::Preconditioner<scalar_type,
                                    local_ordinal_type,
                                    global_ordinal_type,
                                    node_type> prec_type;
  typedef ::Tpetra::RowMatrix<scalar_type,
                              local_ordinal_type,
                              global_ordinal_type,
                              node_type> row_matrix_type;

  static_assert (std::is_same<MatrixType, row_matrix_type>::value,
                 "Ifpack2::Details::OneLevelFactory: MatrixType must be a "
                 "Tpetra::RowMatrix specialization.");

  /** \brief Create an instance of Preconditioner given the string
   * name of the preconditioner type.
   *
   * \param precType [in] Name of preconditioner type to be created.
   * \param matrix [in] Matrix used to define the preconditioner
   *
   * Throw an exception if the preconditioner with that input name
   * does not exist.  Otherwise, return a newly created preconditioner
   * object.
   */
  Teuchos::RCP<prec_type>
  create (const std::string& precType,
          const Teuchos::RCP<const row_matrix_type>& matrix) const;
};

} // namespace Details
} // namespace Ifpack2

#endif // IFPACK2_DETAILS_ONELEVELFACTORY_DECL_HPP
