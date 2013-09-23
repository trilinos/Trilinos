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

/// \file Ifpack2_SupportGraph_decl.hpp
/// \brief Declaration of SupportGraph preconditioner

#ifndef IFPACK2_SUPPORTGRAPH_DECL_HPP
#define IFPACK2_SUPPORTGRAPH_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Heap.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/config.hpp>
#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {

/// \class SupportGraph
/// \brief SupportGraph of a Tpetra sparse matrix.
/// \tparam Specialization of Tpetra::CrsMatrix or Tpetra::RowMatrix.
///
/// This class computes a maximum weight spanning tree
/// or multiple trees (forest), of a given sparse matrix
/// represented as a Tpetra::RowMatrix.
///
/// \warning If the matrix is distributed over multiple MPI processes,
///   this class will not work correctly by itself.  You must use it
///   as a subdomain solver inside of a domain decomposition method
///   like AdditiveSchwarz (which see).  If you use Factory to create
///   an SupportGraph preconditioner, the Factory will automatically wrap
///   SupportGraph in AdditiveSchwarz for you, if the matrix's communicator
///   contains multiple processes.
///
/// See the documentation of setParameters() for a list of valid
/// parameters.
///
template<class MatrixType>
class SupportGraph :
    virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,
                                           typename MatrixType::local_ordinal_type,
                                           typename MatrixType::global_ordinal_type,
                                           typename MatrixType::node_type>
{

public:
  //! \name Typedefs
  //@{

  //! The type of the entries of the input MatrixType.
  typedef typename MatrixType::scalar_type scalar_type;

  //! Preserved only for backwards compatibility.  Please use "scalar_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;


  //! The type of local indices in the input MatrixType.
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;


  //! The type of global indices in the input MatrixType.
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;

  //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;


  //! The type of the Kokkos Node used by the input MatrixType.
  typedef typename MatrixType::node_type node_type;

  //! Preserved only for backwards compatibility.  Please use "node_type".
  TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


  //! The type of the magnitude (absolute value) of a matrix entry.
  typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;

  //! Preserved only for backwards compatibility.  Please use "magnitude_type".
  TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

  //@}
  //! \name Constructors and Destructors
  //@{

  /// \brief Constructor
  ///
  /// \param A [in] The sparse matrix to factor, as a
  ///   Tpetra::RowMatrix.  (Tpetra::CrsMatrix inherits from this, so
  ///   you may use a Tpetra::CrsMatrix here instead.)
  ///
  /// The factorization will <i>not</i> modify the input matrix.  It
  /// stores the L and U factors in the incomplete factorization
  /// separately.
  explicit SupportGraph(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &A);

  //! Destructor
  virtual ~SupportGraph();

  //@}
  //! \name Methods for setting up and computing the incomplete factorization
  //@{

  /// \brief Set preconditioner parameters.
  ///
  /// SupportGraph implements the following parameters:
  /// <ul>
  /// <li> "fact: absolute threshold" (\c magnitude_type)
  /// <li> "fact: relative threshold" (\c magnitude_type)
  /// </ul>
  /// The absolute and relative threshold parameters affect how this
  /// code modifies the diagonal entry of the output factor.
  ///
 
  void setParameters (const Teuchos::ParameterList& params);


  void findSupport ();


  /// \brief Clear any previously computed factors.
  ///
  /// You may call this before calling compute().  The compute()
  /// method will call this automatically if it has not yet been
  /// called.  If you call this after calling compute(), you must
  /// recompute the factorization (by calling compute() again) before
  /// you may call apply().
  void initialize ();

  //! Returns \c true if the preconditioner has been successfully initialized.
  inline bool isInitialized() const {
    return(IsInitialized_);
  }

  //! Compute factor by having Cholmod do a complete factorization.
  /*! This function computes the SupportGraph factor using the current:
    <ol>
    <li> Value for the \e a \e priori diagonal threshold values.
    </ol>
   */
  void compute();

  //! If compute() is completed, this query returns true, otherwise it returns false.
  inline bool isComputed() const {
    return(IsComputed_);
  }

  //@}

  //! @name Methods implementing Tpetra::Operator.
  //@{

  //! Returns the result of a SupportGraph forward/back solve on a Tpetra::MultiVector X in Y.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
    \param
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(
      const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS,
               scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Tpetra::Map representing the domain of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getDomainMap() const;

  //! Tpetra::Map representing the range of this operator.
  Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> > getRangeMap() const;

  //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
  bool hasTransposeApply() const;

  //@}

  //@{
  //! \name Mathematical functions.

  //! Returns the result of a SupportGraph forward/back solve on a Tpetra::MultiVector X in Y.
  /*!
    \param
    X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
    \param
    Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing result.
  */
  template<class DomainScalar, class RangeScalar>
  void applyTempl(
      const Tpetra::MultiVector<DomainScalar,local_ordinal_type,global_ordinal_type,node_type>& X,
            Tpetra::MultiVector<RangeScalar,local_ordinal_type,global_ordinal_type,node_type>& Y,
            Teuchos::ETransp mode = Teuchos::NO_TRANS,
               RangeScalar alpha = Teuchos::ScalarTraits<scalar_type>::one(),
               RangeScalar beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;

  //! Computes the estimated condition number and returns the value.
  magnitude_type computeCondEst(CondestType CT = Cheap,
                               local_ordinal_type MaxIters = 1550,
                               magnitude_type Tol = 1e-9,
                               const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &Matrix_in = Teuchos::null);

  //! Returns the computed estimated condition number, or -1.0 if no computed.
  magnitude_type getCondEst() const { return Condest_; }

  //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
  Teuchos::RCP<const Teuchos::Comm<int> > getComm() const;

  //! Returns a reference to the matrix to be preconditioned.
  Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const;


  //! Returns the number of calls to Initialize().
  int getNumInitialize() const;

  //! Returns the number of calls to Compute().
  int getNumCompute() const;

  //! Returns the number of calls to apply().
  int getNumApply() const;

  //! Returns the time spent in Initialize().
  double getInitializeTime() const;

  //! Returns the time spent in Compute().
  double getComputeTime() const;

  //! Returns the time spent in apply().
  double getApplyTime() const;

  
  //! Get absolute threshold value
  inline magnitude_type getAbsoluteThreshold() const {
    return(Athresh_);
  }

  //! Get relative threshold value
  inline magnitude_type getRelativeThreshold() const {
    return(Rthresh_);
  }

  
  
  // @}

  //! @name Overridden from Teuchos::Describable
  //@{

  /** \brief Return a simple one-line description of this object. */
  std::string description() const;

  /** \brief Print the object with some verbosity level to an FancyOStream object. */
  void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;

  //@}

private:
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;
  typedef typename Teuchos::Array<local_ordinal_type>::size_type size_type;

  //@{ Internal methods

  //! Copy constructor (declared private and undefined; may not be used)
  SupportGraph(const SupportGraph<MatrixType>& RHS);

  //! operator= (declared private and undefined; may not be used)
  SupportGraph<MatrixType>& operator=(const SupportGraph<MatrixType>& RHS);

  //@}
  // \name Internal data
  //@{

  //! reference to the matrix to be preconditioned.
  const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > A_;
  
  Teuchos::RCP<MatrixType> Support_;

  Teuchos::RCP<Amesos2::Solver<MatrixType, Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > solver;

  //@}
  // \name Parameters (set by the input ParameterList)
  //@{

  magnitude_type Athresh_; //!< Absolute threshold
  magnitude_type Rthresh_; //!< Relative threshold
  //! Condition number estimate
  magnitude_type Condest_;
  int Randomize_;
  int NumForests_;
  double KeepDiag_;
  //@}
  // \name Other internal data
  //@{

  
  //! Total time in seconds for all successful calls to initialize().
  double InitializeTime_;
  //! Total time in seconds for all successful calls to compute().
  double ComputeTime_;
  //! Total timer in seconds for all successful calls to apply().
  mutable double ApplyTime_;
  //! The number of successful calls to initialize().
  int NumInitialize_;
  //! The number of successful call to compute().
  int NumCompute_;
  //! The number of successful call to apply().
  mutable int NumApply_;
  //! \c true if \c this object has been initialized
  bool IsInitialized_;
  //! \c true if \c this object has been computed
  bool IsComputed_;

  //@}
}; // class SupportGraph

} // namespace Ifpack2

#endif /* IFPACK2_SUPPORTGRAPH_HPP */
