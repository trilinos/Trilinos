/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
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

/// \file Ifpack2_Details_FastILU_Base_decl.hpp
/// \brief Defines Ifpack2::Details::FastILU_Base,
///   the base class that contains all common data and behavior
///   of Ifpack2's FastILU, FastILDL and FastIC preconditioners.

#ifndef __IFPACK2_FASTILU_BASE_DECL_HPP__ 
#define __IFPACK2_FASTILU_BASE_DECL_HPP__ 

#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Kokkos_DefaultNode.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include "Ifpack2_Details_CanChangeMatrix.hpp"

namespace Ifpack2
{
namespace Details
{

/// \class FastILU_Base
/// \brief The base class of the Ifpack2 FastILU wrappers (Filu, Fildl and Fic)
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class FastILU_Base : public Ifpack2::Preconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>,
                       public Ifpack2::Details::CanChangeMatrix<
                         Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
{
  public:
    //! Kokkos device type
    typedef typename Node::device_type device_type;
    //! Kokkos execution space
    typedef typename device_type::execution_space execution_space;
    //! Tpetra row matrix
    typedef Tpetra::RowMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TRowMatrix;
    //! Tpetra CRS matrix
    typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> TCrsMatrix;
    //! Tpetra multivector
    typedef Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> TMultiVec;
    //! Kokkos CRS matrix 
    typedef KokkosSparse::CrsMatrix<Scalar, LocalOrdinal, execution_space> KCrsMatrix;
    //! Array of LocalOrdinal on device
    typedef Kokkos::View<LocalOrdinal *, execution_space> OrdinalArray;
    //! Array of LocalOrdinal on host
    typedef Kokkos::View<LocalOrdinal *, Kokkos::HostSpace> OrdinalArrayHost;
    //! Array of Scalar on device
    typedef Kokkos::View<Scalar *, execution_space> ScalarArray;
    //! Array of Scalar on host
    typedef Kokkos::View<Scalar *, Kokkos::HostSpace> ScalarArrayHost;

    //! Constructor
    FastILU_Base(Teuchos::RCP<const TRowMatrix> mat_);

    //! Get the domain map of the matrix
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getDomainMap () const;

    //! Get the range map of the matrix
    Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
    getRangeMap () const;

    //! Apply the preconditioner
    void
    apply (const TMultiVec& X,
           TMultiVec& Y,
           Teuchos::ETransp mode = Teuchos::NO_TRANS,
           Scalar alpha = Teuchos::ScalarTraits<Scalar>::one(),
           Scalar beta = Teuchos::ScalarTraits<Scalar>::zero()) const;
    //! Validate parameters, and set defaults when parameters are not provided.
    /*! Available parameters:
     *    Parameter Name | Parameter Type | Default | Description
     *    ---------------|----------------|---------|------------
     *    \"sweeps\" | int | 5 | The number of iterations of the preconditioner during apply()
     *    \"triangular solve iterations\" | int | 1 | The number of iterations of the block Jacobi triangular solver during apply()
     *    \"level\" | int | 0 | The level of fill
     *    \"damping factor\" | double | 0.5 | The damping factor used during apply() -- must be between 0 (inclusive) and 1 (exclusive)
     *    \"shift\" | double | 0 | The Manteuffel shift parameter
     *    \"guess\" | bool | true | Whether to create multiple preconditioners at successively lower levels of fill to create the initial guess
     *    \"block size\" | int | 1 | The block size for the block Jacobi iterations
     */
    void setParameters (const Teuchos::ParameterList& List);

    //! Initialize the preconditioner
    void initialize();

    //! Whether initialize() has been called since the last time the matrix's structure was changed
    bool isInitialized() const;

    //! Compute the preconditioner
    void compute();

    //! Whether compute() has been called since the last time the matrix's values or structure were changed
    bool isComputed() const;

    //! Get the current matrix
    Teuchos::RCP<const TRowMatrix> getMatrix() const;

    //! Get the number of times initialize() was called.
    int getNumInitialize() const;

    //! Get the number of times compute() was called.
    int getNumCompute() const;

    //! Get the number of times apply() was called.
    int getNumApply() const;

    //! Get the time spent in the last initialize() call.
    double getInitializeTime() const;

    //! Get the time spent in the last compute() call.
    double getComputeTime() const;

    //! Get the time spent in the last apply() call.
    double getApplyTime() const;

    //! Get the time spent deep copying local 3-array CRS out of the matrix.
    double getCopyTime() const;

    //! Get the "sweeps" parameter
    virtual int getSweeps() const = 0;

    //! Get the "triangular solve iterations" parameter
    virtual int getNTrisol() const = 0;

    //! Verify and print debug information about the underlying ILU preconditioner (only supported if \c this is an Ifpack2::Details::Filu)
    virtual void checkLocalILU() const;

    //! Verify and print debug information about the underlying IC preconditioner
    virtual void checkLocalIC() const;

    //! Return a brief description of the preconditioner, in YAML format
    std::string description() const;

    //! Provide a new matrix
    //! If the A's graph is different from that of the existing matrix, initialize() and compute() will need to be called again before the next apply().
    //! If A's values are different than those of the existing matrix, compute() will need to be called again before the next apply().
    void setMatrix(const Teuchos::RCP<const TRowMatrix>& A);

  protected:
    Teuchos::RCP<const TRowMatrix> mat_;
    bool initFlag_;
    bool computedFlag_;
    int nInit_;
    int nComputed_;
    mutable int nApply_;
    //store the local CRS components
    ScalarArray localValues_;     //set at beginning of compute()
    OrdinalArray localRowPtrs_;   //set in initialize()
    OrdinalArray localColInds_;   //set in initialize()
    OrdinalArrayHost localRowPtrsHost_; //set in initialize() and used to get localValues_ in compute()
    double initTime_;
    double computeTime_;
    mutable double applyTime_;
    double crsCopyTime_;         //total time spent deep copying values, rowptrs, colinds out of mat

    //Store validated parameter values (instead of keeping a ParameterList around)
    struct Params
    {
      Params() {}
      Params(const Teuchos::ParameterList& pL, std::string precType);
      int nFact;
      int nTrisol;
      int level;
      double omega;
      double shift;
      bool guessFlag;
      int blockSize;
      static Params getDefaults();
    };

    Params params_;

    //! Construct the underlying preconditioner (localPrec_) using given params and then call localPrec_->initialize()
    // \pre !mat_.is_null()
    virtual void initLocalPrec() = 0;
    //! Get values array from the matrix and then call compute() on the underlying preconditioner
    virtual void computeLocalPrec() = 0;
    //! Apply the local preconditioner with 1-D views of the local parts of X and Y (one vector only)
    virtual void applyLocalPrec(ScalarArray x, ScalarArray y) const = 0;
    //! Get the name of the underlying preconditioner ("Filu", "Fildl" or "Fic")
    virtual std::string getName() const = 0;
};

} //namespace Details
} //namespace Ifpack2

#endif

