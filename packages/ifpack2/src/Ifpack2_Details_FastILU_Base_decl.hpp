// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// \file Ifpack2_Details_FastILU_Base_decl.hpp
/// \brief Defines Ifpack2::Details::FastILU_Base,
///   the base class that contains all common data and behavior
///   of Ifpack2's FastILU, FastILDL and FastIC preconditioners.

#ifndef __IFPACK2_FASTILU_BASE_DECL_HPP__ 
#define __IFPACK2_FASTILU_BASE_DECL_HPP__ 

#include <Tpetra_RowMatrix.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_KokkosCompat_DefaultNode.hpp>
#include <KokkosSparse_CrsMatrix.hpp>
#include <Ifpack2_Preconditioner.hpp>
#include <Ifpack2_Details_CanChangeMatrix.hpp>
#include <shylu_fastutil.hpp>

#ifdef HAVE_IFPACK2_METIS
#include "metis.h"
#endif

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
    //! Kokkos scalar type
    typedef typename Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>::impl_scalar_type ImplScalar;
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
    typedef typename Kokkos::View<LocalOrdinal *, execution_space>::HostMirror OrdinalArrayHost;
    //! Array of Scalar on device
    typedef Kokkos::View<  ImplScalar *, execution_space>  ImplScalarArray;
    typedef Kokkos::View<      Scalar *, execution_space>      ScalarArray;
    typedef Kokkos::View<const Scalar *, execution_space> ConstScalarArray;
    #ifdef HAVE_IFPACK2_METIS
    typedef Kokkos::View<idx_t*, Kokkos::HostSpace> MetisArrayHost;
    #endif

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

    bool isBlockCrs() const;

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

    //! Get the name of triangular solve algorithm
    virtual std::string getSpTrsvType() const = 0;

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
    ImplScalarArray localValues_; //set at beginning of compute()
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
      bool use_metis;
      FastILU::SpTRSV sptrsv_algo;
      int nFact;
      int nTrisol;
      int level;
      int blkSize;
      double omega;
      double shift;
      bool guessFlag;
      int blockSizeILU;
      int blockSize;
      bool blockCrs;
      int blockCrsSize;
      bool fillBlocks;
      static Params getDefaults();
    };

    Params params_;

    #ifdef HAVE_IFPACK2_METIS
    MetisArrayHost metis_perm_;
    MetisArrayHost metis_iperm_;
    #endif

    //! Construct the underlying preconditioner (localPrec_) using given params and then call localPrec_->initialize()
    // \pre !mat_.is_null()
    virtual void initLocalPrec() = 0;
    //! Get values array from the matrix and then call compute() on the underlying preconditioner
    virtual void computeLocalPrec() = 0;
    //! Apply the local preconditioner with 1-D views of the local parts of X and Y (one vector only)
    virtual void applyLocalPrec(ImplScalarArray x, ImplScalarArray y) const = 0;
    //! Get the name of the underlying preconditioner ("Filu", "Fildl" or "Fic")
    virtual std::string getName() const = 0;
};

} //namespace Details
} //namespace Ifpack2

#endif

