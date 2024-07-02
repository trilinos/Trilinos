// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/// @file Ifpack2_Details_FastILU_Base_def.hpp

#ifndef __IFPACK2_FASTILU_BASE_DEF_HPP__
#define __IFPACK2_FASTILU_BASE_DEF_HPP__

#include <Ifpack2_Details_CrsArrays.hpp>
#include "Tpetra_BlockCrsMatrix.hpp"
#include "Tpetra_BlockCrsMatrix_Helpers.hpp"
#include "Ifpack2_Details_getCrsMatrix.hpp"
#include <KokkosKernels_Utils.hpp>
#include <Kokkos_Timer.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <stdexcept>

namespace Ifpack2
{
namespace Details
{

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
FastILU_Base(Teuchos::RCP<const TRowMatrix> A) :
  mat_(A),
  initFlag_(false),
  computedFlag_(false),
  nInit_(0),
  nComputed_(0),
  nApply_(0),
  initTime_(0.0),
  computeTime_(0.0),
  applyTime_(0.0),
  crsCopyTime_(0.0),
  params_(Params::getDefaults()) {}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getDomainMap () const
{
  return mat_->getDomainMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getRangeMap () const
{
  return mat_->getRangeMap();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
apply (const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &X,
       Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node> &Y,
       Teuchos::ETransp mode,
       Scalar alpha,
       Scalar beta) const
{
  const std::string timerName ("Ifpack2::FastILU::apply");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  Teuchos::TimeMonitor timeMon (*timer);

  if(!isInitialized() || !isComputed())
  {
    throw std::runtime_error(std::string("Called ") + getName() + "::apply() without first calling initialize() and/or compute().");
  }
  if(X.getNumVectors() != Y.getNumVectors())
  {
    throw std::invalid_argument(getName() + "::apply: X and Y have different numbers of vectors (pass X and Y with exactly matching dimensions)");
  }
  if(X.getLocalLength() != Y.getLocalLength())
  {
    throw std::invalid_argument(getName() + "::apply: X and Y have different lengths (pass X and Y with exactly matching dimensions)");
  }
  //zero out applyTime_ now, because the calls to applyLocalPrec() will add to it
  applyTime_ = 0;
  int  nvecs = X.getNumVectors();
  auto nrowsX = X.getLocalLength();
  auto nrowsY = Y.getLocalLength();
  if(nvecs == 1)
  {
    auto x2d = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto y2d = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
    ImplScalarArray x1d (const_cast<ImplScalar*>(x2d.data()), nrowsX);
    ImplScalarArray y1d (const_cast<ImplScalar*>(y2d.data()), nrowsY);

    applyLocalPrec(x1d, y1d);
  }
  else
  {
    //Solve each vector one at a time (until FastILU supports multiple RHS)
    auto x2d = X.getLocalViewDevice(Tpetra::Access::ReadOnly);
    auto y2d = Y.getLocalViewDevice(Tpetra::Access::ReadWrite);
    for(int i = 0; i < nvecs; i++)
    {
      auto xColView1d = Kokkos::subview(x2d, Kokkos::ALL(), i);
      auto yColView1d = Kokkos::subview(y2d, Kokkos::ALL(), i);
      ImplScalarArray x1d (const_cast<ImplScalar*>(xColView1d.data()), nrowsX);
      ImplScalarArray y1d (const_cast<ImplScalar*>(yColView1d.data()), nrowsY);

      applyLocalPrec(x1d, y1d);
    }
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setParameters (const Teuchos::ParameterList& List)
{
  //Params constructor does all parameter validation, and sets default values
  params_ = Params(List, getName());
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
isBlockCrs() const
{
  return params_.blockCrs && params_.blockCrsSize > 1;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
initialize()
{
  const std::string timerName ("Ifpack2::FastILU::initialize");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  Teuchos::TimeMonitor timeMon (*timer);

  if(mat_.is_null())
  {
    throw std::runtime_error(std::string("Called ") + getName() + "::initialize() but matrix was null (call setMatrix() with a non-null matrix first)");
  }

  if (isBlockCrs()) {
    auto crs_matrix = Ifpack2::Details::getCrsMatrix(this->mat_);

    if (params_.fillBlocks) {
      // Create new TCrsMatrix with the new filled data and conver to Bcrs
      auto crs_matrix_block_filled = Tpetra::fillLogicalBlocks(*crs_matrix, params_.blockCrsSize);
      auto bcrs_matrix = Tpetra::convertToBlockCrsMatrix(*crs_matrix_block_filled, params_.blockCrsSize);
      mat_ = bcrs_matrix;
    }
    else {
      // Assume input is already filled, just convert to Bcrs
      auto bcrs_matrix = Tpetra::convertToBlockCrsMatrix(*crs_matrix, params_.blockCrsSize);
      mat_ = bcrs_matrix;
    }
  }

  Kokkos::Timer copyTimer;
  CrsArrayReader<Scalar, ImplScalar, LocalOrdinal, GlobalOrdinal, Node>::getStructure(mat_.get(), localRowPtrsHost_, localRowPtrs_, localColInds_);
  CrsArrayReader<Scalar, ImplScalar, LocalOrdinal, GlobalOrdinal, Node>::getValues(mat_.get(), localValues_, localRowPtrsHost_);
  crsCopyTime_ = copyTimer.seconds();

  if (params_.use_metis)
  {
    assert(!params_.blockCrs);
    const std::string timerNameMetis ("Ifpack2::FastILU::Metis");
    Teuchos::RCP<Teuchos::Time> timerMetis = Teuchos::TimeMonitor::lookupCounter (timerNameMetis);
    if (timerMetis.is_null ()) {
      timerMetis = Teuchos::TimeMonitor::getNewCounter (timerNameMetis);
    }
    Teuchos::TimeMonitor timeMonMetis (*timerMetis);
    #ifdef HAVE_IFPACK2_METIS
    idx_t nrows = localRowPtrsHost_.size() - 1;
    if (nrows > 0) {
      // reorder will convert both graph and perm/iperm to the internal METIS integer type
      metis_perm_  = MetisArrayHost(Kokkos::ViewAllocateWithoutInitializing("metis_perm"),  nrows);
      metis_iperm_ = MetisArrayHost(Kokkos::ViewAllocateWithoutInitializing("metis_iperm"), nrows);

      // copy ColInds to host
      auto localColIndsHost_ = Kokkos::create_mirror_view(localColInds_);
      Kokkos::deep_copy(localColIndsHost_, localColInds_);

      // prepare for calling metis
      idx_t nnz = localColIndsHost_.size();
      MetisArrayHost metis_rowptr;
      MetisArrayHost metis_colidx;

      bool metis_symmetrize = true;
      if (metis_symmetrize) {
        // symmetrize
        using OrdinalArrayMirror = typename OrdinalArray::host_mirror_type;
        KokkosKernels::Impl::symmetrize_graph_symbolic_hashmap<
          OrdinalArrayHost, OrdinalArrayMirror, MetisArrayHost, MetisArrayHost, Kokkos::HostSpace::execution_space>
          (nrows, localRowPtrsHost_, localColIndsHost_, metis_rowptr, metis_colidx);

        // remove diagonals
        idx_t old_nnz = nnz = 0;
        for (idx_t i = 0; i < nrows; i++) {
          for (LocalOrdinal k = old_nnz; k < metis_rowptr(i+1); k++) {
            if (metis_colidx(k) != i) {
              metis_colidx(nnz) = metis_colidx(k);
              nnz++;
            }
          }
          old_nnz = metis_rowptr(i+1);
          metis_rowptr(i+1) = nnz;
        }
      } else {
        // copy and remove diagonals
        metis_rowptr = MetisArrayHost(Kokkos::ViewAllocateWithoutInitializing("metis_rowptr"), nrows+1);
        metis_colidx = MetisArrayHost(Kokkos::ViewAllocateWithoutInitializing("metis_colidx"), nnz);
        nnz = 0;
        metis_rowptr(0) = 0;
        for (idx_t i = 0; i < nrows; i++) {
          for (LocalOrdinal k = localRowPtrsHost_(i); k < localRowPtrsHost_(i+1); k++) {
            if (localColIndsHost_(k) != i) {
              metis_colidx(nnz) = localColIndsHost_(k);
              nnz++;
            }
          }
          metis_rowptr(i+1) = nnz;
        }
      }

      // call metis
      int info = METIS_NodeND(&nrows, metis_rowptr.data(), metis_colidx.data(),
                              NULL, NULL, metis_perm_.data(), metis_iperm_.data());
      if (METIS_OK != info) {
        throw std::runtime_error(std::string("METIS_NodeND returned info = " + info));
      }
    }
    #else
    throw std::runtime_error(std::string("TPL METIS is not enabled"));
    #endif
  }

  initLocalPrec();  //note: initLocalPrec updates initTime
  initFlag_ = true;
  nInit_++;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
isInitialized() const
{
  return initFlag_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
compute()
{
  if(!initFlag_)
  {
    throw std::runtime_error(getName() + ": initialize() must be called before compute()");
  }

  const std::string timerName ("Ifpack2::FastILU::compute");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  Teuchos::TimeMonitor timeMon (*timer);

  //get copy of values array from matrix
  Kokkos::Timer copyTimer;
  CrsArrayReader<Scalar, ImplScalar, LocalOrdinal, GlobalOrdinal, Node>::getValues(mat_.get(), localValues_, localRowPtrsHost_);
  crsCopyTime_ += copyTimer.seconds(); //add to the time spent getting rowptrs/colinds
  computeLocalPrec(); //this updates computeTime_
  computedFlag_ = true;
  nComputed_++;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
isComputed() const
{
  return computedFlag_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getMatrix() const
{
  return mat_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getNumInitialize() const
{
  return nInit_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getNumCompute() const
{
  return nComputed_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getNumApply() const
{
  return nApply_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getInitializeTime() const
{
  return initTime_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getComputeTime() const
{
  return computeTime_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getApplyTime() const
{
  return applyTime_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
double FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
getCopyTime() const
{
  return crsCopyTime_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
checkLocalILU() const
{
  //if the underlying type of this doesn't implement checkLocalILU, it's an illegal operation
  throw std::runtime_error(std::string("Preconditioner type Ifpack2::Details::") + getName() + " doesn't support checkLocalILU().");
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
checkLocalIC() const
{
  //if the underlying type of this doesn't implement checkLocalIC, it's an illegal operation
  throw std::runtime_error(std::string("Preconditioner type Ifpack2::Details::") + getName() + " doesn't support checkLocalIC().");
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
std::string FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const
{
  std::ostringstream os;
  //Output is a YAML dictionary
  os << "\"Ifpack2::Details::" << getName() << "\": {";
  os << "Initialized: " << (isInitialized() ? "true" : "false") << ", ";
  os << "Computed: " << (isComputed() ? "true" : "false") << ", ";
  os << "Sweeps: " << getSweeps() << ", ";
  os << "Triangular solve type: " << getSpTrsvType() << ", ";
  if (getSpTrsvType() == "Fast") {
    os << "# of triangular solve iterations: " << getNTrisol() << ", ";
  }
  if(mat_.is_null())
  {
    os << "Matrix: null";
  }
  else
  {
    os << "Global matrix dimensions: [" << mat_->getGlobalNumRows() << ", " << mat_->getGlobalNumCols() << "]";
    os << ", Global nnz: " << mat_->getGlobalNumEntries();
  }
  return os.str();
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
setMatrix(const Teuchos::RCP<const TRowMatrix>& A)
{
  if(A.is_null())
  {
    throw std::invalid_argument(std::string("Ifpack2::Details::") + getName() + "::setMatrix() called with a null matrix. Pass a non-null matrix.");
  }
  //bmk note: this modeled after RILUK::setMatrix
  if(mat_.get() != A.get())
  {
    mat_ = A;
    initFlag_ = false;
    computedFlag_ = false;
  }
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
typename FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Params
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Params::getDefaults()
{
  Params p;
  p.use_metis = false;
  p.sptrsv_algo = FastILU::SpTRSV::Fast;
  p.nFact = 5;          // # of sweeps for computing fastILU
  p.nTrisol = 5;        // # of sweeps for applying fastSpTRSV
  p.level = 0;          // level of ILU
  p.omega = 1.0;        // damping factor for fastILU
  p.shift = 0;
  p.guessFlag = true;
  p.blockSizeILU = 1;   // # of nonzeros / thread, for fastILU
  p.blockSize = 1;      // # of rows / thread, for SpTRSV
  p.blockCrs = false;   // whether to use block CRS for fastILU
  p.blockCrsSize = 1;   // block size for block CRS
  p.fillBlocks = false; // whether input matrix needs to be filled
  return p;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Params::Params(const Teuchos::ParameterList& pL, std::string precType)
{
  *this = getDefaults();
  //For each parameter, check that if the parameter exists, it has the right type
  //Then get the value and sanity check it
  //If the parameter doesn't exist, leave it as default (from Params::getDefaults())
  //"sweeps" aka nFact
  #define TYPE_ERROR(name, correctTypeName) {throw std::invalid_argument(precType + "::setParameters(): parameter \"" + name + "\" has the wrong type (must be " + correctTypeName + ")");}
  #define CHECK_VALUE(param, member, cond, msg) {if(cond) {throw std::invalid_argument(precType + "::setParameters(): parameter \"" + param + "\" has value " + std::to_string(member) + " but " + msg);}}

  //metis
  if(pL.isParameter("metis"))
  {
    if(pL.isType<bool>("metis"))
      use_metis = pL.get<bool>("metis");
    else
      TYPE_ERROR("metis", "bool");
  }

  if(pL.isParameter("sweeps"))
  {
    if(pL.isType<int>("sweeps"))
    {
      nFact = pL.get<int>("sweeps");
      CHECK_VALUE("sweeps", nFact, nFact < 1, "must have a value of at least 1");
    }
    else
      TYPE_ERROR("sweeps", "int");
  }
  std::string sptrsv_type = "Fast";
  if(pL.isParameter("triangular solve type")) {
    sptrsv_type = pL.get<std::string> ("triangular solve type");
  }
  if (sptrsv_type == "Standard Host") {
    sptrsv_algo = FastILU::SpTRSV::StandardHost;
  } else if (sptrsv_type == "Standard") {
    sptrsv_algo = FastILU::SpTRSV::Standard;
  }

  //"triangular solve iterations" aka nTrisol
  if(pL.isParameter("triangular solve iterations"))
  {
    if(pL.isType<int>("triangular solve iterations"))
    {
      nTrisol = pL.get<int>("triangular solve iterations");
      CHECK_VALUE("triangular solve iterations", nTrisol, nTrisol < 1, "must have a value of at least 1");
    }
    else
      TYPE_ERROR("triangular solve iterations", "int");
  }
  //"level"
  if(pL.isParameter("level"))
  {
    if(pL.isType<int>("level"))
    {
      level = pL.get<int>("level");
    }
    else if(pL.isType<double>("level"))
    {
      //Level can be read as double (like in ILUT), but must be an exact integer
      //Any integer used for level-of-fill can be represented exactly in double (so use exact compare)
      double dval = pL.get<double>("level");
      double ipart;
      double fpart = modf(dval, &ipart);
      level = ipart;
      CHECK_VALUE("level", level, fpart != 0, "must be an integral value");
    }
    else
    {
      TYPE_ERROR("level", "int");
    }
    CHECK_VALUE("level", level, level < 0, "must be nonnegative");
  }
  if(pL.isParameter("damping factor"))
  {
    if(pL.isType<double>("damping factor"))
      omega = pL.get<double>("damping factor");
    else
      TYPE_ERROR("damping factor", "double");
  }
  if(pL.isParameter("shift"))
  {
    if(pL.isType<double>("shift"))
      shift = pL.get<double>("shift");
    else
      TYPE_ERROR("shift", "double");
  }
  //"guess" aka guessFlag
  if(pL.isParameter("guess"))
  {
    if(pL.isType<bool>("guess"))
      guessFlag = pL.get<bool>("guess");
    else
      TYPE_ERROR("guess", "bool");
  }
  //"block size" aka blkSz
  if(pL.isParameter("block size for ILU"))
  {
    if(pL.isType<int>("block size for ILU"))
    {
      blockSizeILU = pL.get<int>("block size for ILU");
      CHECK_VALUE("block size for ILU", blockSizeILU, blockSizeILU < 1, "must have a value of at least 1");
    }
    else
      TYPE_ERROR("block size for ILU", "int");
  }
  //"block size" aka blkSz
  if(pL.isParameter("block size for SpTRSV"))
  {
    if(pL.isType<int>("block size for SpTRSV"))
      blockSize = pL.get<int>("block size for SpTRSV");
    else
      TYPE_ERROR("block size for SpTRSV", "int");
  }
  //"block crs" aka blockCrs
  if(pL.isParameter("block crs"))
  {
    if(pL.isType<bool>("block crs"))
      blockCrs = pL.get<bool>("block crs");
    else
      TYPE_ERROR("block crs", "bool");
  }
  //"block crs block size" aka blockCrsSize
  if(pL.isParameter("block crs block size"))
  {
    if(pL.isType<int>("block crs block size"))
      blockCrsSize = pL.get<int>("block crs block size");
    else
      TYPE_ERROR("block crs block size", "int");
  }
  //"fill blocks for input" aka fillBlocks
  if(pL.isParameter("fill blocks for input"))
  {
    if(pL.isType<bool>("fill blocks for input"))
      blockCrsSize = pL.get<bool>("fill blocks for input");
    else
      TYPE_ERROR("fill blocks for input", "bool");
  }

  #undef CHECK_VALUE
  #undef TYPE_ERROR
}

#define IFPACK2_DETAILS_FASTILU_BASE_INSTANT(S, L, G, N) \
template class Ifpack2::Details::FastILU_Base<S, L, G, N>;

} //namespace Details
} //namespace Ifpack2

#endif
