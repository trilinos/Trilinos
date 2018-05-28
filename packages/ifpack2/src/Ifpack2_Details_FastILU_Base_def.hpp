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

/// @file Ifpack2_Details_FastILU_Base_def.hpp

#ifndef __IFPACK2_FASTILU_BASE_DEF_HPP__ 
#define __IFPACK2_FASTILU_BASE_DEF_HPP__ 

#include <Ifpack2_Details_CrsArrays.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <stdexcept>
#include "Teuchos_TimeMonitor.hpp"


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
  int nvecs = X.getNumVectors();
  if(nvecs == 1)
  {
    auto x2d = X.template getLocalView<execution_space>();
    auto y2d = Y.template getLocalView<execution_space>();
    auto x1d = Kokkos::subview(x2d, Kokkos::ALL(), 0);
    auto y1d = Kokkos::subview(y2d, Kokkos::ALL(), 0);
    applyLocalPrec(x1d, y1d);
  }
  else
  {
    //Solve each vector one at a time (until FastILU supports multiple RHS)
    for(int i = 0; i < nvecs; i++)
    {
      auto Xcol = X.getVector(i);
      auto Ycol = Y.getVector(i);
      auto xColView2d = Xcol->template getLocalView<execution_space>();
      auto yColView2d = Ycol->template getLocalView<execution_space>();
      ScalarArray xColView1d = Kokkos::subview(xColView2d, Kokkos::ALL(), 0);
      ScalarArray yColView1d = Kokkos::subview(yColView2d, Kokkos::ALL(), 0);
      applyLocalPrec(xColView1d, yColView1d);
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
void FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
initialize()
{

  const std::string timerName ("Ifpack2::FastILU::initialize");
  Teuchos::RCP<Teuchos::Time> timer = Teuchos::TimeMonitor::lookupCounter (timerName);
  if (timer.is_null ()) {
    timer = Teuchos::TimeMonitor::getNewCounter (timerName);
  }
  
  if(mat_.is_null())
  {
    throw std::runtime_error(std::string("Called ") + getName() + "::initialize() but matrix was null (call setMatrix() with a non-null matrix first)");
  }
  Kokkos::Impl::Timer copyTimer;
  CrsArrayReader<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getStructure(mat_.get(), localRowPtrsHost_, localRowPtrs_, localColInds_);
  crsCopyTime_ = copyTimer.seconds();
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


  //get copy of values array from matrix
  Kokkos::Impl::Timer copyTimer;
  CrsArrayReader<Scalar, LocalOrdinal, GlobalOrdinal, Node>::getValues(mat_.get(), localValues_, localRowPtrsHost_);
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
  os << "# of triangular solve iterations: " << getNTrisol() << ", ";
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
  typedef Tpetra::RowGraph<LocalOrdinal, GlobalOrdinal, Node> RGraph;
  Teuchos::RCP<const RGraph> aGraph;    //graph of A
  Teuchos::RCP<const RGraph> matGraph;  //graph of current mat_
  try
  {
    aGraph = A->getGraph();
  }
  catch(...)
  {
    aGraph = Teuchos::null;
  }
  if(!mat_.is_null())
  {
    try
    {
      matGraph = mat_->getGraph();
    }
    catch(...)
    {
      matGraph = Teuchos::null;
    }
  }
  //bmk note: this modeled after RILUK::setMatrix
  if(mat_.get() != A.get())
  {
    mat_ = A;
    if(matGraph.is_null() || (matGraph.getRawPtr() != aGraph.getRawPtr()))
    {
      //must assume that matrix's graph changed, so need to copy the structure again in initialize()
      initFlag_ = false;
    }
    computedFlag_ = false;
  }
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
typename FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Params
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Params::getDefaults()
{
  Params p;
  p.nFact = 5;
  p.nTrisol = 1;
  p.level = 0;
  p.omega = 0.5;
  p.shift = 0;
  p.guessFlag = true;
  p.blockSize = 1;
  return p;
}

template<typename Scalar, typename LocalOrdinal, typename GlobalOrdinal, typename Node>
FastILU_Base<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Params::Params(const Teuchos::ParameterList& pL, std::string precType)
{
  *this = getDefaults();
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude;
  //For each parameter, check that if the parameter exists, it has the right type
  //Then get the value and sanity check it 
  //If the parameter doesn't exist, leave it as default (from Params::getDefaults())
  //"sweeps" aka nFact
  #define TYPE_ERROR(name, correctTypeName) {throw std::invalid_argument(precType + "::setParameters(): parameter \"" + name + "\" has the wrong type (must be " + correctTypeName + ")");}
  #define CHECK_VALUE(param, member, cond, msg) {if(cond) {throw std::invalid_argument(precType + "::setParameters(): parameter \"" + param + "\" has value " + std::to_string(member) + " but " + msg);}}
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
  //"damping factor" aka omega -- try both double and magnitude as type
  if(pL.isParameter("damping factor"))
  {
    if(pL.isType<double>("damping factor"))
      omega = pL.get<double>("damping factor");
    else if(pL.isType<magnitude>("damping factor"))
      omega = pL.get<magnitude>("damping factor");
    else
      TYPE_ERROR("damping factor", "double or magnitude_type");
    CHECK_VALUE("damping factor", omega, omega <= 0 || omega > 1, "must be in the range (0, 1]");
  }
  //"shift" -- also try both double and magnitude
  if(pL.isParameter("shift"))
  {
    if(pL.isType<double>("shift"))
      shift = pL.get<double>("shift");
    else if(pL.isType<magnitude>("shift"))
      shift = pL.get<magnitude>("shift");
    else
      TYPE_ERROR("shift", "double or magnitude_type");
    //no hard requirements for shift value so don't
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
  if(pL.isParameter("block size"))
  {
    if(pL.isType<int>("block size"))
      blockSize = pL.get<int>("block size");
    else
      TYPE_ERROR("block size", "int");
  }
  #undef CHECK_VALUE
  #undef TYPE_ERROR
} 

#define IFPACK2_DETAILS_FASTILU_BASE_INSTANT(S, L, G, N) \
template class Ifpack2::Details::FastILU_Base<S, L, G, N>;

} //namespace Details
} //namespace Ifpack2

#endif

