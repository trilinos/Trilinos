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

#ifndef IFPACK2_ADDITIVESCHWARZ_DEF_HPP
#define IFPACK2_ADDITIVESCHWARZ_DEF_HPP

#include "Ifpack2_AdditiveSchwarz_decl.hpp"

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
#include "Xpetra_RowMatrix.hpp"
#include "Xpetra_TpetraRowMatrix.hpp"
#include "Zoltan2_XpetraRowMatrixInput.hpp"
#include "Zoltan2_OrderingProblem.hpp"
#include "Zoltan2_OrderingSolution.hpp"
#endif

#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Details_CanChangeMatrix.hpp"
#include "Ifpack2_LocalFilter_def.hpp"
#include "Ifpack2_OverlappingRowMatrix_def.hpp"
#include "Ifpack2_Parameters.hpp"
#include "Ifpack2_ReorderFilter_def.hpp"
#include "Ifpack2_SingletonFilter_def.hpp"

#ifdef HAVE_MPI
#include "Teuchos_DefaultMpiComm.hpp"
#endif

#include <locale> // std::ctype::toupper

namespace Ifpack2 {

template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A) :
  Matrix_ (A),
  IsInitialized_(false),
  IsComputed_(false),
  IsOverlapping_(false),
  OverlapLevel_ (0),
  CombineMode_(Tpetra::ADD),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  ComputeCondest_(true),
  UseReordering_(false),
  ReorderingAlgorithm_("none"),
  UseSubdomain_(false),
  FilterSingletons_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0)
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;

  RCP<const Teuchos::Comm<int> > comm = Matrix_->getComm ();
  RCP<const map_type> rowMap = Matrix_->getRowMap ();
  RCP<node_type> node = Matrix_->getNode ();
  const global_size_t INVALID =
    Teuchos::OrdinalTraits<global_size_t>::invalid ();

  // If there's only one process in the matrix's communicator,
  // then there's no need to compute overlap.
  if (comm->getSize () == 1) {
    OverlapLevel_ = 0;
    IsOverlapping_ = false;
  } else if (OverlapLevel_ != 0) {
    IsOverlapping_ = true;
  }

  if (OverlapLevel_ == 0) {
    const global_ordinal_type indexBase = rowMap->getIndexBase ();

    // FIXME (mfh 28 Sep 2013) I don't understand why this is called a
    // "serial Map."  It's the same Map as the input matrix's row Map!
    // It's also the same Map as "distributed Map"!  I would change it
    // myself, but I don't want to break anything, so I just
    // reformatted the code to comply better with Ifpack2 standards
    // and left the names alone.
    SerialMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));
    DistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));

    RCP<const SerialComm<int> > localComm (new SerialComm<int> ());

    LocalDistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeNumElements (),
                         indexBase, localComm, node));
  }

  // Set parameters to default values
  Teuchos::ParameterList plist;
  setParameters (plist);
}

template<class MatrixType, class LocalInverseType>
AdditiveSchwarz<MatrixType, LocalInverseType>::
AdditiveSchwarz (const Teuchos::RCP<const row_matrix_type>& A,
                 const int overlapLevel) :
  Matrix_ (A),
  IsInitialized_(false),
  IsComputed_(false),
  IsOverlapping_(false),
  OverlapLevel_ (overlapLevel),
  CombineMode_(Tpetra::ADD),
  Condest_ (-Teuchos::ScalarTraits<magnitude_type>::one ()),
  ComputeCondest_(true),
  UseReordering_(false),
  ReorderingAlgorithm_("none"),
  UseSubdomain_(false),
  FilterSingletons_(false),
  NumInitialize_(0),
  NumCompute_(0),
  NumApply_(0),
  InitializeTime_(0.0),
  ComputeTime_(0.0),
  ApplyTime_(0.0)
{
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::SerialComm;

  RCP<const Teuchos::Comm<int> > comm = Matrix_->getComm ();
  RCP<const map_type> rowMap = Matrix_->getRowMap ();
  RCP<node_type> node = Matrix_->getNode ();
  const global_size_t INVALID =
    Teuchos::OrdinalTraits<global_size_t>::invalid ();

  // If there's only one process in the matrix's communicator,
  // then there's no need to compute overlap.
  if (comm->getSize () == 1) {
    OverlapLevel_ = 0;
    IsOverlapping_ = false;
  } else if (OverlapLevel_ != 0) {
    IsOverlapping_ = true;
  }

  if (OverlapLevel_ == 0) {
    const global_ordinal_type indexBase = rowMap->getIndexBase ();

    // FIXME (mfh 28 Sep 2013) I don't understand why this is called a
    // "serial Map."  It's the same Map as the input matrix's row Map!
    // It's also the same Map as "distributed Map"!  I would change it
    // myself, but I don't want to break anything, so I just
    // reformatted the code to comply better with Ifpack2 standards
    // and left the names alone.
    SerialMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));
    DistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeElementList (),
                         indexBase, comm, node));

    RCP<const SerialComm<int> > localComm (new SerialComm<int> ());

    LocalDistributedMap_ =
      rcp (new map_type (INVALID, rowMap->getNodeNumElements (),
                         indexBase, localComm, node));
  }

  // Set parameters to default values
  Teuchos::ParameterList plist;
  setParameters (plist);
}


template<class MatrixType,class LocalInverseType>
AdditiveSchwarz<MatrixType,LocalInverseType>::~AdditiveSchwarz () {}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type > >
AdditiveSchwarz<MatrixType,LocalInverseType>::getDomainMap() const
{
  return Matrix_->getDomainMap ();
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::Map<typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> >
AdditiveSchwarz<MatrixType,LocalInverseType>::getRangeMap () const
{
  return Matrix_->getRangeMap ();
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Tpetra::RowMatrix<typename MatrixType::scalar_type, typename MatrixType::local_ordinal_type, typename MatrixType::global_ordinal_type, typename MatrixType::node_type> > AdditiveSchwarz<MatrixType,LocalInverseType>::getMatrix() const
{
  return Matrix_;
}


template<class MatrixType,class LocalInverseType>
void
AdditiveSchwarz<MatrixType,LocalInverseType>::
apply (const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &X,
       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> &Y,
       Teuchos::ETransp mode,
       scalar_type alpha,
       scalar_type beta) const
{
  using Teuchos::RCP;
  using Teuchos::rcp;

  const scalar_type ZERO = Teuchos::ScalarTraits<scalar_type>::zero ();

  TEUCHOS_TEST_FOR_EXCEPTION(
    ! IsComputed_, std::runtime_error,
    "Ifpack2::AdditiveSchwarz::apply: "
    "isComputed() must be true before you may call apply().");

  TEUCHOS_TEST_FOR_EXCEPTION(
    X.getNumVectors() != Y.getNumVectors(), std::invalid_argument,
    "Ifpack2::AdditiveSchwarz::apply: "
    "X and Y must have the same number of columns.  X has "
    << X.getNumVectors() << " columns, but Y has " << Y.getNumVectors() << ".");

  const size_t numVectors = X.getNumVectors ();

  Time_->reset ();
  Time_->start ();

  RCP<MV> OverlappingX,OverlappingY,Xtmp;

  if (IsOverlapping_) {
    // Setup if we're overlapping
    OverlappingX = rcp (new MV (OverlappingMatrix_->getRowMap (), numVectors));
    OverlappingY = rcp (new MV (OverlappingMatrix_->getRowMap (), numVectors));
    // FIXME (mfh 28 Sep 2013) MV's constructor fills with zeros by default,
    // so there is no need to call putScalar().
    OverlappingY->putScalar (ZERO);
    OverlappingX->putScalar (ZERO);
    OverlappingMatrix_->importMultiVector (X, *OverlappingX, Tpetra::INSERT);
    // FIXME from Ifpack1: Will not work with non-zero starting solutions.
  }
  else {
    Xtmp = rcp (new MV (X));

    MV Serial (SerialMap_, numVectors);
    // Create Import object on demand, if necessary.
    if (SerialImporter_.is_null ()) {
      SerialImporter_ =
        rcp (new import_type (SerialMap_, Matrix_->getDomainMap ()));
    }
    Serial.doImport (*Xtmp, *SerialImporter_, Tpetra::INSERT);

    OverlappingX = rcp (new MV (LocalDistributedMap_, numVectors));
    OverlappingY = rcp (new MV (LocalDistributedMap_, numVectors));

    //OverlappingX->putScalar(0.0);
    //OverlappingY->putScalar(0.0);

    MV Distributed (DistributedMap_, numVectors);
    // Create Import object on demand, if necessary.
    if (DistributedImporter_.is_null ()) {
      DistributedImporter_ =
        rcp (new import_type (DistributedMap_, Matrix_->getDomainMap ()));
    }
    Distributed.doImport (*Xtmp, *DistributedImporter_, Tpetra::INSERT);

    // FIXME (mfh 28 Sep 2013) Please don't call replaceLocalValue()
    // for every entry.  It's important to understand how MultiVector
    // views work.
    Teuchos::ArrayRCP<const scalar_type> values = Distributed.get1dView();
    size_t index = 0;

    for (size_t v = 0; v < numVectors; v++) {
      for (size_t i = 0; i < Matrix_->getRowMap()->getNodeNumElements(); i++) {
        OverlappingX->replaceLocalValue(i, v, values[index]);
        index++;
      }
    }
  }

  if (FilterSingletons_) {
    // process singleton filter
    MV ReducedX (SingletonMatrix_->getRowMap (), numVectors);
    MV ReducedY (SingletonMatrix_->getRowMap (), numVectors);
    SingletonMatrix_->SolveSingletons (*OverlappingX, *OverlappingY);
    SingletonMatrix_->CreateReducedRHS (*OverlappingY, *OverlappingX, ReducedX);

    // process reordering
    if (! UseReordering_) {
      Inverse_->apply (ReducedX, ReducedY);
    }
    else {
      MV ReorderedX (ReducedX);
      MV ReorderedY (ReducedY);
      ReorderedLocalizedMatrix_->permuteOriginalToReordered (ReducedX, ReorderedX);
      Inverse_->apply (ReorderedX, ReorderedY);
      ReorderedLocalizedMatrix_->permuteReorderedToOriginal (ReorderedY, ReducedY);
    }

    // finish up with singletons
    SingletonMatrix_->UpdateLHS (ReducedY, *OverlappingY);
  }
  else {

    // process reordering
    if (! UseReordering_) {
      Inverse_->apply (*OverlappingX, *OverlappingY);
    }
    else {
      MV ReorderedX (*OverlappingX);
      MV ReorderedY (*OverlappingY);
      ReorderedLocalizedMatrix_->permuteOriginalToReordered (*OverlappingX, ReorderedX);
      Inverse_->apply (ReorderedX, ReorderedY);
      ReorderedLocalizedMatrix_->permuteReorderedToOriginal (ReorderedY, *OverlappingY);
    }
  }

  if (IsOverlapping_) {
    OverlappingMatrix_->exportMultiVector (*OverlappingY, Y, CombineMode_);
  }
  else {
    Teuchos::ArrayRCP<const scalar_type> values = OverlappingY->get1dView();
    size_t index = 0;

    // FIXME (mfh 28 Sep 2013) Please don't call replaceLocalValue()
    // for every entry.  It's important to understand how MultiVector
    // views work.
    for (size_t v = 0; v < numVectors; v++) {
      for (size_t i = 0; i < Matrix_->getRowMap()->getNodeNumElements(); i++) {
        Y.replaceLocalValue(i, v, values[index]);
        index++;
      }
    }
  }

  NumApply_++;
  ApplyTime_ += Time_->stop();
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::
setParameters (const Teuchos::ParameterList& plist)
{
  using Teuchos::getIntegralValue;

  // mfh 18 Nov 2013: Ifpack2's setParameters() method annoyingly
  // passes in the input list as const.  This means that we have to
  // copy it before validation.
  List_ = plist;
  List_.validateParameters (* getValidParameters ());

  // mfh 18 Nov 2013: Supplying the current value as the default value
  // ensures "delta" behavior when users pass in new parameters: any
  // unspecified parameters in the new list retain their values in the
  // old list.  This preserves backwards compatiblity with this class'
  // previous behavior.  Note that validateParametersAndSetDefaults()
  // would have different behavior: any parameters not in the new list
  // would get default values, which could be different than their
  // values in the original list.

  ComputeCondest_ = List_.get ("schwarz: compute condest", ComputeCondest_);

  bool gotCombineMode = false;
  try {
    CombineMode_ = Teuchos::getIntegralValue<Tpetra::CombineMode> (List_, "schwarz: combine mode");
    gotCombineMode = true;
  }
  catch (Teuchos::Exceptions::InvalidParameterName&) {
    // The caller didn't provide that parameter.  Just keep the
    // existing value of CombineMode_.
    gotCombineMode = true;
  }
  catch (Teuchos::Exceptions::InvalidParameterType&) {
    // The user perhaps supplied it as an Tpetra::CombineMode enum
    // value.  Let's try again (below).  If it doesn't succeed, we
    // know that the type is wrong, so we can let it throw whatever
    // exception it would throw.
  }
  if (! gotCombineMode) {
    CombineMode_ = List_.get ("schwarz: combine mode", CombineMode_);
  }

  OverlapLevel_ = List_.get ("schwarz: overlap level", OverlapLevel_);
  if (OverlapLevel_ != 0 && Matrix_->getComm ()->getSize () > 1) {
    IsOverlapping_ = true;
  }

  // Will we be doing reordering?  Unlike Ifpack, we'll use a
  // "schwarz: reordering list" to give to Zoltan2.
  UseReordering_ = List_.get ("schwarz: use reordering", UseReordering_);

#if !defined(HAVE_IFPACK2_XPETRA) || !defined(HAVE_IFPACK2_ZOLTAN2)
  TEUCHOS_TEST_FOR_EXCEPTION(
    UseReordering_, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setParameters: You specified \"schwarz: use reordering\" = true.  "
    "This is only valid when Trilinos was built with Ifpack2, Xpetra, and "
    "Zoltan2 enabled.  Either Xpetra or Zoltan2 was not enabled in your build "
    "of Trilinos.");
#endif

  // FIXME (mfh 18 Nov 2013) Now would be a good time to validate the
  // "schwarz: reordering list" parameter list.  Currently, that list
  // gets extracted in setup().

  // Subdomain check
  if (List_.isParameter ("schwarz: subdomain id") && List_.get ("schwarz: subdomain id", -1) > 0) {
    UseSubdomain_ = true;
  }
  TEUCHOS_TEST_FOR_EXCEPTION(
    UseSubdomain_, std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setParameters: You specified the \"schwarz: subdomain id\" parameter, "
    "with a value other than -1.  This parameter is not yet supported.");

  // if true, filter singletons. NOTE: the filtered matrix can still have
  // singletons! A simple example: upper triangular matrix, if I remove
  // the lower node, I still get a matrix with a singleton! However, filter
  // singletons should help for PDE problems with Dirichlet BCs.
  FilterSingletons_ = List_.get ("schwarz: filter singletons", FilterSingletons_);

  // FIXME (mfh 18 Nov 2013) If the inner solver exists, now might be
  // a good time to validate its parameters.
}


template<class MatrixType,class LocalInverseType>
Teuchos::RCP<const Teuchos::ParameterList>
AdditiveSchwarz<MatrixType,LocalInverseType>::
getValidParameters () const
{
  using Teuchos::ParameterList;
  using Teuchos::parameterList;
  using Teuchos::RCP;
  using Teuchos::rcp_const_cast;

  if (validParams_.is_null ()) {
    const int overlapLevel = 0;
    const bool useReordering = false;
    const bool computeCondest = false;
    const bool filterSingletons = false;

    RCP<ParameterList> plist = parameterList ("Ifpack2::AdditiveSchwarz");

    Tpetra::setCombineModeParameter (*plist, "schwarz: combine mode");
    plist->set ("schwarz: overlap level", overlapLevel);
    plist->set ("schwarz: use reordering", useReordering);
    plist->set ("schwarz: compute condest", computeCondest);
    plist->set ("schwarz: filter singletons", filterSingletons);

    // FIXME (mfh 18 Nov 2013) Get valid parameters from inner solver.
    //
    // FIXME (mfh 18 Nov 2013) Get valid parameters from Zoltan2, if
    // Zoltan2 was enabled in the build.

    validParams_ = rcp_const_cast<const ParameterList> (plist);
  }
  return validParams_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::initialize()
{
  using Teuchos::rcp;

  IsInitialized_ = false;
  IsComputed_ = false; // values required
  Condest_ = -Teuchos::ScalarTraits<magnitude_type>::one ();

  if (Time_.is_null ()) {
    Time_ = rcp (new Teuchos::Time ("Ifpack2::AdditiveSchwarz"));
  }

  Time_->start ();

  // compute the overlapping matrix if necessary
  if (IsOverlapping_) {
    if (UseSubdomain_) {
      const int sid = List_.get ("subdomain id", -1);
      OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_, sid));
    } else {
      OverlappingMatrix_ = rcp (new OverlappingRowMatrix<row_matrix_type> (Matrix_, OverlapLevel_));
    }
  }

  // Setup
  setup ();

  // Initialize inverse
  //
  // FIXME (mfh 28 Sep 2013) The "inverse" should have its own sublist
  // in the input ParameterList.  We shouldn't pass AdditiveSchwarz's
  // parameters directly to the "inverse."
  //
  // FIXME (mfh 28 Sep 2013) Why don't we call these methods in setup()?
  Inverse_->setParameters (List_);
  Inverse_->initialize ();

  IsInitialized_ = true;
  NumInitialize_++;
  InitializeTime_ += Time_->stop();
}


template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isInitialized() const
{
  return IsInitialized_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::compute ()
{
  if (! IsInitialized_) {
    initialize ();
  }

  Time_->reset();
  Time_->start();
  IsComputed_ = false;
  Condest_ = -Teuchos::ScalarTraits<magnitude_type>::one ();

  Inverse_->compute ();

  // FIXME (mfh 28 Sep 2013) It's really a bad idea to use barriers to
  // "make the timings consistent."  This hides load imbalance.  If
  // you want good timing information, you should collect min and max
  // timings from all processes and report the entire range.
#if defined (HAVE_IFPACK2_TIMER_BARRIER)
  Matrix_->getDomainMap ()->getComm ()->barrier ();
#endif

  IsComputed_ = true;
  ++NumCompute_;
  ComputeTime_ += Time_->stop();
}

//==============================================================================
// Returns true if the  preconditioner has been successfully computed, false otherwise.
template<class MatrixType,class LocalInverseType>
bool AdditiveSchwarz<MatrixType,LocalInverseType>::isComputed() const
{
  return IsComputed_;
}


template<class MatrixType,class LocalInverseType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
AdditiveSchwarz<MatrixType,LocalInverseType>::
computeCondEst (CondestType CT,
                local_ordinal_type MaxIters,
                magnitude_type Tol,
                const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &Matrix_in)
{
  // The preconditioner must have been computed in order to estimate
  // its condition number.
  if (! isComputed ()) {
    return -Teuchos::ScalarTraits<magnitude_type>::one ();
  }

  Condest_ = Ifpack2::Condest (*this, CT, MaxIters, Tol, Matrix_in);
  return Condest_;
}


template<class MatrixType,class LocalInverseType>
typename Teuchos::ScalarTraits<typename MatrixType::scalar_type>::magnitudeType
AdditiveSchwarz<MatrixType,LocalInverseType>::getCondEst() const
{
  return Condest_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumInitialize() const
{
  return NumInitialize_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumCompute() const
{
  return NumCompute_;
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getNumApply() const
{
  return NumApply_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getInitializeTime() const
{
  return InitializeTime_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getComputeTime() const
{
  return ComputeTime_;
}


template<class MatrixType,class LocalInverseType>
double AdditiveSchwarz<MatrixType,LocalInverseType>::getApplyTime() const
{
  return ApplyTime_;
}


template<class MatrixType,class LocalInverseType>
std::string AdditiveSchwarz<MatrixType,LocalInverseType>::description() const
{
  std::ostringstream oss;
  oss << Teuchos::Describable::description();
  if (isInitialized()) {
    if (isComputed()) {
      oss << "{status = initialized, computed";
    }
    else {
      oss << "{status = initialized, not computed";
    }
  }
  else {
    oss << "{status = not initialized, not computed";
  }
  oss << ", overlap level =" << OverlapLevel_;
  oss << ", subdomain reordering =" << ReorderingAlgorithm_;
  oss << "}";
  return oss.str();
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::describe(Teuchos::FancyOStream &os, const Teuchos::EVerbosityLevel verbLevel) const
{
  using std::endl;

  os << endl;
  os << "================================================================================" << endl;
  os << "Ifpack2::AdditiveSchwarz, overlap level = " << OverlapLevel_ << endl;
  if (CombineMode_ == Tpetra::INSERT)
    os << "Combine mode                          = Insert" << endl;
  else if (CombineMode_ == Tpetra::ADD)
    os << "Combine mode                          = Add" << endl;
  else if (CombineMode_ == Tpetra::REPLACE)
    os << "Combine mode                          = Replace" << endl;
  else if (CombineMode_ == Tpetra::ABSMAX)
    os << "Combine mode                          = AbsMax" << endl;

  os << "Condition number estimate             = " << Condest_ << endl;
  os << "Global number of rows                 = " << Matrix_->getGlobalNumRows() << endl;

  os << endl;
  os << "Phase           # calls   Total Time (s)       Total MFlops     MFlops/s" << endl;
  os << "-----           -------   --------------       ------------     --------" << endl;
  os << "Initialize()    "   << std::setw(5) << getNumInitialize()
     << "  " << std::setw(15) << getInitializeTime() << endl;
  os << "Compute()       "   << std::setw(5) << getNumCompute()
     << "  " << std::setw(15) << getComputeTime() << endl;
  os << "Apply()         "   << std::setw(5) << getNumApply()
     << "  " << std::setw(15) << getApplyTime() <<endl;
  os << "================================================================================" << endl;
  os << endl;
}


template<class MatrixType,class LocalInverseType>
std::ostream& AdditiveSchwarz<MatrixType,LocalInverseType>::print(std::ostream& os) const
{
  Teuchos::FancyOStream fos(Teuchos::rcp(&os,false));
  fos.setOutputToRootOnly(0);
  describe(fos);
  return(os);
}


template<class MatrixType,class LocalInverseType>
int AdditiveSchwarz<MatrixType,LocalInverseType>::getOverlapLevel() const
{
  return OverlapLevel_;
}


template<class MatrixType,class LocalInverseType>
void AdditiveSchwarz<MatrixType,LocalInverseType>::setup ()
{
#ifdef HAVE_MPI
  using Teuchos::MpiComm;
#endif // HAVE_MPI
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::rcpFromRef;

#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
  typedef Xpetra::RowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> XpetraMatrixType;
  typedef Xpetra::TpetraRowMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> XpetraTpetraMatrixType;
#endif

  // Localized version of Matrix_ or OverlappingMatrix_.
  RCP<row_matrix_type> LocalizedMatrix;

  // The "most current local matrix."  At the end of this method, this
  // will be handed off to the inner solver.
  RCP<row_matrix_type> ActiveMatrix;

  // Create localized matrix.
  if (! OverlappingMatrix_.is_null ()) {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Ifpack2::AdditiveSchwarz::setup: subdomain code not yet supported.");
      //
      // FIXME (mfh 18 Nov 2013) btw what's the difference between
      // Ifpack_NodeFilter and Ifpack_LocalFilter?  The former's
      // documentation sounds a lot like what Ifpack2::LocalFilter
      // does.
      //
      //Ifpack2_NodeFilter *tt = new Ifpack2_NodeFilter(OverlappingMatrix_,nodeID); //FIXME
      //LocalizedMatrix = Teuchos::rcp(tt);
    }
    else
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (OverlappingMatrix_));
  }
  else {
    if (UseSubdomain_) {
      //      int sid = List_.get("subdomain id",-1);
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::logic_error,
        "Ifpack2::AdditiveSchwarz::setup: subdomain code not yet supported.");
    }
    else {
      LocalizedMatrix = rcp (new LocalFilter<row_matrix_type> (Matrix_));
    }
  }

  // Sanity check; I don't trust the logic above to have created LocalizedMatrix.
  TEUCHOS_TEST_FOR_EXCEPTION(
    LocalizedMatrix.is_null (), std::logic_error,
    "Ifpack2::AdditiveSchwarz::setup: LocalizedMatrix is null, after the code "
    "that claimed to have created it.  This should never be the case.  Please "
    "report this bug to the Ifpack2 developers.");

  // Mark localized matrix as active
  ActiveMatrix = LocalizedMatrix;

  // Singleton Filtering
  if (FilterSingletons_) {
    SingletonMatrix_ = rcp (new SingletonFilter<row_matrix_type> (LocalizedMatrix));
    ActiveMatrix = SingletonMatrix_;
  }

  // Do reordering
  if (UseReordering_) {
#if defined(HAVE_IFPACK2_XPETRA) && defined(HAVE_IFPACK2_ZOLTAN2)
    // Unlike Ifpack, Zoltan2 does all the dirty work here.
    Teuchos::ParameterList zlist = List_.sublist ("schwarz: reordering list");
    ReorderingAlgorithm_ = List_.get<std::string> ("order_method", "rcm");
    XpetraTpetraMatrixType XpetraMatrix (ActiveMatrix);
    Zoltan2::XpetraRowMatrixInput<XpetraMatrixType> Zoltan2Matrix (rcpFromRef (XpetraMatrix));

    typedef Zoltan2::OrderingProblem<Zoltan2::XpetraRowMatrixInput<XpetraMatrixType> > ordering_problem_type;
#ifdef HAVE_MPI
    // Grab the MPI Communicator and build the ordering problem with that
    MPI_Comm myRawComm;

    RCP<const MpiComm<int> > mpicomm =
      rcp_dynamic_cast<const MpiComm<int> > (ActiveMatrix->getComm ());
    if (mpicomm == Teuchos::null) {
      myRawComm = MPI_COMM_SELF;
    } else {
      myRawComm = * (mpicomm->getRawMpiComm ());
    }
    ordering_problem_type MyOrderingProblem (&Zoltan2Matrix, &zlist, myRawComm);
#else
    ordering_problem_type MyOrderingProblem (&Zoltan2Matrix, &zlist);
#endif
    MyOrderingProblem.solve ();

    // Now create the reordered matrix & mark it as active

    typedef ReorderFilter<row_matrix_type> reorder_filter_type;
    typedef Zoltan2::OrderingSolution<global_ordinal_type, local_ordinal_type> ordering_solution_type;
    ReorderedLocalizedMatrix_ = rcp (new reorder_filter_type (ActiveMatrix, rcp (new ordering_solution_type (*MyOrderingProblem.getSolution ()))));

    ActiveMatrix = ReorderedLocalizedMatrix_;
#else
    // This is a logic_error, not a runtime_error, because
    // setParameters() should have excluded this case already.
    TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::logic_error, "Ifpack2::AdditiveSchwarz::setup: "
      "The Zoltan2 and Xpetra packages must be enabled in order "
      "to support reordering.");
#endif
  }

  innerMatrix_ = ActiveMatrix;

  TEUCHOS_TEST_FOR_EXCEPTION(
    innerMatrix_.is_null (), std::logic_error, "Ifpack2::AdditiveSchwarz::"
    "setup: Inner matrix is null right before constructing inner solver.  "
    "Please report this bug to the Ifpack2 developers.");

  // Construct the inner solver if necessary.  We go through a bit
  // more trouble than usual to do so, because we want to exercise the
  // new setInnerPreconditioner feature.
  if (Inverse_.is_null ()) {
    setInnerPreconditioner (rcp (new LocalInverseType (Teuchos::null)));
  } else {
    setInnerPreconditioner (Inverse_);
  }
}


template<class MatrixType, class LocalInverseType>
void AdditiveSchwarz<MatrixType, LocalInverseType>::
setInnerPreconditioner (const Teuchos::RCP<Preconditioner<scalar_type,
                                                          local_ordinal_type,
                                                          global_ordinal_type,
                                                          node_type> >& innerPrec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(
    innerPrec.is_null (), std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setInnerPreconditioner: Inner preconditioner must be nonnull.");

  // Make sure that the new inner solver knows how to have its matrix changed.
  typedef Details::CanChangeMatrix<row_matrix_type> can_change_type;
  can_change_type* innerSolver = dynamic_cast<can_change_type*> (&*innerPrec);
  TEUCHOS_TEST_FOR_EXCEPTION(
    innerSolver == NULL, std::invalid_argument, "Ifpack2::AdditiveSchwarz::"
    "setInnerPreconditioner: The input preconditioner does not implement the "
    "setMatrix() feature.  Only input preconditioners that inherit from "
    "Ifpack2::Details::CanChangeMatrix implement this feature.");

  // Give the local matrix to the new inner solver.
  innerSolver->setMatrix (innerMatrix_);

  // Set the new inner solver.
  Inverse_ = innerPrec;
}



} // namespace Ifpack2

#endif // IFPACK2_ADDITIVESCHWARZ_DECL_HPP
