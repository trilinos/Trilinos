// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_TESTING_TOOLS_HPP
#define THYRA_TESTING_TOOLS_HPP

#include "Thyra_TestingToolsDecl.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_AssertOp.hpp"
#include "Teuchos_as.hpp"


template <class Scalar>
typename Teuchos::ScalarTraits<Scalar>::magnitudeType
Thyra::relVectorErr( const VectorBase<Scalar> &v1, const VectorBase<Scalar> &v2 )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
#ifdef TEUCHOS_DEBUG
  THYRA_ASSERT_VEC_SPACES( "relErr(v1,v2)", *v1.space(), *v2.space() );
#endif
  RCP<VectorBase<Scalar> >
    diff = createMember(v1.space());
  V_VmV( diff.ptr(), v1, v2 );
  const ScalarMag
    nrm_v1 = norm(v1),
    nrm_v2 = norm(v2),
    nrm_diff = norm(*diff);
  return
    ( nrm_diff
      / (
        ST::magnitude(
          Teuchos::RelErrSmallNumber<ST::hasMachineParameters,Scalar>::smallNumber()
          )
        + std::max( nrm_v1, nrm_v2 )
        )
      );
}


template<class Scalar1, class Scalar2, class ScalarMag>
bool Thyra::testRelErrors(
  const std::string &v1_name,
  const ArrayView<const Scalar1> &v1,
  const std::string &v2_name,
  const ArrayView<const Scalar2> &v2,
  const std::string &maxRelErr_error_name,
  const ScalarMag &maxRelErr_error,
  const std::string &maxRelErr_warning_name,
  const ScalarMag &maxRelErr_warning,
  const Ptr<std::ostream> &out,
  const std::string &li
  )
{
  using std::setw;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  typedef typename Teuchos::PromotionTraits<Scalar1,Scalar2>::promote Scalar;

  TEUCHOS_ASSERT_EQUALITY(v1.size(), v2.size());
  const int num_scalars = v1.size();
  
  if(num_scalars==1) {
    return Teuchos::testRelErr<Scalar>(
      v1_name, v1[0], v2_name, v2[0],
      maxRelErr_error_name, maxRelErr_error,
      maxRelErr_warning_name, maxRelErr_warning,
      out
      );
  }

  bool success = true;

  if (nonnull(out)) *out
    << std::endl
    << li << "Check: rel_err(" << v1_name << "," << v2_name << ") <= " << maxRelErr_error_name << " ?\n";

  for( int i = 0; i < num_scalars; ++i ) {
    const ScalarMag rel_err = relErr<Scalar>( v1[i], v2[i] );
    const bool result = ( !SMT::isnaninf(rel_err) && !SMT::isnaninf(maxRelErr_error) && rel_err <= maxRelErr_error );
    if(!result) success = false;
    if(nonnull(out)) {
      *out
        << li << "  "<<setw(2)<<i<<": rel_err("<<v1[i]<<","<<v2[i]<<") "<<"= "<<rel_err
        << " <= " << maxRelErr_error << " : " << passfail(result) << std::endl;
      if( result && rel_err >= maxRelErr_warning ) {
        *out
          << li << "      Warning! rel_err(...) >= " << maxRelErr_warning_name << " = " << maxRelErr_warning << "!\n";
      }
    }
  }

  return success;

}


template<class Scalar>
bool Thyra::testRelNormDiffErr(
  const std::string &v1_name,
  const VectorBase<Scalar> &v1,
  const std::string &v2_name,
  const VectorBase<Scalar> &v2,
  const std::string &maxRelErr_error_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &maxRelErr_error,
  const std::string &maxRelErr_warning_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &maxRelErr_warning,
  std::ostream *out_inout,
  const Teuchos::EVerbosityLevel verbLevel,
  const std::string &li
  )
{
  using std::endl;
  using Teuchos::as;
  using Teuchos::OSTab;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  const RCP<FancyOStream> out =
    Teuchos::fancyOStream(Teuchos::rcp(out_inout, false));
  const ScalarMag
    nrm_v1 = norm(v1),
    nrm_v2 = norm(v2);
  const ScalarMag rel_err = relVectorErr(v1,v2);
  const bool success =
    (
      !SMT::isnaninf(rel_err)
      && !SMT::isnaninf(maxRelErr_error)
      && rel_err <= maxRelErr_error
      );
  if (nonnull(out)) {
    *out
      << endl
      << li << "Testing relative error between vectors "
      << v1_name << " and " << v2_name << ":\n";
    OSTab tab(out);
    *out
      << li << "||"<<v1_name<<"|| = " << nrm_v1 << endl
      << li << "||"<<v2_name<<"|| = " << nrm_v2 << endl;
    if (as<int>(verbLevel) >= as<int>(Teuchos::VERB_HIGH)) {
      *out
        << li << v1_name << " = " << describe(v1,verbLevel)
        << li << v2_name << " = " << describe(v2,verbLevel);
      RCP<VectorBase<Scalar> > diff = createMember(v1.space());
      V_VmV( diff.ptr(), v1, v2 );
      *out
        << li << v1_name << " - " << v2_name << " = " << describe(*diff,verbLevel);
    }
    *out
      << li << "Check: rel_err(" << v1_name << "," << v2_name << ") = "
      << rel_err << " <= " << maxRelErr_error_name << " = "
      << maxRelErr_error << " : " << passfail(success) << endl;
    if( success && rel_err >= maxRelErr_warning ) {
      *out
        << li << "Warning! rel_err(" << v1_name << "," << v2_name << " >= "
        << maxRelErr_warning_name << " = " << maxRelErr_warning << "!\n";
    }
  }
  return success;
}


template<class Scalar>
bool Thyra::testMaxErr(
  const std::string                                             &error_name
  ,const Scalar                                                 &error
  ,const std::string                                            &max_error_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &max_error
  ,const std::string                                            &max_warning_name
  ,const typename Teuchos::ScalarTraits<Scalar>::magnitudeType  &max_warning
  ,std::ostream                                                 *out
  ,const std::string                                            &li
  )
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;
  const ScalarMag error_mag = ST::magnitude(error);
  const bool success = (
    !SMT::isnaninf(error_mag)
    && !SMT::isnaninf(max_error)
    && error_mag <= max_error );
  if(out) {
    *out
      << std::endl
      << li << "Check: |" << error_name << "| = " << error_mag
      << " <= " << max_error_name << " = " << max_error << " : "
      << passfail(success) << std::endl;
    if( success && error_mag >= max_warning ) {
      *out
        << li << "Warning! " << error_name << " = " << error_mag
        << " >= " << max_warning_name << " = " << max_warning << "!\n";
    }
  }
  return success;
}


template<class Scalar>
bool Thyra::testMaxErrors(
  const std::string &error_name,
  const ArrayView<const typename Teuchos::ScalarTraits<Scalar>::magnitudeType> &errors,
  const std::string &max_error_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &max_error,
  const std::string &max_warning_name,
  const typename Teuchos::ScalarTraits<Scalar>::magnitudeType &max_warning,
  const Ptr<std::ostream> &out,
  const std::string &li
  )
{
  using std::setw;
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType ScalarMag;
  typedef Teuchos::ScalarTraits<ScalarMag> SMT;

  const int num_scalars = errors.size();

  if(num_scalars==1) {
    return testMaxErr<Scalar>(
      error_name, errors[0],
      max_error_name, max_error,
      max_warning_name, max_warning,
      out.get(), li
      );
  }

  bool success = true;

  if (nonnull(out)) *out
    << std::endl
    << li << "Check: |"<<error_name<<"| <= "<<max_error_name<<" ?\n";
  for( int i = 0; i < num_scalars; ++i ) {
    const ScalarMag error_mag = ST::magnitude(errors[i]);
    const bool result = (
      !SMT::isnaninf(error_mag)
      && !SMT::isnaninf(max_error)
      && error_mag <= max_error );
    if(!result) success = false;
    if(nonnull(out)) {
      *out
        << li << "  "<<setw(2)<<i<<": |"<<errors[i]<<"| = "<<error_mag<<" <= "
        <<max_error<<" : "<<passfail(success)<<"\n";
      if( result && error_mag >= max_warning ) {
        *out
          << li << "      Warning! |...| >= "<<max_warning_name<<" = "<<max_warning<<"!\n";
      }
    }
  }
  return success;
}


template<class Scalar>
std::ostream& Thyra::operator<<( std::ostream& o, const VectorBase<Scalar>& v )
{
  return o << Teuchos::describe(v, Teuchos::VERB_EXTREME);
}


template<class Scalar>
std::ostream& Thyra::operator<<( std::ostream& o, const LinearOpBase<Scalar>& M )
{
  return o << Teuchos::describe(M, Teuchos::VERB_EXTREME);
}


#endif // THYRA_TESTING_TOOLS_HPP
