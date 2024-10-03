// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_ASSERT_OP_HPP
#define THYRA_ASSERT_OP_HPP


#include "Thyra_OperatorVectorTypes.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Teuchos_Assert.hpp"


namespace Thyra {


/* Utility struct for dumping vector space names, dimension etc.
 */
template<class Scalar>
struct dump_vec_spaces_t {
public:
  dump_vec_spaces_t(
    const Thyra::VectorSpaceBase<Scalar>& _vec_space1,
    const std::string &_vec_space1_name,
    const Thyra::VectorSpaceBase<Scalar>& _vec_space2,
    const std::string &_vec_space2_name
    )
    :vec_space1(_vec_space1),vec_space1_name(_vec_space1_name)
    ,vec_space2(_vec_space2),vec_space2_name(_vec_space2_name)
    {}
  const Thyra::VectorSpaceBase<Scalar> &vec_space1;
  const std::string                    vec_space1_name;
  const Thyra::VectorSpaceBase<Scalar> &vec_space2;
  const std::string                    vec_space2_name;
}; // end dum_vec_spaces


/* Utility function for dumping vector space names, dimension etc.
 */
template<class Scalar>
inline dump_vec_spaces_t<Scalar> dump_vec_spaces(
  const Thyra::VectorSpaceBase<Scalar>& vec_space1,
  const std::string &vec_space1_name,
  const Thyra::VectorSpaceBase<Scalar>& vec_space2,
  const std::string &vec_space2_name
  )
{
  return dump_vec_spaces_t<Scalar>(
    vec_space1,vec_space1_name,vec_space2,vec_space2_name);
}


// Notice!!!!!!!  Place a breakpoint in following function in order to halt the
// program just before an exception is thrown!


/* Utility ostream operator for dumping vector space names, dimension etc.
 */
template<class Scalar>
std::ostream& operator<<( std::ostream& o, const dump_vec_spaces_t<Scalar>& d )
{

  using Teuchos::OSTab;
  const Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_MEDIUM;
  o << "Error, the following vector spaces are not compatible:\n\n";
  OSTab(o).o()
    << d.vec_space1_name << " : "
    << Teuchos::describe(d.vec_space1,verbLevel);
  o << "\n";
  OSTab(o).o()
    << d.vec_space2_name << " : "
    << Teuchos::describe(d.vec_space2,verbLevel);
  return o;
}


/* Utility enum for selecting domain or range spaces
 */
enum EM_VS { VS_RANGE, VS_DOMAIN };


/** \brief Utility function for selecting domain or range spaces
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
template<class Scalar>
const Thyra::VectorSpaceBase<Scalar>& linear_op_op(
  const Thyra::LinearOpBase<Scalar>& M,
  Thyra::EOpTransp M_trans,
  EM_VS M_VS
  )
{
  if(real_trans(M_trans) == NOTRANS && M_VS == VS_RANGE)
    return *M.range();
  if(real_trans(M_trans) == TRANS && M_VS == VS_RANGE)
    return *M.domain();
  if(real_trans(M_trans) == NOTRANS && M_VS == VS_DOMAIN)
    return *M.domain();
  // real_trans(M_trans) == TRANS && M_VS == VS_DOMAIN
  return *M.range();
}


} // end namespace Thyra


/** \brief This macro just asserts that a LHS argument is set
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_LHS_ARG(FUNC_NAME,LHS_ARG) \
  TEUCHOS_TEST_FOR_EXCEPTION( \
    (&*LHS_ARG) == NULL, std::invalid_argument, \
    FUNC_NAME << " : Error!" \
    );


// Notice!!!!!!!  Setting a breakpoint at the line number that is printed by this macro
// and then trying to set the condition !isCompatible does not work (at least not
// in gdb).


/** \brief Helper assertion macro
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_VEC_SPACES_NAMES(FUNC_NAME,VS1,VS1_NAME,VS2,VS2_NAME) \
{ \
  const bool l_isCompatible = (VS1).isCompatible(VS2); \
  TEUCHOS_TEST_FOR_EXCEPTION( \
    !l_isCompatible, ::Thyra::Exceptions::IncompatibleVectorSpaces, \
    FUNC_NAME << "\n\n" \
    << ::Thyra::dump_vec_spaces(VS1,VS1_NAME,VS2,VS2_NAME) \
    ) \
}


/** \brief This is a very useful macro that should be used to validate
 * that two vector spaces are compatible.
 *
 * If the vector spaces are not compatible then a very helpful error
 * message is generated (in the std::exception class) along with the
 * file name and line number where this macro is called.  The error
 * message string embedded in the thrown exception gives the concrete
 * types of the two vector spaces involved as well as their dimensions.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_VEC_SPACES(FUNC_NAME,VS1,VS2)\
THYRA_ASSERT_VEC_SPACES_NAMES(FUNC_NAME,VS1,#VS1,VS2,#VS2)


/** \brief This macro validates that a linear operator and a vector space
 * for the domain vector are compatible.
 *
 * This macro is not recommended for casual users.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_MAT_VEC_SPACES(FUNC_NAME,M,M_T,M_VS,VS) \
{ \
  std::ostringstream M_VS_name; \
  M_VS_name << "(" #M << ( (M_T) == Thyra::NOTRANS ? "" : "^T" ) << ")" \
            << "." << ( (M_VS) == Thyra::VS_RANGE ? "range()" : "domain()" ); \
  THYRA_ASSERT_VEC_SPACES_NAMES( \
    FUNC_NAME, \
    ::Thyra::linear_op_op(M,M_T,M_VS),M_VS_name.str().c_str(), \
    (VS),#VS \
    ) \
}


/** \brief This is a very useful macro that should be used to validate
 * that the spaces for the vector version of the
 * <tt>LinearOpBase::apply()</tt> function (or related operations).
 *
 * If the vector spaces are not compatible then a very helpful error
 * message is generated (in the std::exception class) along with the
 * file name and line number where this macro is called.  The error
 * message string embedded in the thrown exception gives the concrete
 * types of the vector spaces involved as well as their dimensions.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_LINEAR_OP_VEC_APPLY_SPACES(FUNC_NAME,M,M_T,X,Y) \
  { \
    std::ostringstream headeross; \
    headeross \
      << FUNC_NAME << ":\n" \
      << "Spaces check failed for " \
      << #M << ( (M_T) == Thyra::NOTRANS ? "" : "^T" ) << " * " \
      << #X << " and " << #Y; \
    const std::string &header = headeross.str(); \
    THYRA_ASSERT_MAT_VEC_SPACES(header,M,M_T,::Thyra::VS_RANGE,*(Y)->space()); \
    THYRA_ASSERT_MAT_VEC_SPACES(header,M,M_T,::Thyra::VS_DOMAIN,*(X).space()); \
  }


/** \brief This is a very useful macro that should be used to validate
 * that the spaces for the multi-vector version of the
 * <tt>LinearOpBase::apply()</tt> function (or related operations).
 *
 * If the vector spaces are not compatible then a very helpful error
 * message is generated (in the std::exception class) along with the
 * file name and line number where this macro is called.  The error
 * message string embedded in the thrown exception gives the concrete
 * types of the vector spaces involved as well as their dimensions.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_LINEAR_OP_MULTIVEC_APPLY_SPACES(FUNC_NAME,M,M_T,X,Y) \
  { \
    std::ostringstream headeross; \
    headeross \
      << FUNC_NAME << ":\n\n" \
      << "Spaces check failed for " \
      << #M << ( (M_T) == Thyra::NOTRANS ? "" : "^T" ) << " * " \
      << #X << " and " << #Y << ":\n\n"; \
    const std::string &header = headeross.str(); \
    THYRA_ASSERT_VEC_SPACES(header,*(X).domain(),*(Y)->domain()); \
    THYRA_ASSERT_MAT_VEC_SPACES(header,M,M_T,::Thyra::VS_RANGE,*(Y)->range()); \
    THYRA_ASSERT_MAT_VEC_SPACES(header,M,M_T,::Thyra::VS_DOMAIN,*(X).range()); \
  }


namespace Thyra {


template<class Scalar>
void assertLinearOpPlusLinearOpNames(
  const std::string &funcName,
  const LinearOpBase<Scalar> &M1, const EOpTransp M1_trans_in, const std::string &M1_name,
  const LinearOpBase<Scalar> &M2, const EOpTransp M2_trans_in, const std::string &M2_name
  )
{
  const EOpTransp M1_trans = real_trans(M1_trans_in);
  const EOpTransp M2_trans = real_trans(M2_trans_in);
  std::ostringstream headeross;
  headeross
    << funcName << ":\n\n"
    << "Spaces check failed for "
    << "(" << M1_name << ")" << ( M1_trans == NOTRANS ? "" : "^T" )
    << " + "
    << "(" << M2_name << ")" << ( M2_trans == NOTRANS ? "" : "^T" )
    << " where:\n\n"
    << "  " << M1_name << ": " << M1.description() << "\n\n"
    << "  " << M2_name << ": " << M2.description();
  const std::string &header = headeross.str();
  if ( M1_trans == M2_trans ) {
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.range(), M1_name + ".range()",
      *M2.range(), M2_name + ".range()" );
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.domain(), M1_name + ".domain()",
      *M2.domain(), M2_name + ".domain()" );
  }
  else { // M1_trans != M2_trans
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.domain(), M1_name + ".domain()",
      *M2.range(), M2_name + ".range()" );
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.range(), M1_name + ".range()",
      *M2.domain(), M2_name + ".domain()" );
  }
}


template<class Scalar>
void assertLinearOpTimesLinearOpNames(
  const std::string &funcName,
  const LinearOpBase<Scalar> &M1, const EOpTransp M1_trans_in, const std::string &M1_name,
  const LinearOpBase<Scalar> &M2, const EOpTransp M2_trans_in, const std::string &M2_name
  )
{
  const EOpTransp M1_trans = real_trans(M1_trans_in);
  const EOpTransp M2_trans = real_trans(M2_trans_in);
  std::ostringstream headeross;
  headeross
    << funcName << ":\n\n"
    << "Spaces check failed for "
    << "(" << M1_name << ")" << ( M1_trans == NOTRANS ? "" : "^T" )
    << " * "
    << "(" << M2_name << ")" << ( M2_trans == NOTRANS ? "" : "^T" )
    << " where:\n\n"
    << "  " << M1_name << ": " << M1.description() << "\n\n"
    << "  " << M2_name << ": " << M2.description();
  const std::string &header = headeross.str();
  if ( M1_trans == NOTRANS &&  M2_trans == NOTRANS ) {
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.domain(), M1_name + ".domain()",
      *M2.range(), M2_name + ".range()" );
  }
  else if ( M1_trans == NOTRANS && M2_trans == TRANS ) {
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.domain(), M1_name + ".domain()",
      *M2.domain(), M2_name + ".domain()" );
  }
  else if ( M1_trans == TRANS &&  M2_trans == NOTRANS ) {
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.domain(), M1_name + ".range()",
      *M2.range(), M2_name + ".range()" );
  }
  else if ( M1_trans == TRANS &&  M2_trans == TRANS ) {
    THYRA_ASSERT_VEC_SPACES_NAMES( header,
      *M1.domain(), M1_name + ".range()",
      *M2.range(), M2_name + ".domain()" );
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION( true, std::logic_error,
      header << "\n\n" << "Error, invalid value for trasponse enums!" );
  }
}


} // namespace Thyra


/** \brief Assert that a linear operator addition matches up.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_LINEAR_OP_PLUS_LINEAR_OP_SPACES_NAMES(FUNC_NAME,M1,M1_T,M1_N,M2,M2_T,M2_N) \
  ::Thyra::assertLinearOpPlusLinearOpNames(FUNC_NAME,M1,M1_T,M1_N,M2,M2_T,M2_N)


/** \brief Assert that a linear operator multiplication matches up.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_LINEAR_OP_TIMES_LINEAR_OP_SPACES_NAMES(FUNC_NAME,M1,M1_T,M1_N,M2,M2_T,M2_N) \
  ::Thyra::assertLinearOpTimesLinearOpNames(FUNC_NAME,M1,M1_T,M1_N,M2,M2_T,M2_N)


/** \brief Helper assertion macro
 *
 * This macro is not recommended for casual users.
 *
 * \ingroup Thyra_Op_Vec_ANA_Development_grp
 */
#define THYRA_ASSERT_MAT_MAT_SPACES(FUNC_NAME,M1,M1_T,M1_VS,M2,M2_T,M2_VS) \
  { \
    std::ostringstream headeross; \
    headeross \
      << FUNC_NAME << "\n" \
      << "Spaces check failed for " \
      << #M1 << ( (M1_T) == Thyra::NOTRANS ? "" : "^T" ) << " and " \
      << #M2 << ( (M2_T) == Thyra::NOTRANS ? "" : "^T" ); \
    const std::string &header = headeross.str(); \
    std::ostringstream M1_VS_name, M2_VS_name; \
    M1_VS_name << "(" #M1 << ( M1_T == ::Thyra::NOTRANS ? "" : "^T" ) << ")" \
               << "." << ( M1_VS == ::Thyra::VS_RANGE ? "range()" : "domain()" ); \
    M2_VS_name << "(" #M2 << ( M2_T == ::Thyra::NOTRANS ? "" : "^T" ) << ")" \
               << "." << ( M2_VS == ::Thyra::VS_RANGE ? "range()" : "domain()" ); \
    THYRA_ASSERT_VEC_SPACES_NAMES( \
      header, \
      ::Thyra::linear_op_op(M1,M1_T,M1_VS),M1_VS_name.str().c_str() \
      ::Thyra::linear_op_op(M2,M2_T,M2_VS),M2_VS_name.str().c_str() \
      ); \
  }


#endif // THYRA_ASSERT_OP_HPP
