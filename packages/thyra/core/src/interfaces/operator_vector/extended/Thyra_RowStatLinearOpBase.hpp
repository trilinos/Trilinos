// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_ROW_SUM_LINEAR_OP_BASE_HPP
#define THYRA_ROW_SUM_LINEAR_OP_BASE_HPP

#include "Thyra_LinearOpBase_decl.hpp"


namespace Thyra {


namespace RowStatLinearOpBaseUtils {


/** \brief Rows statistic requested. */
enum ERowStat {
  /** \brief Inverse absolute row sums. */
  ROW_STAT_INV_ROW_SUM,
  /** \brief Absolute row sums. */
  ROW_STAT_ROW_SUM,
  /** \brief Inverse absolute column sums. */
  ROW_STAT_INV_COL_SUM,
  /** \brief Absolute column sums. */
  ROW_STAT_COL_SUM
};


} // namespace RowStatLinearOpBaseUtils


/** \brief Interface for exxtracting row statistics as a <tt>VectorBase</tt>
 * from a supporting <tt>LinearOpBase</tt> object.
 *
 * \ingroup Thyra_Op_Vec_extended_interfaces_code_grp
 */
template<class Scalar>
class RowStatLinearOpBase : virtual public LinearOpBase<Scalar> {
public:

  /** @name Non-virtual public interface functions. */
  //@{

  /** \brief Determine if a given row stat is supported. */
  bool rowStatIsSupported(
    const RowStatLinearOpBaseUtils::ERowStat rowStat
    ) const
    { return rowStatIsSupportedImpl(rowStat); }

  /** \brief Get some statistics about a supported row.
   *
   * \precondition <tt>this->rowStatIsSupported(rowStat)==true</tt>
   */
  void getRowStat(
    const RowStatLinearOpBaseUtils::ERowStat rowStat,
    const Ptr<VectorBase<Scalar> > &rowStatVec
    ) const
    {
      TEUCHOS_ASSERT(rowStatIsSupported(rowStat));
      getRowStatImpl(rowStat, rowStatVec);
    }

  //@}

protected:

  /** \name Protected virtual functions to be overridden by subclasses. */
  //@{

  /** \brief . */
  virtual bool rowStatIsSupportedImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat) const = 0;

  /** \brief . */
  virtual void getRowStatImpl(
    const RowStatLinearOpBaseUtils::ERowStat rowStat,
    const Ptr<VectorBase<Scalar> > &rowStatVec) const = 0;

  //@}

};


}	// end namespace Thyra


#endif	// THYRA_SCALED_LINEAR_OP_BASE_HPP
