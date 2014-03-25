// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
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
