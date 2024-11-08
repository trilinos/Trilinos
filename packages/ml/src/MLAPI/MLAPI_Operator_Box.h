#ifndef MLAPI_OPERATOR_BOX_H
#define MLAPI_OPERATOR_BOX_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_Operator_Box.h

\brief Wrapper for ML_Operator to make it suitable for RefCountPtr.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#include "ml_operator.h"

namespace MLAPI {

/*!
\class ML_Operator_Box

\brief Simple wrapper for ML_Operator struct.

\author Marzio Sala, SNL 9214.

\date Last updated on Feb-05.
*/

class ML_Operator_Box {

public:
  //! Constructor.
  ML_Operator_Box(ML_Operator* Op, bool Ownership = true)
  {
    Op_ = Op;
    Ownership_ = Ownership;
  }

  //! Destructor.
  ~ML_Operator_Box()
  {
    if (Op_ && Ownership_)
      ML_Operator_Destroy(&Op_);
  }

  //! Returns a pointer to the internally stored ML_Operator.
  ML_Operator* GetData() const
  {
    return(Op_);
  }

private:

  //! Pointer to ML_Operator struct.
  ML_Operator* Op_;
  //! If \c true, the destructor will delete \c Op_.
  bool Ownership_;

};

} // namespace MLAPI

#endif // MLAPI_OPERATOR_BOX_H
