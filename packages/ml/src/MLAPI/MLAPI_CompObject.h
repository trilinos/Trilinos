#ifndef MLAPI_COMPOBJECT_H
#define MLAPI_COMPOBJECT_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
\file MLAPI_CompObject.h

\brief Class to count flops.

\author Marzio Sala, D-INFK/ETHZ.

\date Last updated on Mar-06.
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

namespace MLAPI {

/*!
\class CompObject

\brief Class to count flops.

\author Marzio Sala, SNL 9214

\date Last updated on Feb-05.

*/

class CompObject {

public:

  //! Constructor, set counter to 0.0.
  CompObject()
  {
    Flops_ = 0.0;
  }

  //! Destructor.
  ~CompObject() {};

  //! Returns the internal counter of flops.
  inline double GetFlops() const
  {
    return(Flops_);
  }

  //! Sets internal counter to \c Flops.
  inline void SetFlops(double Flops) const
  {
    Flops_ = Flops;
  }

  //! Updates internal counter by summing \c Flops.
  inline void UpdateFlops(double Flops) const
  {
    Flops_ += Flops;
  }

private:

  mutable double Flops_;

}; // class CompObject

} // namespace MLPI

#endif // MLAPI_COMPOBJECT_H
