// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef   __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__
#define   __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__

// Panzer
#include "Panzer_KokkosUtils_VectorToView.hpp"
#include "Panzer_GlobalEvaluationData.hpp"

// Teuchos
#include "Teuchos_RCP.hpp"

// Thyra
#include "Thyra_VectorBase.hpp"

namespace panzer {

/** This class encapsulates the needs of a gather operation to do a halo
  * exchange.
  */
class ReadOnlyVector_GlobalEvaluationData : public GlobalEvaluationData {
public:

  //! Virtual d
  virtual ~ReadOnlyVector_GlobalEvaluationData() {}

  //! Is this object initialized
  virtual bool isInitialized() const = 0;

  /** For this class, this method does the halo exchange for the
    * vector.
    */
  virtual void globalToGhost(int mem) = 0;

  /** For this class, this method does nothing.
    */
  virtual void ghostToGlobal(int /* mem */) {}

  //! Set the owned vector
  virtual void setOwnedVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & ownedVector) = 0;

  //! Get the owned vector
  virtual Teuchos::RCP<const Thyra::VectorBase<double> > getOwnedVector() const = 0;

  //! Get the ghosted vector
  virtual Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const = 0;

  /**
   *  \brief Element access.
   *
   *  Get the `lid`-th element in this `GlobalEvaluationData`.
   *
   *  \note This will pull the appropriate element out of either the owned or
   *        ghosted vector, depending on the value of `lid`.
   *
   *  \param[in] lid The local ID of the element you'd like to get.
   *
   *  \returns The `lid`-th element in this `GlobalEvaluationData`.
   */
  const double&
  operator[](
    const int& lid) const
  {
    if (lid < static_cast<int>(ownedView_.extent(0)))
      return ownedView_(lid);
    else // if (lid >= static_cast<int>(ownedView_.extent(0)))
      return ghostedView_(lid - ownedView_.extent(0));
  } // end of operator[]()

protected:

  /**
   *  \brief The `Kokkos::View` of the owned vector.
   */
  typename panzer::kokkos_utils::VectorToViewTraits<const Epetra_Vector>::
  View ownedView_;

  /**
   *  \brief The `Kokkos::View` of the ghosted vector.
   */
  typename panzer::kokkos_utils::VectorToViewTraits<Epetra_Vector>::View
  ghostedView_;

}; // end of class ReadOnlyVector_GlobalEvaluationData

} // end of namespace panzer

#endif // __Panzer_ReadOnlyVector_GlobalEvaluationData_hpp__
