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

#ifndef   __Panzer_BlockedVector_ReadOnly_GlobalEvaluationData_hpp__
#define   __Panzer_BlockedVector_ReadOnly_GlobalEvaluationData_hpp__

///////////////////////////////////////////////////////////////////////////////
//
//  Include Files
//
///////////////////////////////////////////////////////////////////////////////

// Panzer
#include "Panzer_ReadOnlyVector_GlobalEvaluationData.hpp"

///////////////////////////////////////////////////////////////////////////////
//
//  Forward Declarations
//
///////////////////////////////////////////////////////////////////////////////

namespace Thyra
{
  template<typename> class DefaultProductVectorSpace;
}

namespace panzer
{
  /**
   *  \brief This class encapsulates the needs of a gather operation to do a
   *         halo exchange for blocked vectors.
   */
  class BlockedVector_ReadOnly_GlobalEvaluationData
    :
    public ReadOnlyVector_GlobalEvaluationData
  {
    public:

      /**
       *  \brief Default Constructor.
       */
      BlockedVector_ReadOnly_GlobalEvaluationData();

      /**
       *  \brief Copy Constructor.
       *
       *  \param[in] src The object to be copied.
       */
      BlockedVector_ReadOnly_GlobalEvaluationData(
        const BlockedVector_ReadOnly_GlobalEvaluationData& src);

      /**
       *  \brief Initializing Constructor.
       *
       *  \param[in] ghostedSpace A `DefaultProductVectorSpace` corresponding
       *                          to the ghosted vector.
       *  \param[in] ownedSpace   A `DefaultProductVectorSpace` corresponding
       *                          to the owned vector.  It's currently ignored,
       *                          but it's included for future changes.
       *  \param[in] gedBlocks    `GlobalEvaluationData` objects that handle
       *                          each block of the vector.
       */
      BlockedVector_ReadOnly_GlobalEvaluationData(
        const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ghostedSpace,
        const Teuchos::RCP<const Thyra::VectorSpaceBase<double>> ownedSpace,
        const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>>&
          gedBlocks);

      /**
       *  \brief Destructor.
       */
      virtual ~BlockedVector_ReadOnly_GlobalEvaluationData()
      {
      } // end of Destructor

      /**
       *  \brief Initialize this object using the sub-`GlobalEvaluationData`
       *         objects.
       *
       *  You must specify the owned and ghosted spaces.  At completion,
       *  `isInitialized_` will be set to `true`.
       *
       *  \param[in] ghostedSpace A `DefaultProductVectorSpace` corresponding
       *                          to the ghosted vector.
       *  \param[in] ownedSpace   A `DefaultProductVectorSpace` corresponding
       *                          to the owned vector.  It's currently ignored,
       *                          but it's included for future changes.
       *  \param[in] gedBlocks    `GlobalEvaluationData` objects that handle
       *                          each block of the vector.
       */
      void
      initialize(
        const Teuchos::RCP<const Thyra::VectorSpaceBase<double>>& ghostedSpace,
        const Teuchos::RCP<const Thyra::VectorSpaceBase<double>>& ownedSpace,
        const std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>>&
          gedBlocks);

      /**
       *  \brief Is this object initialized?
       *
       *  \returns Whether or not the object is initialized.
       */
      virtual bool
      isInitialized() const
      {
        return isInitialized_;
      } // end of isInitialized()

      /**
       *  \brief Communicate the owned data to the ghosted vector.
       *
       *  For this class, this method does the halo exchange for the vector.
       *
       *  \param[in] mem Not needed for this class, but part of the
       *                 `GlobalEvaluationData` interface.
       */
      virtual void
      globalToGhost(
        int mem);

      /**
       *  \brief Initialize internal data for communication.
       *
       *  This clears out the ghosted vector.
       */
      virtual void
      initializeData();

      /**
       * 	\brief Set the owned vector.
       *
       * 	\param[in] ownedVector A `Thyra::VectorBase<double>` that you would
       * 	                       like to set as the owned vector.
       */
      virtual void
      setOwnedVector(
        const Teuchos::RCP<const Thyra::VectorBase<double>>& ownedVector);

      /**
       *  \brief Get the owned vector.
       *
       *  \returns The owned vector as a `Thyra::VectorBase<double>`.
       */
      virtual Teuchos::RCP<const Thyra::VectorBase<double>>
      getOwnedVector() const;

      /**
       *  \brief Get the ghosted vector.
       *
       *  \returns The ghosted vector as a `Thyra::VectorBase<double>`.
       */
      virtual Teuchos::RCP<Thyra::VectorBase<double>>
      getGhostedVector() const;

      /**
       *  \brief How many blocks do we have?
       *
       *  \returns The number of blocks in this `GlobalEvaluationData` object.
       */
      size_t
      getBlockCount() const
      {
        return gedBlocks_.size();
      } // end of getBlockCount()

      /**
       *  \brief Get the `i`-th block (non const version).
       *
       *  \returns This object's `i`-th `GlobalEvaluationData` block.
       */
      Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>
      getGEDBlock(
        int i)
      {
        return gedBlocks_[i];
      } // end of getGEDBlock()

      /**
       *  \brief Get the `i`-th block (const version).
       *
       *  \returns This object's `i`-th `GlobalEvaluationData` block.
       */
      Teuchos::RCP<const ReadOnlyVector_GlobalEvaluationData>
      getGEDBlock(
        int i) const
      {
        return gedBlocks_[i];
      } // end of getGEDBlock()

      /**
       *  \brief Determine if a Dirichlet adjustment is necessary.
       *
       *  For this class, there's nothing to do because it's read-only.
       *
       *  \returns False.
       */
      bool
      requiresDirichletAdjustment() const
      {
        return false;
      } // end of requiresDirichletAdjustment()

    private:

      /**
       *  \brief A flag indicating whether or not the object has been
       *         initialized.
       */
      bool isInitialized_;

      /**
       *  \brief A `vector` of the `GlobalEvaluationData` blocks.
       */
      std::vector<Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>>
      gedBlocks_;

      /**
       *  \brief The owned vector.
       */
      Teuchos::RCP<const Thyra::VectorBase<double>> ownedVector_;

      /**
       *  \brief The vector space corresponding to the ghosted vector.
       */
      Teuchos::RCP<const Thyra::DefaultProductVectorSpace<double>>
      ghostedSpace_;

  }; // end of class BlockedVector_ReadOnly_GlobalEvaluationData

} // end of namespace panzer

#endif // __Panzer_BlockedVector_ReadOnly_GlobalEvaluationData_hpp__
