// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BlockedReordering_hpp__
#define __Teko_BlockedReordering_hpp__

#include <string>
#include <vector>

#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpDefaultBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_ProductMultiVectorBase.hpp"

namespace Teko {

/** \brief Class that describes how a flat blocked operator
 *        should be reordered.
 *
 * Class that describes how a flat blocked operator should
 * be reordered. The semenatics are very similar to a
 * Teuchos::ParameterList. Each level expects you to set the
 * size of the level and define what each subblock is. For
 * example, to change the blocked 3x3 matrix from
 *
 * \f$ Z = \left[\begin{array}{ccc}
 *       A & B & C \\
 *       D & E & F \\
 *       G & H & I
 *     \end{array}\right]\f$
 *
 * to
 *
 * \f$ Z' = \left[\begin{array}{cc}
 *     \left[\begin{array}{cc}
 *        A & B \\
 *        D & E
 *     \end{array}\right] &
 *     \left[\begin{array}{c}
 *     C \\ F \end{array}\right] \\
 *     \left[\begin{array}{cc}
 *       G & H
 *     \end{array}\right] & I
 *     \end{array}\right]\f$
 *
 * the algorithm to build the appropriate block manager is
 *
 * <code>
 *    RCP<BlockReorderManager> bmm = rcp(new BlockReorderManager(2));\n
 *    bmm->GetBlock(0)->SetNumBlocks(2);\n
 *    bmm->GetBlock(0)->SetBlock(0,0);\n
 *    bmm->GetBlock(0)->SetBlock(1,1);\n
 *    bmm->SetBlock(1,2);
 * </code>
 *
 * Suppose now you want to take your \f$Z\f$ operator and build
 * \f$Z'\f$. The relevant code is
 *
 * <code>
 *    BlockedLinearOp Z = ... \n
 *    LinearOp Zprime = buildReorderedLinearOp(*bmm,Z);
 * </code>
 */
class BlockReorderManager {
 public:
  //! @name Constructors
  //@{

  //! Basic empty constructor
  BlockReorderManager() : children_(0) {}

  //! Set this level to have size sz
  BlockReorderManager(int sz) : children_(sz, Teuchos::null) {}

  //! Copy constructor
  BlockReorderManager(const BlockReorderManager& bmm) : children_(bmm.children_.size()) {
    for (unsigned int i = 0; i < children_.size(); i++) children_[i] = bmm.children_[i]->Copy();
  }

  //! Do nothing destructor
  virtual ~BlockReorderManager() {}

  //@}

  //! Returns a copy of this object
  virtual Teuchos::RCP<BlockReorderManager> Copy() const {
    return Teuchos::rcp(new BlockReorderManager(*this));
  }

  //! Sets the number of subblocks
  virtual void SetNumBlocks(int sz) {
    children_.clear();
    children_.resize(sz);
  }

  //! Gets the number of subblocks
  virtual int GetNumBlocks() const { return children_.size(); }

  /** \brief Sets the sublock to a specific index value
   *
   * Sets the sublock to a specific index value
   * \param[in] blockIndex Subblock to be set
   * \param[in] reorder    The value of the index in this subblock
   *
   * \pre <code>blockIndex<this->GetNumBlocks()</code>
   */
  virtual void SetBlock(int blockIndex, int reorder);

  /** \brief Set the subblock to a use a particular reorder
   *        manager.
   *
   * Set the subblock to a use a particular reorder
   * manager. This facilitates recursive construction algorithms
   * where the base case is not differentiated.
   *
   * \param[in] blockIndex Subblock to be set
   * \param[in] reorder Reorder manager to be used in this block
   *
   * \pre <code>blockIndex<this->GetNumBlocks()</code>
   */
  virtual void SetBlock(int blockIndex, const Teuchos::RCP<BlockReorderManager>& reorder);

  /** \brief Get a particular block. If there is no block at this
   *        index location return a new one.
   *
   * Get a particular block. If there is no block at this
   * index location return a new one.
   *
   * \param[in] blockIndex The index queried.
   *
   * \returns A pointer to the BlockReorderManager object at this
   *          location, or if there is none, a new BlockReorderManager
   *          is created and returned.
   *
   * \pre <code>blockIndex<this->GetNumBlocks()</code>
   * \post return value is not <code>null</code>
   */
  virtual const Teuchos::RCP<BlockReorderManager> GetBlock(int blockIndex);

  /** \brief Get a particular block. If there is no block at this
   *        index location return <code>Teuchos::null</code>
   *
   * Get a particular block. If there is no block at this
   * index location return <code>Teuchos::null</code>
   *
   * \param[in] blockIndex The index queried.
   *
   * \returns A pointer to the BlockReorderManager object at this
   *          location, or if there is none <code>Teuchos::null</code>
   *          is returned.
   *
   * \pre <code>blockIndex<this->GetNumBlocks()</code>
   */
  virtual const Teuchos::RCP<const BlockReorderManager> GetBlock(int blockIndex) const;

  //! For sanities sake, print a readable string
  virtual std::string toString() const;

  //! Largest index in this manager
  virtual int LargestIndex() const;

 protected:
  //! Definitions of the subblocks.
  std::vector<Teuchos::RCP<BlockReorderManager> > children_;
};

/** A class that corresponds to the leaf, or stopping critera
 * for BlockReorderManager. This class should not be used
 * directly.
 */
class BlockReorderLeaf : public BlockReorderManager {
 public:
  //! @name Constructors
  //@{

  //! Simple constructor that sets the index
  BlockReorderLeaf(int ind) : value_(ind) {}

  //! Copy constructor
  BlockReorderLeaf(const BlockReorderLeaf& brl) : value_(brl.value_) {}
  //@}

  //! Make a copy of this object
  virtual Teuchos::RCP<BlockReorderManager> Copy() const {
    return Teuchos::rcp(new BlockReorderLeaf(*this));
  }

  //! Get the number of subblocks (this one returns 0 b/c its a leaf)
  virtual int GetNumBlocks() const { return 0; }

  //! Set the number of subblocks (this one does nothing b/c its a leaf)
  virtual void SetNumBlocks(int /* sz */) {}

  //! Set the sub block, this does nothing b/c its a leaf
  virtual void SetBlock(int /* blockIndex */, int /* reorder */) {}

  //! Get a particular subblock...this returns null
  virtual const Teuchos::RCP<BlockReorderManager> GetBlock(int /* blockIndex */) {
    return Teuchos::null;
  }

  //! Get a particular subblock...this returns null
  virtual const Teuchos::RCP<const BlockReorderManager> GetBlock(int /* blockIndex */) const {
    return Teuchos::null;
  }

  //! Get the the index that is stored in this block
  int GetIndex() const { return value_; }

  //! Return a string description of this leaf class
  virtual std::string toString() const {
    std::stringstream ss;
    ss << value_;
    return ss.str();
  }

  //! Largest index in this manager
  virtual int LargestIndex() const { return value_; }

 protected:
  using BlockReorderManager::SetBlock;

  //! The value of the index for this leaf
  int value_;

 private:
  BlockReorderLeaf();  // hidden from users
};

/** \brief Use the BlockReorderManager to change a flat square blocked operator
 *        into a composite operator.
 *
 * Use the BlockReorderManager to chanage a flat square blocked operator into
 * a more complex composite structure. The manager should not have any indicies
 * larger then the size of the blocked operator.
 *
 * \param[in] bmm   BlockReorderManager that specifies how the blocked operator
 *                  is to be restructured.
 * \param[in] blkOp The block operator to be reordered and restructured. Only the
 *                  first level of the operator will be considered. Each subblock
 *                  (even if it is itself blocked) will be handed as an individual
 *                  operator.
 *
 * \returns The reordered blocked linear operator.
 *
 * \pre The largest index in <code>bmm</code> is smaller then the dimension of the
 *      <code>blkOp</code>.
 * \pre The opertor is square.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<const Thyra::LinearOpBase<double> > buildReorderedLinearOp(
    const BlockReorderManager& bmm,
    const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> >& blkOp);

/** \brief Use the BlockReorderManager to change a flat blocked operator
 *        into a composite operator.
 *
 * Use the BlockReorderManager to chanage a flat square blocked operator into
 * a more complex composite structure. The manager should not have any indicies
 * larger then the size of the blocked operator.
 *
 * \param[in] rowMgr BlockReorderManager that specifies how the rows are to
 *                   be restructured.
 * \param[in] colMgr BlockReorderManager that specifies how the columns are to
 *                   be restructured.
 * \param[in] blkOp  The block operator to be reordered and restructured. Only the
 *                   first level of the operator will be considered. Each subblock
 *                   (even if it is itself blocked) will be handed as an individual
 *                   operator.
 *
 * \returns The reordered blocked linear operator.
 *
 * \pre The largest index in <code>bmm</code> is smaller then the dimension of the
 *      <code>blkOp</code>.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<const Thyra::LinearOpBase<double> > buildReorderedLinearOp(
    const BlockReorderManager& rowMgr, const BlockReorderManager& colMgr,
    const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> >& blkOp);

/** \brief Use the BlockReorderManager to change a flat vector space
 *        into a composite vector space.
 *
 * Use the BlockReorderManager to chanage a flat vector space
 * a more complex composite structure. The manager should not have any indicies
 * larger then the size of the blocked operator.
 *
 * \param[in] mgr BlockReorderManager that specifies how the space is to
 *                be restructured.
 * \param[in] blkSpc  The block space to be reordered and restructured. Only the
 *                    first level of the space will be considered. Each subspace
 *                    (even if it is itself blocked) will be handed as an individual
 *                    space.
 *
 * \returns The reordered blocked vector space.
 *
 * \pre The largest index in <code>bmm</code> is smaller then the dimension of the
 *      <code>blkSpc</code>.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > buildReorderedVectorSpace(
    const BlockReorderManager& mgr,
    const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> >& blkSpc);

/** \brief Convert a flat multi vector into a reordered multivector.
 *
 * Convert a flat multi vector into a reordered multivector.
 *
 * \param[in] mgr    Block manager describing the reordered multivector.
 * \param[in] blkVec The flat multivector.
 *
 * \returns A reordered multivector structured to be consistent with <code>mgr</code>.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<Thyra::MultiVectorBase<double> > buildReorderedMultiVector(
    const BlockReorderManager& mgr,
    const Teuchos::RCP<Thyra::ProductMultiVectorBase<double> >& blkVec);

/** \brief Convert a flat multi vector into a reordered multivector.
 *
 * Convert a flat multi vector into a reordered multivector.
 *
 * \param[in] mgr    Block manager describing the reordered multivector.
 * \param[in] blkVec The flat multivector.
 *
 * \returns A reordered multivector structured to be consistent with <code>mgr</code>.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<const Thyra::MultiVectorBase<double> > buildReorderedMultiVector(
    const BlockReorderManager& mgr,
    const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> >& blkVec);

/** \brief Convert a reordered multivector into a flat multivector.
 *
 * Convert a reordered multivector into a flat multivector.
 *
 * \param[in] mgr    Block manager describing the reordered multivector.
 * \param[in] blkVec The reordered multivector structured in a way that
 *                   is consistent with by <code>mgr</code>
 *
 * \returns A flattened multivector.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<Thyra::MultiVectorBase<double> > buildFlatMultiVector(
    const BlockReorderManager& mgr,
    const Teuchos::RCP<Thyra::ProductMultiVectorBase<double> >& blkVec);

/** \brief Convert a reordered multivector into a flat multivector.
 *
 * Convert a reordered multivector into a flat multivector.
 *
 * \param[in] mgr    Block manager describing the reordered multivector.
 * \param[in] blkVec The reordered multivector structured in a way that
 *                   is consistent with by <code>mgr</code>
 *
 * \returns A flattened multivector.
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<const Thyra::MultiVectorBase<double> > buildFlatMultiVector(
    const BlockReorderManager& mgr,
    const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> >& blkVec);

/** \brief Convert a reordered vector space into a flat vector space
 */
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > buildFlatVectorSpace(
    const BlockReorderManager& mgr,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double> >& blkSpc);

/** \brief Convert a string to a block reorder manager object
 *
 * Convert a string to a block reorder manager object. These
 * strings have numbers delimted by [,]. For example,
 * the string "[[2 1] 0]" will give a manager with [2 1] in the
 * first block and 0 in the second block.
 *
 * \param[in] reorder Block structure corresponding to the manager
 *
 * \returns A block reorder manager with the requested structure
 *
 * \relates BlockReorderManager
 */
Teuchos::RCP<const BlockReorderManager> blockedReorderFromString(std::string& reorder);

}  // end namespace Teko

#endif
