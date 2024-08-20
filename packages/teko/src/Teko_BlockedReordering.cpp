// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include <string>
#include <vector>
#include <stack>

#include "Teko_BlockedReordering.hpp"
#include "Teko_Utilities.hpp"

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_StrUtils.hpp"

#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"

using Teuchos::Array;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;

namespace Teko {

void BlockReorderManager::SetBlock(int blockIndex, int reorder) {
  TEUCHOS_ASSERT(blockIndex < (int)children_.size());

  RCP<BlockReorderManager> child = rcp(new BlockReorderLeaf(reorder));

  children_[blockIndex] = child;
}

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
void BlockReorderManager::SetBlock(int blockIndex, const RCP<BlockReorderManager> &reorder) {
  TEUCHOS_ASSERT(blockIndex < (int)children_.size());

  children_[blockIndex] = reorder;
}

const Teuchos::RCP<BlockReorderManager> BlockReorderManager::GetBlock(int blockIndex) {
  TEUCHOS_ASSERT(blockIndex < (int)children_.size());

  if (children_[blockIndex] == Teuchos::null)
    children_[blockIndex] = rcp(new BlockReorderManager());

  return children_[blockIndex];
}

const Teuchos::RCP<const BlockReorderManager> BlockReorderManager::GetBlock(int blockIndex) const {
  TEUCHOS_ASSERT(blockIndex < (int)children_.size());

  return children_[blockIndex];
}

std::string BlockReorderManager::toString() const {
  // build the string by recursively calling each child
  std::stringstream ss;
  ss << "[";
  for (unsigned int i = 0; i < children_.size(); i++) {
    if (children_[i] == Teuchos::null)
      ss << " <NULL> ";
    else
      ss << " " << children_[i]->toString() << " ";
  }
  ss << "]";

  return ss.str();
}

int BlockReorderManager::LargestIndex() const {
  int max = 0;
  for (unsigned int i = 0; i < children_.size(); i++) {
    // see if current child is larger
    if (children_[i] != Teuchos::null) {
      int subMax = children_[i]->LargestIndex();
      max        = max > subMax ? max : subMax;
    }
  }

  return max;
}

Teuchos::RCP<const Thyra::LinearOpBase<double> > buildReorderedLinearOp(
    const BlockReorderManager &bmm,
    const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > &blkOp) {
  return buildReorderedLinearOp(bmm, bmm, blkOp);
}

Teuchos::RCP<const Thyra::LinearOpBase<double> > buildReorderedLinearOp(
    const BlockReorderManager &rowMgr, const BlockReorderManager &colMgr,
    const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > &blkOp) {
  typedef RCP<const BlockReorderManager> BRMptr;

  int rowSz = rowMgr.GetNumBlocks();
  int colSz = colMgr.GetNumBlocks();

  if (rowSz == 0 && colSz == 0) {
    // both are leaf nodes
    const BlockReorderLeaf &rowLeaf = dynamic_cast<const BlockReorderLeaf &>(rowMgr);
    const BlockReorderLeaf &colLeaf = dynamic_cast<const BlockReorderLeaf &>(colMgr);

    // simply return entry in matrix
    Teko::LinearOp linOp = blkOp->getBlock(rowLeaf.GetIndex(), colLeaf.GetIndex());

    // somehow we need to set this operator up
    if (linOp == Teuchos::null) {
      linOp = Thyra::zero(blkOp->productRange()->getBlock(rowLeaf.GetIndex()),
                          blkOp->productDomain()->getBlock(colLeaf.GetIndex()));
    }

    return linOp;
  } else if (rowSz == 0) {
    Teko::BlockedLinearOp reBlkOp = Teko::createBlockedOp();

    // operator will be rowSz by colSz
    reBlkOp->beginBlockFill(1, colSz);

    // fill the column entries
    for (int col = 0; col < colSz; col++) {
      BRMptr colPtr = colMgr.GetBlock(col);

      reBlkOp->setBlock(0, col, buildReorderedLinearOp(rowMgr, *colPtr, blkOp));
    }

    // done building
    reBlkOp->endBlockFill();

    return reBlkOp;
  } else if (colSz == 0) {
    Teko::BlockedLinearOp reBlkOp = Teko::createBlockedOp();

    // operator will be rowSz by colSz
    reBlkOp->beginBlockFill(rowSz, 1);

    // fill the row entries
    for (int row = 0; row < rowSz; row++) {
      BRMptr rowPtr = rowMgr.GetBlock(row);

      reBlkOp->setBlock(row, 0, buildReorderedLinearOp(*rowPtr, colMgr, blkOp));
    }

    // done building
    reBlkOp->endBlockFill();

    return reBlkOp;
  } else {
    Teko::BlockedLinearOp reBlkOp = Teko::createBlockedOp();

    // this is the general case
    TEUCHOS_ASSERT(rowSz > 0);
    TEUCHOS_ASSERT(colSz > 0);

    // operator will be rowSz by colSz
    reBlkOp->beginBlockFill(rowSz, colSz);

    for (int row = 0; row < rowSz; row++) {
      BRMptr rowPtr = rowMgr.GetBlock(row);

      for (int col = 0; col < colSz; col++) {
        BRMptr colPtr = colMgr.GetBlock(col);

        reBlkOp->setBlock(row, col, buildReorderedLinearOp(*rowPtr, *colPtr, blkOp));
      }
    }

    // done building
    reBlkOp->endBlockFill();

    return reBlkOp;
  }
}

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
    const BlockReorderManager &mgr,
    const Teuchos::RCP<const Thyra::ProductVectorSpaceBase<double> > &blkSpc) {
  typedef RCP<const BlockReorderManager> BRMptr;

  int sz = mgr.GetNumBlocks();

  if (sz == 0) {
    // its a  leaf nodes
    const BlockReorderLeaf &leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);

    // simply return entry in matrix
    return blkSpc->getBlock(leaf.GetIndex());
  } else {
    Array<RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces;

    // loop over each row
    for (int i = 0; i < sz; i++) {
      BRMptr blkMgr = mgr.GetBlock(i);

      const RCP<const Thyra::VectorSpaceBase<double> > lvs =
          buildReorderedVectorSpace(*blkMgr, blkSpc);

      vecSpaces.push_back(lvs);
    }

    // build a vector space
    const RCP<const Thyra::DefaultProductVectorSpace<double> > vs =
        Thyra::productVectorSpace<double>(vecSpaces);

    // build the vector
    return vs;
  }
}

/** \brief Convert a flat multi vector into a reordered multivector.
 *
 * Convert a flat multi vector into a reordered multivector.
 */
Teuchos::RCP<Thyra::MultiVectorBase<double> > buildReorderedMultiVector(
    const BlockReorderManager &mgr,
    const Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > &blkVec) {
  typedef RCP<const BlockReorderManager> BRMptr;

  int sz = mgr.GetNumBlocks();

  if (sz == 0) {
    // its a  leaf nodes
    const BlockReorderLeaf &leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);

    // simply return entry in matrix
    return blkVec->getNonconstMultiVectorBlock(leaf.GetIndex());
  } else {
    Array<RCP<Thyra::MultiVectorBase<double> > > multiVecs;
    Array<RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces;

    // loop over each row
    for (int i = 0; i < sz; i++) {
      BRMptr blkMgr = mgr.GetBlock(i);

      const RCP<Thyra::MultiVectorBase<double> > lmv = buildReorderedMultiVector(*blkMgr, blkVec);
      const RCP<const Thyra::VectorSpaceBase<double> > lvs = lmv->range();

      multiVecs.push_back(lmv);
      vecSpaces.push_back(lvs);
    }

    // build a vector space
    const RCP<const Thyra::DefaultProductVectorSpace<double> > vs =
        Thyra::productVectorSpace<double>(vecSpaces);

    // build the vector
    return Thyra::defaultProductMultiVector<double>(vs, multiVecs);
  }
}

/** \brief Convert a flat multi vector into a reordered multivector.
 *
 * Convert a flat multi vector into a reordered multivector.
 */
Teuchos::RCP<const Thyra::MultiVectorBase<double> > buildReorderedMultiVector(
    const BlockReorderManager &mgr,
    const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> > &blkVec) {
  typedef RCP<const BlockReorderManager> BRMptr;

  int sz = mgr.GetNumBlocks();

  if (sz == 0) {
    // its a  leaf nodes
    const BlockReorderLeaf &leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);

    // simply return entry in matrix
    return blkVec->getMultiVectorBlock(leaf.GetIndex());
  } else {
    Array<RCP<const Thyra::MultiVectorBase<double> > > multiVecs;
    Array<RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces;

    // loop over each row
    for (int i = 0; i < sz; i++) {
      BRMptr blkMgr = mgr.GetBlock(i);

      const RCP<const Thyra::MultiVectorBase<double> > lmv =
          buildReorderedMultiVector(*blkMgr, blkVec);
      const RCP<const Thyra::VectorSpaceBase<double> > lvs = lmv->range();

      multiVecs.push_back(lmv);
      vecSpaces.push_back(lvs);
    }

    // build a vector space
    const RCP<const Thyra::DefaultProductVectorSpace<double> > vs =
        Thyra::productVectorSpace<double>(vecSpaces);

    // build the vector
    return Thyra::defaultProductMultiVector<double>(vs, multiVecs);
  }
}

/** Helper function to assist with the non-constant
 * version of buildFlatMultiVector.
 */
void buildNonconstFlatMultiVector(const BlockReorderManager &mgr,
                                  const RCP<Thyra::MultiVectorBase<double> > &blkVec,
                                  Array<RCP<Thyra::MultiVectorBase<double> > > &multivecs,
                                  Array<RCP<const Thyra::VectorSpaceBase<double> > > &vecspaces) {
  int sz = mgr.GetNumBlocks();

  if (sz == 0) {
    // its a  leaf nodes
    const BlockReorderLeaf &leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);
    int index                    = leaf.GetIndex();

    // simply return entry in matrix
    multivecs[index] = blkVec;
    vecspaces[index] = blkVec->range();
  } else {
    const RCP<Thyra::ProductMultiVectorBase<double> > prodMV =
        rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(blkVec);

    // get flattened elements from each child
    for (int i = 0; i < sz; i++) {
      const RCP<Thyra::MultiVectorBase<double> > mv = prodMV->getNonconstMultiVectorBlock(i);
      buildNonconstFlatMultiVector(*mgr.GetBlock(i), mv, multivecs, vecspaces);
    }
  }
}

/** Helper function to assist with the function
 * of the same name.
 */
void buildFlatMultiVector(const BlockReorderManager &mgr,
                          const RCP<const Thyra::MultiVectorBase<double> > &blkVec,
                          Array<RCP<const Thyra::MultiVectorBase<double> > > &multivecs,
                          Array<RCP<const Thyra::VectorSpaceBase<double> > > &vecspaces) {
  int sz = mgr.GetNumBlocks();

  if (sz == 0) {
    // its a  leaf nodes
    const BlockReorderLeaf &leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);
    int index                    = leaf.GetIndex();

    // simply return entry in matrix
    multivecs[index] = blkVec;
    vecspaces[index] = blkVec->range();
  } else {
    const RCP<const Thyra::ProductMultiVectorBase<double> > prodMV =
        rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(blkVec);

    // get flattened elements from each child
    for (int i = 0; i < sz; i++) {
      RCP<const Thyra::MultiVectorBase<double> > mv = prodMV->getMultiVectorBlock(i);
      buildFlatMultiVector(*mgr.GetBlock(i), mv, multivecs, vecspaces);
    }
  }
}

/** \brief Convert a reordered multivector into a flat multivector.
 *
 * Convert a reordered multivector into a flat multivector.
 */
Teuchos::RCP<Thyra::MultiVectorBase<double> > buildFlatMultiVector(
    const BlockReorderManager &mgr,
    const Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > &blkVec) {
  int numBlocks = mgr.LargestIndex() + 1;

  Array<RCP<Thyra::MultiVectorBase<double> > > multivecs(numBlocks);
  Array<RCP<const Thyra::VectorSpaceBase<double> > > vecspaces(numBlocks);

  // flatten everything into a vector first
  buildNonconstFlatMultiVector(mgr, blkVec, multivecs, vecspaces);

  // build a vector space
  const RCP<Thyra::DefaultProductVectorSpace<double> > vs =
      Thyra::productVectorSpace<double>(vecspaces);

  // build the vector
  return Thyra::defaultProductMultiVector<double>(vs, multivecs);
}

/** \brief Convert a reordered multivector into a flat multivector.
 *
 * Convert a reordered multivector into a flat multivector.
 */
Teuchos::RCP<const Thyra::MultiVectorBase<double> > buildFlatMultiVector(
    const BlockReorderManager &mgr,
    const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> > &blkVec) {
  int numBlocks = mgr.LargestIndex() + 1;

  Array<RCP<const Thyra::MultiVectorBase<double> > > multivecs(numBlocks);
  Array<RCP<const Thyra::VectorSpaceBase<double> > > vecspaces(numBlocks);

  // flatten everything into a vector first
  buildFlatMultiVector(mgr, blkVec, multivecs, vecspaces);

  // build a vector space
  const RCP<const Thyra::DefaultProductVectorSpace<double> > vs =
      Thyra::productVectorSpace<double>(vecspaces);

  // build the vector
  return Thyra::defaultProductMultiVector<double>(vs, multivecs);
}

/** Helper function to assist with the function
 * of the same name.
 */
void buildFlatVectorSpace(const BlockReorderManager &mgr,
                          const RCP<const Thyra::VectorSpaceBase<double> > &blkSpc,
                          Array<RCP<const Thyra::VectorSpaceBase<double> > > &vecspaces) {
  int sz = mgr.GetNumBlocks();

  if (sz == 0) {
    // its a  leaf nodes
    const BlockReorderLeaf &leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);
    int index                    = leaf.GetIndex();

    // simply return entry in matrix
    vecspaces[index] = blkSpc;
  } else {
    const RCP<const Thyra::ProductVectorSpaceBase<double> > prodSpc =
        rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<double> >(blkSpc);

    // get flattened elements from each child
    for (int i = 0; i < sz; i++) {
      RCP<const Thyra::VectorSpaceBase<double> > space = prodSpc->getBlock(i);
      buildFlatVectorSpace(*mgr.GetBlock(i), space, vecspaces);
    }
  }
}

/** \brief Convert a reordered vector space into a flat vector space
 */
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > buildFlatVectorSpace(
    const BlockReorderManager &mgr,
    const Teuchos::RCP<const Thyra::VectorSpaceBase<double> > &blkSpc) {
  int numBlocks = mgr.LargestIndex() + 1;

  Array<RCP<const Thyra::VectorSpaceBase<double> > > vecspaces(numBlocks);

  // flatten everything into a vector first
  buildFlatVectorSpace(mgr, blkSpc, vecspaces);

  // build a vector space
  return Thyra::productVectorSpace<double>(vecspaces);
}

////////////////////////////////////////////////////////////////
// The next three functions are useful for parsing the string
// description of a BlockReorderManager.
////////////////////////////////////////////////////////////////

// this function tokenizes a string, breaking out whitespace but saving the
// brackets [,] as special tokens.
static void tokenize(std::string srcInput, std::string whitespace, std::string prefer,
                     std::vector<std::string> &tokens) {
  std::string input = srcInput;
  std::vector<std::string> wsTokens;
  std::size_t endPos = input.length() - 1;
  while (endPos < input.length()) {
    std::size_t next = input.find_first_of(whitespace);

    // get the sub string
    std::string s;
    if (next != std::string::npos) {
      s = input.substr(0, next);

      // break out the old substring
      input = input.substr(next + 1, endPos);
    } else {
      s     = input;
      input = "";
    }

    endPos = input.length() - 1;

    // add it to the WS tokens list
    if (s == "") continue;
    wsTokens.push_back(s);
  }

  for (unsigned int i = 0; i < wsTokens.size(); i++) {
    // get string to break up
    input = wsTokens[i];

    endPos = input.length() - 1;
    while (endPos < input.length()) {
      std::size_t next = input.find_first_of(prefer);

      std::string s = input;
      if (next > 0 && next < input.length()) {
        // get the sub string
        s = input.substr(0, next);

        input = input.substr(next, endPos);
      } else if (next == 0) {
        // get the sub string
        s = input.substr(0, next + 1);

        input = input.substr(next + 1, endPos);
      } else
        input = "";

      // break out the old substring
      endPos = input.length() - 1;

      // add it to the tokens list
      tokens.push_back(s);
    }
  }
}

// this function takes a set of tokens and returns the first "block", i.e. those
// values (including) brackets that correspond to the first block
static std::vector<std::string>::const_iterator buildSubBlock(
    std::vector<std::string>::const_iterator begin, std::vector<std::string>::const_iterator end,
    std::vector<std::string> &subBlock) {
  std::stack<std::string> matched;
  std::vector<std::string>::const_iterator itr;
  for (itr = begin; itr != end; ++itr) {
    subBlock.push_back(*itr);

    // push/pop brackets as they are discovered
    if (*itr == "[")
      matched.push("[");
    else if (*itr == "]")
      matched.pop();

    // found all matching brackets
    if (matched.empty()) return itr;
  }

  TEUCHOS_ASSERT(matched.empty());

  return itr - 1;
}

// This function takes a tokenized vector and converts it to a block reorder manager
static RCP<BlockReorderManager> blockedReorderFromTokens(const std::vector<std::string> &tokens) {
  // base case
  if (tokens.size() == 1)
    return rcp(new Teko::BlockReorderLeaf(Teuchos::StrUtils::atoi(tokens[0])));

  // check first and last character
  TEUCHOS_ASSERT(*(tokens.begin()) == "[")
  TEUCHOS_ASSERT(*(tokens.end() - 1) == "]");

  std::vector<RCP<Teko::BlockReorderManager> > vecRMgr;
  std::vector<std::string>::const_iterator itr = tokens.begin() + 1;
  while (itr != tokens.end() - 1) {
    // figure out which tokens are relevant for this block
    std::vector<std::string> subBlock;
    itr = buildSubBlock(itr, tokens.end() - 1, subBlock);

    // build the child block reorder manager
    vecRMgr.push_back(blockedReorderFromTokens(subBlock));

    // move the iterator one more
    itr++;
  }

  // build the parent reorder manager
  RCP<Teko::BlockReorderManager> rMgr = rcp(new Teko::BlockReorderManager(vecRMgr.size()));
  for (unsigned int i = 0; i < vecRMgr.size(); i++) rMgr->SetBlock(i, vecRMgr[i]);

  return rMgr;
}

////////////////////////////////////////////////////////////////

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
 */
Teuchos::RCP<const BlockReorderManager> blockedReorderFromString(std::string &reorder) {
  // vector of tokens to use
  std::vector<std::string> tokens;

  // manager to be returned

  // build tokens vector
  tokenize(reorder, " \t\n", "[]", tokens);

  // parse recursively and build reorder manager
  Teuchos::RCP<BlockReorderManager> mgr = blockedReorderFromTokens(tokens);

  return mgr;
}

}  // end namespace Teko
