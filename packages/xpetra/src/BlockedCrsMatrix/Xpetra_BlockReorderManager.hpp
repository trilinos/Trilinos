// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef XPETRA_BLOCKREORDERMANAGER_HPP_
#define XPETRA_BLOCKREORDERMANAGER_HPP_

#include <stack>

#include <Teuchos_StrUtils.hpp>

namespace Xpetra {
class BlockReorderManager {
 public:
  //! @name Constructors

  //! Basic empty constructor
  BlockReorderManager()
    : children_(0) {}

  //! Copy constructor
  BlockReorderManager(const BlockReorderManager& bmm)
    : children_(bmm.children_.size()) {
    for (size_t i = 0; i < children_.size(); i++) children_[i] = bmm.children_[i]->Copy();
  }

  //! empty destructor
  virtual ~BlockReorderManager() {}

  //@}

  //! returns copy of this object
  virtual Teuchos::RCP<BlockReorderManager> Copy() const {
    return Teuchos::rcp(new BlockReorderManager(*this));
  }

  //! Sets number of subblocks
  virtual void SetNumBlocks(size_t sz) {
    children_.clear();
    children_.resize(sz);
  }

  //! Returns number of subblocks
  virtual size_t GetNumBlocks() const {
    return children_.size();
  }

  //! \brief Sets the subblock to a specific index value
  /**
   * Sets the subblock to a specific index value
   * \param[in] blockIndex: the subblock to be set
   * \param[in] reorder:    the value of the index of this subblock
   */
  virtual void SetBlock(int blockIndex, int reorder); /* {
     TEUCHOS_ASSERT(blockIndex < (int) children_.size());
     Teuchos::RCP<BlockReorderManager> child = Teuchos::rcp(new BlockReorderLeaf(reorder));
     children_[blockIndex] = child;
   }*/

  //! \brief Sets the subblock to a specific index value
  /**
   * Sets the subblock to a specific index value
   * \param[in] blockIndex: the subblock to be set
   * \param[in] reorder:    reorder manager for nested reordering
   */
  virtual void SetBlock(int blockIndex, const Teuchos::RCP<BlockReorderManager>& reorder); /* {
    TEUCHOS_ASSERT(blockIndex < (int) children_.size());
    children_[blockIndex] = reorder;
  }*/

  //! \brief Get a particular block. If there is no block at this index location return a new one
  virtual const Teuchos::RCP<BlockReorderManager> GetBlock(int blockIndex) {
    TEUCHOS_ASSERT(blockIndex < (int)children_.size());

    if (children_[blockIndex] == Teuchos::null)
      children_[blockIndex] = Teuchos::rcp(new BlockReorderManager());
    return children_[blockIndex];
  }

  //! \brief Get a particular block. If there is no block at this index location return a new one
  virtual const Teuchos::RCP<const BlockReorderManager> GetBlock(int blockIndex) const {
    TEUCHOS_ASSERT(blockIndex < (int)children_.size());
    return children_[blockIndex];
  }

  //! for sanities sake, print a readable string
  virtual std::string toString() const {
    // build the string by recursively calling each child
    std::stringstream ss;
    ss << "[";
    for (size_t i = 0; i < children_.size(); i++) {
      if (children_[i] == Teuchos::null)
        ss << " <NULL> ";
      else
        ss << " " << children_[i]->toString() << " ";
    }
    ss << "]";
    return ss.str();
  }

  //! returns largest index in this BlockReorderManager class
  virtual int LargestIndex() const {
    int max = 0;
    for (size_t i = 0; i < children_.size(); i++) {
      if (children_[i] != Teuchos::null) {
        int subMax = children_[i]->LargestIndex();
        max        = max > subMax ? max : subMax;
      }
    }
    return max;
  }

 protected:
  //! definitions of the subblocks
  std::vector<Teuchos::RCP<BlockReorderManager> > children_;
};

class BlockReorderLeaf : public BlockReorderManager {
 public:
  BlockReorderLeaf(int ind)
    : value_(ind) {}
  BlockReorderLeaf(const BlockReorderLeaf& brl)
    : value_(brl.value_) {}

  virtual Teuchos::RCP<BlockReorderManager> Copy() const {
    return Teuchos::rcp(new BlockReorderLeaf(*this));
  }

  virtual size_t GetNumBlocks() const { return 0; }
  virtual void SetNumBlocks(size_t /* sz */) {}
  virtual void SetBlock(int /* blockIndex */, int /* reorder */) {}
  virtual void SetBlock(int /* blockIndex */, const Teuchos::RCP<BlockReorderManager>& /* reorder */) {}
  virtual const Teuchos::RCP<BlockReorderManager> GetBlock(int /* blockIndex */) {
    return Teuchos::null;
  }
  virtual const Teuchos::RCP<const BlockReorderManager> GetBlock(int /* blockIndex */) const {
    return Teuchos::null;
  }
  //! Get the index that is stored in this block/leaf
  int GetIndex() const { return value_; }
  virtual std::string toString() const {
    std::stringstream ss;
    ss << value_;
    return ss.str();
  }
  virtual int LargestIndex() const { return value_; }

 protected:
  //! The value of the index for this leaf
  int value_;

 private:
  BlockReorderLeaf()
    : value_(0){};  // hidden from users
};

void tokenize(std::string srcInput, std::string whitespace, std::string prefer, std::vector<std::string>& tokens);

// this function takes a set of tokens and returns the first "block", i.e. those
// values (including) brackets that correspond to the first block
std::vector<std::string>::const_iterator buildSubBlock(
    std::vector<std::string>::const_iterator begin,
    std::vector<std::string>::const_iterator end,
    std::vector<std::string>& subBlock);

// This function takes a tokenized vector and converts it to a block reorder manager
Teuchos::RCP<Xpetra::BlockReorderManager> blockedReorderFromTokens(const std::vector<std::string>& tokens);

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
Teuchos::RCP<const Xpetra::BlockReorderManager> blockedReorderFromString(std::string reorder);
}  // namespace Xpetra

#endif /* XPETRA_BLOCKREORDERMANAGER_HPP_ */
