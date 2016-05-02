// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef XPETRA_BLOCKREORDERMANAGER_HPP_
#define XPETRA_BLOCKREORDERMANAGER_HPP_

namespace Xpetra {
class BlockReorderManager {
public:
  //! @name Constructors

  //! Basic empty constructor
  BlockReorderManager() : children_(0) {}

  //! Copy constructor
  BlockReorderManager(const BlockReorderManager & bmm) :
    children_(bmm.children_.size()) {
    for(size_t i = 0; i < children_.size(); i++) children_[i] = bmm.children_[i]->Copy();
  }

  //! empty destructor
  virtual ~BlockReorderManager() { }

  //@}

  //! returns copy of this object
  virtual Teuchos::RCP<BlockReorderManager> Copy() const {
    return Teuchos::rcp(new BlockReorderManager(*this));
  }

  //! Sets number of subblocks
  virtual void SetNumBlocks( size_t sz ) {
    children_.clear(); children_.resize(sz);
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
  virtual void SetBlock(int blockIndex, int reorder);/* {
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
  virtual void SetBlock(int blockIndex, const Teuchos::RCP<BlockReorderManager> & reorder);/* {
    TEUCHOS_ASSERT(blockIndex < (int) children_.size());
    children_[blockIndex] = reorder;
  }*/

  //! \brief Get a particular block. If there is no block at this index location return a new one
  virtual const Teuchos::RCP<BlockReorderManager> GetBlock(int blockIndex) {
    TEUCHOS_ASSERT(blockIndex<(int) children_.size());

    if(children_[blockIndex]==Teuchos::null)
      children_[blockIndex] = Teuchos::rcp(new BlockReorderManager());
    return children_[blockIndex];
  }

  //! \brief Get a particular block. If there is no block at this index location return a new one
  virtual const Teuchos::RCP<const BlockReorderManager> GetBlock(int blockIndex) const {
    TEUCHOS_ASSERT(blockIndex<(int) children_.size());
    return children_[blockIndex];
  }

  //! for sanities sake, print a readable string
  virtual std::string toString() const {
    // build the string by recursively calling each child
    std::stringstream ss;
    ss << "[";
    for(size_t i = 0; i < children_.size(); i++) {
      if(children_[i] == Teuchos::null)
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
    for(size_t i = 0; i<children_.size(); i++) {
      if(children_[i]!=Teuchos::null) {
        int subMax = children_[i]->LargestIndex();
        max = max > subMax ? max : subMax;
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
  BlockReorderLeaf(int ind) : value_(ind) {}
  BlockReorderLeaf(const BlockReorderLeaf & brl) : value_(brl.value_) {}

  virtual Teuchos::RCP<BlockReorderManager> Copy() const {
    return Teuchos::rcp(new BlockReorderLeaf(*this));
  }

  virtual size_t GetNumBlocks() const { return 0; }
  virtual void SetNumBlocks(size_t sz) {}
  virtual void SetBlock(int blockIndex, int reorder) { }
  virtual void SetBlock(int blockIndex, const Teuchos::RCP<BlockReorderManager> & reorder) {}
  virtual const Teuchos::RCP<BlockReorderManager> GetBlock(int blockIndex) {
    return Teuchos::null;
  }
  virtual const Teuchos::RCP<const BlockReorderManager> GetBlock(int blockIndex)const {
    return Teuchos::null;
  }
  //! Get the index that is stored in this block/leaf
  int GetIndex() const { return value_; }
  virtual std::string toString() const {
    std::stringstream ss; ss << value_; return ss.str();
  }
  virtual int LargestIndex() const { return value_; }
protected:
  //! The value of the index for this leaf
  int value_;

private:
  BlockReorderLeaf() : value_(0) {}; // hidden from users
};

void BlockReorderManager::SetBlock(int blockIndex, int reorder) {
  TEUCHOS_ASSERT(blockIndex < (int) children_.size());
  Teuchos::RCP<BlockReorderManager> child = Teuchos::rcp(new BlockReorderLeaf(reorder));
  children_[blockIndex] = child;
}

void BlockReorderManager::SetBlock(int blockIndex, const Teuchos::RCP<BlockReorderManager> & reorder) {
  TEUCHOS_ASSERT(blockIndex < (int) children_.size());
  children_[blockIndex] = reorder;
}


}

#endif /* XPETRA_BLOCKREORDERMANAGER_HPP_ */
