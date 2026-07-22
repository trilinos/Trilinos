// @HEADER
// *****************************************************************************
//             Xpetra: A linear algebra interface package
//
// Copyright 2012 NTESS and the Xpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Xpetra_BlockReorderManager.hpp>

namespace Xpetra {

void BlockReorderManager::SetBlock(int blockIndex, int reorder) {
  TEUCHOS_ASSERT(blockIndex < (int)children_.size());
  Teuchos::RCP<BlockReorderManager> child = Teuchos::rcp(new BlockReorderLeaf(reorder));
  children_[blockIndex]                   = child;
}

void BlockReorderManager::SetBlock(int blockIndex, const Teuchos::RCP<BlockReorderManager>& reorder) {
  TEUCHOS_ASSERT(blockIndex < (int)children_.size());
  children_[blockIndex] = reorder;
}

// this function tokenizes a string, breaking out whitespace but saving the
// brackets [,] as special tokens.
void tokenize(std::string srcInput, std::string whitespace, std::string prefer, std::vector<std::string>& tokens) {
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
std::vector<std::string>::const_iterator buildSubBlock(
    std::vector<std::string>::const_iterator begin,
    std::vector<std::string>::const_iterator end,
    std::vector<std::string>& subBlock) {
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
    if (matched.empty())
      return itr;
  }

  TEUCHOS_ASSERT(matched.empty());

  return itr - 1;
}

// This function takes a tokenized vector and converts it to a block reorder manager
Teuchos::RCP<Xpetra::BlockReorderManager> blockedReorderFromTokens(const std::vector<std::string>& tokens) {
  // base case
  if (tokens.size() == 1)
    return Teuchos::rcp(new Xpetra::BlockReorderLeaf(Teuchos::StrUtils::atoi(tokens[0])));

  // check first and last character
  TEUCHOS_ASSERT(*(tokens.begin()) == "[")
  TEUCHOS_ASSERT(*(tokens.end() - 1) == "]");

  std::vector<Teuchos::RCP<Xpetra::BlockReorderManager> > vecRMgr;
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
  Teuchos::RCP<Xpetra::BlockReorderManager> rMgr = Teuchos::rcp(new Xpetra::BlockReorderManager());
  rMgr->SetNumBlocks(vecRMgr.size());
  for (unsigned int i = 0; i < vecRMgr.size(); i++)
    rMgr->SetBlock(i, vecRMgr[i]);

  return rMgr;
}

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
Teuchos::RCP<const Xpetra::BlockReorderManager> blockedReorderFromString(std::string reorder) {
  // vector of tokens to use
  std::vector<std::string> tokens;

  // manager to be returned

  // build tokens vector
  tokenize(reorder, " \t\n", "[]", tokens);

  // parse recursively and build reorder manager
  Teuchos::RCP<Xpetra::BlockReorderManager> mgr = blockedReorderFromTokens(tokens);

  return mgr;
}

}  // namespace Xpetra
