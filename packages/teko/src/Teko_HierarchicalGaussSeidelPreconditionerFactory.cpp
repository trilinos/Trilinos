/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#include "Teko_BlockImplicitLinearOp.hpp"
#include "Teko_HierarchicalGaussSeidelPreconditionerFactory.hpp"

using Teuchos::rcp;
using Teuchos::RCP;

namespace Teko {

namespace {

BlockedMultiVector extractSubBlockVector(const BlockedMultiVector& blo,
                                         const std::vector<int>& rows) {
  const int nrows = rows.size();
  TEUCHOS_ASSERT(blockCount(blo) >= nrows);
  std::vector<MultiVector> subVectors;
  for (int row_id = 0; row_id < nrows; ++row_id) {
    const auto row = rows[row_id];
    subVectors.push_back(getBlock(row, blo));
  }

  return buildBlockedMultiVector(subVectors);
}

BlockedLinearOp extractSubblockMatrix(const BlockedLinearOp& blo, const std::vector<int>& rows,
                                      const std::vector<int>& cols) {
  int nOuterBlockRows = blockRowCount(blo);
  for (auto&& row : rows) {
    TEUCHOS_ASSERT(row < nOuterBlockRows);
  }
  int nOuterBlockCols = blockColCount(blo);
  for (auto&& col : cols) {
    TEUCHOS_ASSERT(col < nOuterBlockCols);
  }

  // allocate new operator
  auto subblock = createBlockedOp();

  const int nInnerBlockRows = rows.size();
  const int nInnerBlockCols = cols.size();

  // build new operator
  subblock->beginBlockFill(nInnerBlockRows, nInnerBlockCols);
  for (int innerBlockRow = 0U; innerBlockRow < nInnerBlockRows; ++innerBlockRow) {
    for (int innerBlockCol = 0U; innerBlockCol < nInnerBlockCols; ++innerBlockCol) {
      auto [outerBlockRow, outerBlockCol] =
          std::make_tuple(rows[innerBlockRow], cols[innerBlockCol]);
      auto A_row_col = blo->getBlock(outerBlockRow, outerBlockCol);

      if (A_row_col != Teuchos::null) {
        subblock->setBlock(innerBlockRow, innerBlockCol, A_row_col);
      } else {
        // scan to find first non-null range/domain, construct zero operator
        VectorSpace range;
        VectorSpace domain;

        for (int outerBlock = 0; outerBlock < nOuterBlockCols; ++outerBlock) {
          auto A_ij = blo->getBlock(outerBlockRow, outerBlock);
          if (A_ij != Teuchos::null) {
            range = A_ij->range();
            break;
          }
        }

        for (int outerBlock = 0; outerBlock < nOuterBlockRows; ++outerBlock) {
          auto A_ij = blo->getBlock(outerBlock, outerBlockCol);
          if (A_ij != Teuchos::null) {
            domain = A_ij->domain();
            break;
          }
        }

        TEUCHOS_ASSERT(range != Teuchos::null);
        TEUCHOS_ASSERT(domain != Teuchos::null);

        subblock->setBlock(innerBlockRow, innerBlockCol, zero(range, domain));
      }
    }
  }

  subblock->endBlockFill();

  return subblock;
}

std::vector<std::string> tokenize_input(std::string input, const std::string& delimiter) {
  size_t pos = 0;
  std::vector<std::string> tokens;
  while ((pos = input.find(delimiter)) != std::string::npos) {
    tokens.push_back(input.substr(0, pos));
    input.erase(0, pos + delimiter.length());
  }
  tokens.push_back(input);
  return tokens;
}

// Parse hierarchical block number from:
// <ParameterList name="Hierarchical Block 1">
// </ParameterList>
std::optional<int> extract_hierarchical_block_number(const Teuchos::ParameterList& params,
                                                     const std::string& key) {
  const std::string blockHierarchy = "Hierarchical Block ";
  if (key.find(blockHierarchy) == std::string::npos) return {};
  if (!params.isSublist(key)) return {};

  int blockNumber = -1;
  try {
    auto tokens = tokenize_input(key, " ");
    blockNumber = std::stoi(tokens.back());
  } catch (std::exception& err) {
    std::ostringstream ss;
    ss << "Error occured when trying to parse entry with key " << key << "\n";
    ss << "It said:\n" << err.what() << "\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, ss.str());
  }
  return blockNumber - 1;  // adjust to 0-based indexing
}

std::vector<int> extract_block_numbers(const Teuchos::ParameterList& params) {
  const std::string includedBlocks = "Included Subblocks";

  std::ostringstream ss;
  ss << "Parameter 'Included Subblocks' is missing for params:\n" << params;
  TEUCHOS_TEST_FOR_EXCEPTION(!params.isParameter(includedBlocks), std::runtime_error, ss.str());
  std::vector<int> blocks;
  try {
    auto blockStrs = tokenize_input(params.get<std::string>(includedBlocks), ",");
    for (const auto& blockStr : blockStrs) {
      blocks.emplace_back(std::stoi(blockStr) - 1);  // adjust to 0-based indexing
    }
  } catch (std::exception& err) {
    std::ostringstream errSS;
    errSS << "Error occured when trying to parse 'Included SubBlocks' for params:\n"
          << params << "\n";
    errSS << "It said:\n" << err.what() << "\n";
    TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, errSS.str());
  }
  return blocks;
}

unsigned index(unsigned i, unsigned j, unsigned N) { return i * N + j; }

}  // namespace

NestedBlockGS::NestedBlockGS(const std::map<int, std::vector<int>>& blockToRow_,
                             const std::map<int, LinearOp>& blockToInvOp_, BlockedLinearOp& A_,
                             bool useLowerTriangle_)
    : blockToRow(blockToRow_),
      blockToInvOp(blockToInvOp_),
      A(A_),
      useLowerTriangle(useLowerTriangle_) {
  productRange_    = A->productRange();
  productDomain_   = A->productDomain();
  const auto Nrows = blockToRow.size();
  Ab.resize(Nrows * Nrows);
  for (auto [hierchicalBlockRow, rows] : blockToRow) {
    for (auto [hierchicalBlockCol, cols] : blockToRow) {
      Ab[index(hierchicalBlockRow, hierchicalBlockCol, Nrows)] =
          extractSubblockMatrix(A, rows, cols);
    }
  }
}

void NestedBlockGS::implicitApply(const BlockedMultiVector& r, BlockedMultiVector& z,
                                  const double alpha, const double beta) const {
  const int blocks = blockToRow.size();

  auto src = deepcopy(r);
  std::vector<BlockedMultiVector> rVec;
  std::vector<BlockedMultiVector> zVec;

  for (int b = 0; b < blocks; b++) {
    rVec.push_back(extractSubBlockVector(src, blockToRow.at(b)));
    zVec.push_back(extractSubBlockVector(z, blockToRow.at(b)));
  }

  if (useLowerTriangle) {
    lowerTriangularImplicitApply(rVec, zVec, alpha, beta);
  } else {
    upperTriangularImplicitApply(rVec, zVec, alpha, beta);
  }
}

void NestedBlockGS::upperTriangularImplicitApply(std::vector<BlockedMultiVector>& r,
                                                 std::vector<BlockedMultiVector>& z,
                                                 const double /* alpha */,
                                                 const double /* beta */) const {
  const int blocks = blockToRow.size();

  for (int b = blocks - 1; b >= 0; b--) {
    applyOp(blockToInvOp.at(b), r[b], z[b]);

    for (int i = 0; i < b; i++) {
      auto u_ib = Ab[index(i, b, blocks)];
      if (u_ib != Teuchos::null) {
        applyOp(u_ib, z[b], r[i], -1.0, 1.0);
      }
    }
  }
}

void NestedBlockGS::lowerTriangularImplicitApply(std::vector<BlockedMultiVector>& r,
                                                 std::vector<BlockedMultiVector>& z,
                                                 const double /* alpha */,
                                                 const double /* beta */) const {
  const int blocks = blockToRow.size();

  for (int b = 0; b < blocks; b++) {
    applyOp(blockToInvOp.at(b), r[b], z[b]);

    // loop over each row
    for (int i = b + 1; i < blocks; i++) {
      auto l_ib = Ab[index(i, b, blocks)];
      if (l_ib != Teuchos::null) {
        applyOp(l_ib, z[b], r[i], -1.0, 1.0);
      }
    }
  }
}

HierarchicalGaussSeidelPreconditionerFactory::HierarchicalGaussSeidelPreconditionerFactory() =
    default;

LinearOp HierarchicalGaussSeidelPreconditionerFactory::buildPreconditionerOperator(
    BlockedLinearOp& blo, BlockPreconditionerState& state) const {
  for (auto [hierarchicalBlockNum, blocks] : blockToRow) {
    auto A_bb                          = extractSubblockMatrix(blo, blocks, blocks);
    blockToInvOp[hierarchicalBlockNum] = this->buildBlockInverse(
        *blockToInverse.at(hierarchicalBlockNum), blockToPreconditioner.at(hierarchicalBlockNum),
        A_bb, state, hierarchicalBlockNum);
  }

  return Teuchos::rcp(new NestedBlockGS(blockToRow, blockToInvOp, blo, useLowerTriangle));
}

LinearOp HierarchicalGaussSeidelPreconditionerFactory::buildBlockInverse(
    const InverseFactory& invFact, const Teuchos::RCP<InverseFactory>& precFact,
    const BlockedLinearOp& matrix, BlockPreconditionerState& state,
    int hierarchicalBlockNum) const {
  std::stringstream ss;
  ss << "hierarchical_block_" << hierarchicalBlockNum;

  ModifiableLinearOp& invOp  = state.getModifiableOp(ss.str());
  ModifiableLinearOp& precOp = state.getModifiableOp("prec_" + ss.str());

  if (precFact != Teuchos::null) {
    if (precOp == Teuchos::null) {
      precOp = precFact->buildInverse(matrix);
      state.addModifiableOp("prec_" + ss.str(), precOp);
    } else {
      Teko::rebuildInverse(*precFact, matrix, precOp);
    }
  }

  if (invOp == Teuchos::null)
    if (precOp.is_null())
      invOp = Teko::buildInverse(invFact, matrix);
    else
      invOp = Teko::buildInverse(invFact, matrix, precOp);
  else {
    if (precOp.is_null())
      Teko::rebuildInverse(invFact, matrix, invOp);
    else
      Teko::rebuildInverse(invFact, matrix, precOp, invOp);
  }

  return invOp;
}

void HierarchicalGaussSeidelPreconditionerFactory::initializeFromParameterList(
    const Teuchos::ParameterList& pl) {
  RCP<const InverseLibrary> invLib = getInverseLibrary();

  for (const auto& entry : pl) {
    const auto key               = entry.first;
    const auto value             = entry.second;
    auto hierarchicalBlockNumber = extract_hierarchical_block_number(pl, key);
    if (!hierarchicalBlockNumber) continue;
    const auto& hierarchicalParams       = pl.sublist(key);
    blockToRow[*hierarchicalBlockNumber] = extract_block_numbers(hierarchicalParams);

    std::ostringstream ss;
    ss << "Missing required parameter \"Inverse Type\" for hierarchical block "
       << *hierarchicalBlockNumber << "\n";
    TEUCHOS_TEST_FOR_EXCEPTION(!hierarchicalParams.isParameter("Inverse Type"), std::runtime_error,
                               ss.str());
    auto invStr                              = hierarchicalParams.get<std::string>("Inverse Type");
    blockToInverse[*hierarchicalBlockNumber] = invLib->getInverseFactory(invStr);

    blockToPreconditioner[*hierarchicalBlockNumber] = Teuchos::null;
    if (hierarchicalParams.isParameter("Preconditioner Type")) {
      auto precStr = hierarchicalParams.get<std::string>("Preconditioner Type");
      blockToPreconditioner[*hierarchicalBlockNumber] = invLib->getInverseFactory(precStr);
    }
  }

  useLowerTriangle = false;
  if (pl.isParameter("Use Upper Triangle")) {
    useLowerTriangle = !pl.get<bool>("Use Upper Triangle");
  }
}

}  // namespace Teko
