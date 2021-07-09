// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

//  Build maps for 1D or 2D matrix distribution
//  Assumes square matrix
//  Karen Devine, SNL 
//

#ifndef __TPETRA_DISTRIBUTIONMM_HPP
#define __TPETRA_DISTRIBUTIONMM_HPP

namespace Tpetra 
{

template <typename gno_t, typename scalar_t>
class DistributionMMFile : public Distribution<gno_t,scalar_t> {
// Distribution of nonzeros is determined by nonzero's value 
// as read from Matrix-Market file.
// Vector entry v_i is assigned to the same processor as matrix diagonal a_{ii}.
// For now, we derive the vector entry assignment by accruing information
// about the diagonal entries' assigments during calls to Mine(I,J,V).
// A better scheme might read vector entries' assignments from a file as well,
// or have a separate discovery of vector entries' assignments in the 
// constructor.
// Assumptions include:
// -  If diagonal entries are needed (e.g., for a Laplacian), they are 
//    included in the MMFile
// -  Part assignments are one-based; that is (I,J,V) = (1,1,4) assigns
//    (I,J) to process 3.
// -  Mine(I,J) is undefined; value V must be provided.

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::nrows;

  DistributionMMFile(size_t nrows_,
                     const Teuchos::RCP<const Teuchos::Comm<int> > &comm_,
                     const Teuchos::ParameterList &params) :
                     Distribution<gno_t,scalar_t>(nrows_, comm_, params)
  {
    if (me == 0) std::cout << "\n MMFile Distribution: " 
                           << "\n     np     = " << np << std::endl;
  }

  inline enum DistributionType DistType() { return MMFile; }

  bool Mine(gno_t i, gno_t j) {
    std::cout << "Invalid call to Mine(i,j); " 
         << "MMFile-distribution requires use of Mine(i,j,p) providing "
         << "process assignment p." << std::endl;
    exit(-1);
  }

  bool Mine(gno_t i, gno_t j, int oneBasedRank) {
    // Nonzero (i,j) is Mine if oneBasedRank-1 == me.

    if (oneBasedRank < 1 || oneBasedRank > np) {
      std::cout << "Invalid rank " << oneBasedRank 
           << " provided in user distribution;  "
           << "rank must be in range 1 to " << np << std::endl;
      exit(-1);
    }

    // Keep track of diagonal entries that I own for use in vector map
    if (oneBasedRank-1 == me && i == j) myVecEntries.insert(i);

    return (oneBasedRank-1 == me);
  }

  // myVecEntries keeps track of which diagonal matrix entries are Mine().
  // myVecEntries is not complete until the entire matrix has been viewed
  // by Mine(), so use of VecMine before that point may produce misleading
  // results.
  inline bool VecMine(gno_t i) {
    return (myVecEntries.find(i) != myVecEntries.end());
  }

private:
  std::set<gno_t> myVecEntries;  // vector entries that are assigned to me
};

}
#endif
