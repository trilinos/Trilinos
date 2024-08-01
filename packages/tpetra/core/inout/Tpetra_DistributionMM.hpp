// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
