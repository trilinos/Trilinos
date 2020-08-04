//  Build maps for 1D or 2D matrix distribution
//  Assumes square matrix
//  Karen Devine, SNL 
//

#ifndef __DISTRIBUTION_HPP
#define __DISTRIBUTION_HPP

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <set>
#ifndef __cplusplus
#define __cplusplus
#endif

#include "Teuchos_Comm.hpp"


enum DistributionType {
  TwoDRandom,    // 2D randomly permuted distribution
  TwoDLinear,    // 2D linear distribution
  TwoDVec,       // 2D distribution based on vector assignment in file
  OneDRandom,    // 1D randomly permuted distribution
  OneDLinear,    // 1D linear distribution
  OneDVec,       // 1D distribution based on vector assignment in file
  LowerTriangularBlock, // Seher Acer's lower-triangular block distrib 
                        // for triangle counting
  MMFile         // Use values in matrix-market file as part assignment
};

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution {
public:

  Distribution(gno_t nrows_, 
               const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
               const Teuchos::ParameterList &params) :
               comm(comm_), me(comm_->getRank()), np(comm_->getSize()), 
               nrows(nrows_) { }

  virtual ~Distribution() {};

  // Return the DistributionType for this distribution.
  virtual enum DistributionType DistType() = 0;

  // Return whether this rank owns nonzero (i,j)
  virtual bool Mine(gno_t i, gno_t j) = 0;
  virtual bool Mine(gno_t i, gno_t j, int p) = 0;

  // Return whether this rank owns vector entry i
  virtual bool VecMine(gno_t i) = 0;

  // Map of nonzeros needed for redistribution, handy for other things
  using NZindex_t = std::pair<gno_t, gno_t>;
  struct compareNzIndex {  // sort nonzeros by row, then column
    bool operator() (const NZindex_t &lhs, const NZindex_t &rhs) const
    { if (lhs.first < rhs.first) return true;
      if ((lhs.first == rhs.first) && (lhs.second < rhs.second)) return true;
      return false;
    }
  };

  using LocalNZmap_t = std::map<NZindex_t, scalar_t, compareNzIndex>;

  // Redistribute nonzeros according to the needs of the Distribution
  // Needed only when the final distribution cannot be determined until
  // all nonzeros are known (e.g., final distribution depends on the number 
  // of nonzeros in a row).
  // If the final distribution can be determined before all nonzeros (e.g.,
  // Trilinos' traditional row map), the redistribution routine is a no-op.
  // Thus, the default implementation is a no-op.
  virtual void Redistribute(LocalNZmap_t &localNZ) { };


protected:
  const Teuchos::RCP<const Teuchos::Comm<int> > comm;
  int me;     // my rank
  int np;     // number of ranks
  gno_t nrows;  // global number of rows in the input matrix

  int HashToProc(gno_t i) {
    // TODO  PUT A GOOD HASH FUNCTION HERE!!!
    // TODO  FOR INITIAL TESTING, USE A CYCLIC HASH.  LAME!
    return(i % np);
  }
};

#include "distribution2D.hpp"
#include "distribution1D.hpp"
#include "distributionMM.hpp"
#include "distributionLowerTriangularBlock.hpp"

#endif
