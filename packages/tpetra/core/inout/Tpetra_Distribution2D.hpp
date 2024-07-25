// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//  2D matrix distribution
//  Assumes square matrix
//  Karen Devine, SNL 
//

#ifndef __TPETRA_DISTRIBUTION2D_HPP
#define __TPETRA_DISTRIBUTION2D_HPP

namespace Tpetra
{

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution2D : public Distribution<gno_t,scalar_t> {
// Processors will be laid out logically first down columns then 
// across rows.  For example, assume np = 8, npRows=2, npCols=4.
// Then the processors will be laid out in 2D as
//   0  2  4  6
//   1  3  5  7
//
// The matrix will be distributed using np=8 row blocks:
//   0  2  4  6
//   1  3  5  7
//   0  2  4  6
//   1  3  5  7
//   0  2  4  6
//   1  3  5  7
//   0  2  4  6
//   1  3  5  7
// 
// The vector will be distributed linearly or randomly.  The row and 
// column maps will be built to allow only row- or column-based 
// communication in the matvec.

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::comm;
  using Distribution<gno_t,scalar_t>::nrows;
  using Distribution<gno_t,scalar_t>::Mine;

  Distribution2D(size_t nrows_, 
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                 const Teuchos::ParameterList &params) :
                 Distribution<gno_t,scalar_t>(nrows_, comm_, params),
                 npRows(-1), npCols(-1)
  {

    {
      const Teuchos::ParameterEntry *pe = params.getEntryPtr("nProcessorRows");
      if (pe != NULL)
        npRows = pe->getValue<int>(&npRows);
    }

    {
      const Teuchos::ParameterEntry *pe = params.getEntryPtr("nProcessorCols");
      if (pe != NULL)
        npCols = pe->getValue<int>(&npCols);
    }

    // Compute the processor configuration npRows * npCols

    if (npRows == -1 && npCols == -1) { // Compute both npRows and npCols 
      // First guess
      npRows = (int)(sqrt(np));   
      npCols = np / npRows;
      // Adjust npRows so that npRows * npCols == np
      while (npRows * npCols != np) {
        npRows++;
        npCols = np / npRows;
      }
    }
    else {  // User specified either npRows or npCols
      if (npRows == -1) // npCols specified; compute npRows
        npRows = np / npCols;
      else if (npCols == -1) // npRows specified; compute npCols
        npCols = np / npRows;

      if (npCols * npRows != np) {
        TEUCHOS_TEST_FOR_EXCEPTION(npRows * npCols != np, std::logic_error,
                                   "nProcessorCols " << npCols <<
                                   " * nProcessorRows " << npRows <<
                                   " = " << npCols * npRows <<
                                   " must equal nProcessors " << np <<
                                   " for 2D distribution");
      }
    }
    if (me == 0) 
      std::cout << "\n2D Distribution: " 
                << "\n    npRows = " << npRows
                << "\n    npCols = " << npCols
                << "\n    np     = " << np << std::endl;

    mypCol = this->TWODPCOL(me);
    mypRow = this->TWODPROW(me);
  }

  virtual ~Distribution2D() {};

protected:

  // Return the processor row for rank p
  inline int TWODPROW(int p) {return (p % npRows);}

  // Return the processor column for rank p
  inline int TWODPCOL(int p) {return (p / npRows);}

  // Return the rank for processor row i and processor column j
  inline int TWODPRANK(int i, int j) {return (j * npRows + (j+i) % npRows);}

  int npRows;  // Number of processor rows
  int npCols;  // Number of processor columns
  int mypRow;  // This rank's processor row
  int mypCol;  // This rank's processor column
};

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution2DLinear : public Distribution2D<gno_t,scalar_t> {

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::nrows;
  using Distribution2D<gno_t,scalar_t>::npRows;
  using Distribution2D<gno_t,scalar_t>::npCols;
  using Distribution2D<gno_t,scalar_t>::mypRow;
  using Distribution2D<gno_t,scalar_t>::mypCol;

  Distribution2DLinear(size_t nrows_,
                       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                       const Teuchos::ParameterList &params) :
                       Distribution2D<gno_t,scalar_t>(nrows_, comm_, params)
  {
    // Build vector describing distribution of vector entries across ranks
    entries.assign(np+1, nrows / np);
    int nExtraEntries = nrows % np;

    // Distribute the extra entries evenly among processors.
    // To evenly distribute them extra entries among processor rows and 
    // columns, we distribute them along diagonals of the matrix distribution.
    // For our example, assume we have seven extra values (the max possible
    // with np=8).  Then we give one extra entry each to ranks 
    // [0, 3, 4, 7, 1, 2, 5].  For fewer extra entries, we follow the same
    // order of assignment, and just stop early.
  
    for (int cnt = 0, i = 0; (cnt < nExtraEntries) && (i < npRows); i++) {
      for (int j = 0; (cnt < nExtraEntries) && (j < npCols); cnt++, j++) {
        int rankForExtra = Distribution2D<gno_t,scalar_t>::TWODPRANK(i, j); 
        entries[rankForExtra+1]++;  // Store in rankForExtra+1 to simplify 
                                    // prefix sum.
      }
    }

    // Perform prefix sum of entries.
    entries[0] = 0;
    for (int i = 1; i <= np; i++)
      entries[i] = entries[i-1] + entries[i];
    // Now entries contains the first vector entry for each rank.

    // Column map:  Easy; consecutive entries for all ranks in column.
    int firstRank = mypCol * npRows;  // First rank in my column
    myFirstCol = entries[firstRank];

    gno_t nMyCols = 0;
    for (int i = firstRank; i < firstRank + npRows; i++)
      nMyCols += entries[i+1] - entries[i];
    myLastCol = myFirstCol + nMyCols - 1;
  }

  inline enum DistributionType DistType() { return TwoDLinear; }

  bool Mine(gno_t i, gno_t j) { 
    int idx = int(float(i) * float(np) / float(nrows));
    while (i < entries[idx]) idx--;
    while (i >= entries[idx+1]) idx++;
    return ((mypRow == Distribution2D<gno_t,scalar_t>::TWODPROW(idx)) // RowMine
            && (j >= myFirstCol && j <= myLastCol));  // ColMine
  }
  inline bool Mine(gno_t i, gno_t j, int p) {return Mine(i,j);}

  inline bool VecMine(gno_t i) {
    return(i >= entries[me] && i < entries[me+1]); 
  }


private:
  std::vector<gno_t> entries; // Describes vector entries' distribution to ranks
                              // Organized like vtxdist
  gno_t myFirstCol;       // First column owned by this rank
  gno_t myLastCol;        // Last column owned by this rank
};

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution2DRandom : public Distribution2D<gno_t,scalar_t> {

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution2D<gno_t,scalar_t>::mypRow;
  using Distribution2D<gno_t,scalar_t>::mypCol;

  Distribution2DRandom(size_t nrows_,
                       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                       const Teuchos::ParameterList &params) :
                       Distribution2D<gno_t,scalar_t>(nrows_, comm_, params)
  { if (me == 0) std::cout << "    randomize = true" << std::endl; }

  inline enum DistributionType DistType() { return TwoDRandom; }

  inline bool Mine(gno_t i, gno_t j) { 
    return ((mypRow == this->TWODPROW(this->HashToProc(i))) &&  // RowMine
            (mypCol == this->TWODPCOL(this->HashToProc(j))));   // ColMine
  }
  inline bool Mine(gno_t i, gno_t j, int p) {return Mine(i,j);}

  inline bool VecMine(gno_t i) { return (me == this->HashToProc(i)); }

};

///////////////////////////////////////////////////////////////////////////////

template <typename gno_t, typename scalar_t>
class Distribution2DVec : public Distribution2D<gno_t,scalar_t>
{
// Distribute non-zeros in a 2D manner based on the vector distribution
// and the nprows x npcols configuration;
// rows are assigned to same process owning the vector entry.
public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::comm;
  using Distribution<gno_t,scalar_t>::nrows;
  using Distribution2D<gno_t,scalar_t>::npRows;
  using Distribution2D<gno_t,scalar_t>::npCols;

  Distribution2DVec(size_t nrows_, 
                    const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                    const Teuchos::ParameterList &params,
                    std::string &distributionfile) :
                    Distribution2D<gno_t,scalar_t>(nrows_, comm_, params)
  {
    if (me == 0) std::cout << "\n 2DVec Distribution: "
                           << "\n      np     = " << np << std::endl;
    std::ifstream fpin;
    if (me == 0) {
      fpin.open(distributionfile.c_str(), std::ios::in);
      if (!fpin.is_open()) {
        std::cout << "Error:  distributionfile " << distributionfile
             << " not found" << std::endl;
        exit(-1);
      }
    }

    // Read the vector part assignment and broadcast it to all processes.
    // Broadcast in chunks of bcastsize values.
    // TODO:  Make the vector part assignment more scalable instead of 
    // TODO:  storing entire vector on every process.

    vecpart = new int[nrows];

    const int bcastsize = 1000000;

    gno_t start = 0;
    int cnt = 0;
    for (size_t i = 0; i < nrows; i++) {
      if (me == 0) fpin >> vecpart[i];
      cnt++;
      if (cnt == bcastsize || i == nrows-1) {
        Teuchos::broadcast(*comm, 0, cnt, &(vecpart[start]));
        start += cnt;
        cnt = 0;
      }
    }

    if (me == 0) fpin.close();
  }

  ~Distribution2DVec() {delete [] vecpart;}

  inline enum DistributionType DistType() { return TwoDVec; }

  bool Mine(gno_t i, gno_t j) {
    return (me == (vecpart[i] % npRows + (vecpart[j] / npRows) * npRows));
  }
  inline bool Mine(gno_t i, gno_t j, int p) {return Mine(i,j);}

  inline bool VecMine(gno_t i) { return(vecpart[i] == me); }

private:
  int *vecpart;

};

}
#endif
