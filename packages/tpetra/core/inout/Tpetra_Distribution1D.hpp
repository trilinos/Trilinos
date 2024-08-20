// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//  1D Row-based Distribution class
//  Assumes square matrix
//  Karen Devine, SNL 
//

#ifndef __TPETRA_DISTRIBUTION1D_HPP
#define __TPETRA_DISTRIBUTION1D_HPP

namespace Tpetra {

// Forward definition
template <typename gno_t, typename scalar_t>
class DistributionLowerTriangularBlock;

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution1D : public Distribution<gno_t,scalar_t> {
// 1D row-wise distribution of matrix and vector entries
// Rows and vector entries may be linearly or randomly distributed, or
// read from a file.
// Row map and vector map are identical

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::nrows;
  using Distribution<gno_t,scalar_t>::Mine;

  Distribution1D(size_t nrows_,
                 const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                 const Teuchos::ParameterList &params) :
                 Distribution<gno_t,scalar_t>(nrows_, comm_, params)
  {
    int npRow = -1;  // Number of processors among which to distribute rows;
                     // Will compute if not set by user
    const Teuchos::ParameterEntry *pe = params.getEntryPtr("nProcessorRows");
    if (pe != NULL) npRow = pe->getValue<int>(&npRow);

    TEUCHOS_TEST_FOR_EXCEPTION(npRow != -1 && npRow != np, std::logic_error,
                               " nProcessorRows " << npRow << " must equal" <<
                               " nProcessors " << np <<
                               " for 1D distribution");

    if (me == 0) std::cout << "\n 1D Distribution: "
                           << "\n     np     = " << np << std::endl;
  }

  // Return whether this rank owns vector entry i.
  virtual bool VecMine(gno_t i) = 0; 

  // Return whether this rank owns nonzero (i,j)
  // Vector map and row map are the same in 1D distribution.
  inline bool Mine(gno_t i, gno_t j) {return VecMine(i);}
  inline bool Mine(gno_t i, gno_t j, int p) {return VecMine(i);}
};

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution1DLinear: public Distribution1D<gno_t,scalar_t> {

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::nrows;

  Distribution1DLinear(size_t nrows_,
                       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                       const Teuchos::ParameterList &params) :
                       Distribution1D<gno_t,scalar_t>(nrows_, comm_, params)
  {
    gno_t nMyRows = getNumRow(me);
    myFirstRow = getFirstRow(me);
    myLastRow = myFirstRow + nMyRows - 1;
  }

  inline enum DistributionType DistType() { return OneDLinear; }

  inline bool VecMine(gno_t i) { return (i >= myFirstRow && i <= myLastRow); }

private:
  gno_t myFirstRow;
  gno_t myLastRow;

  inline size_t getNumRow(int p) { 
    return (nrows / np + (int(nrows % np) > p)); 
  }

  inline gno_t getFirstRow(int p) { 
    return (p * (nrows / np) + std::min<int>(int(nrows % np), p));
  }

// DistributionLowerTriangularBlock class needs a 1DLinear distribution
friend class DistributionLowerTriangularBlock<gno_t,scalar_t>;  
  
};

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution1DRandom : public Distribution1D<gno_t,scalar_t> {


public:
  using Distribution<gno_t,scalar_t>::me;

  Distribution1DRandom(size_t nrows_,
                       const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                       const Teuchos::ParameterList &params) :
                       Distribution1D<gno_t,scalar_t>(nrows_, comm_, params) 
  { if (me == 0) std::cout << "    randomize = true" << std::endl; }

  inline enum DistributionType DistType() { return OneDRandom; }

  inline bool VecMine(gno_t i) { return (this->HashToProc(i) == me); }
};

/////////////////////////////////////////////////////////////////////////////
template <typename gno_t, typename scalar_t>
class Distribution1DVec : public Distribution1D<gno_t,scalar_t> {
// Distribution of nonzeros is determined by the distribution of the
// vector entries, as read from a file.
//
// Assumptions include:
// -  Distribution file containing the vector part assignments (N lines) 
//    is provided.  This file is read during the constructor.  
//    Format for an NxN matrix:
//        line 1 to N:  0-based part assignment of vector entry

public:
  using Distribution<gno_t,scalar_t>::me;
  using Distribution<gno_t,scalar_t>::np;
  using Distribution<gno_t,scalar_t>::comm;
  using Distribution<gno_t,scalar_t>::nrows;
  using Distribution<gno_t,scalar_t>::Mine;

  Distribution1DVec(size_t nrows_,
                    const Teuchos::RCP<const Teuchos::Comm<int> > &comm_, 
                    const Teuchos::ParameterList &params,
                    std::string &distributionfile) :
                    Distribution1D<gno_t,scalar_t>(nrows_, comm_, params)
  {
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
        Teuchos::broadcast<int, int>(*comm, 0, cnt, &(vecpart[start]));
        start += cnt;
        cnt = 0;
      }
    }

    if (me == 0) fpin.close();
  }

  ~Distribution1DVec() {delete [] vecpart;}

  inline enum DistributionType DistType() { return OneDVec; }

  // Vector distribution was read in.
  inline bool VecMine(gno_t i) { return (vecpart[i] == me); }

protected:
  int *vecpart;             // part assignment of vector entries

};

}

#endif
