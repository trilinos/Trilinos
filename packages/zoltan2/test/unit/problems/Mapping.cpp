// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
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
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Test the MappingProblem and MappingSolution classes.
//

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_VectorAdapter.hpp>

//#include <Zoltan2_MappingProblem.hpp>
//#include <Zoltan2_MappingSolution.hpp>

typedef Zoltan2::BasicUserTypes<zscalar_t, zlno_t, zgno_t> zzuser_t;

/////////////////////////////////////////////////////////////////////////////
template <typename User>
class VerySimpleVectorAdapter : public Zoltan2::VectorAdapter<User>
{
public:
  // Defines an (np x nCoordPerRank) grid of points
  // Total number of parts is np * nPartsPerRow; 
  // Half of each row of points belongs to a single part.
  // Lowest part number is lowestPartNum.
  typedef typename Zoltan2::VectorAdapter<User>::lno_t lno_t;
  typedef typename Zoltan2::VectorAdapter<User>::gno_t gno_t;
  typedef typename Zoltan2::VectorAdapter<User>::scalar_t scalar_t;
  typedef typename Zoltan2::VectorAdapter<User>::part_t part_t;

  // Problem dimensions are fixed
  static const int nCoordDim = 2;
  static const int nCoordPerRank = 6;

  // Constructor
  VerySimpleVectorAdapter(
    const Teuchos::Comm<int> &comm_,
    int nPartsPerRow_, 
    int lowestPartNum_, 
    bool useInputParts_=false)
  :
    me(comm_.getRank()),
    np(comm_.getSize()),
    nPartsPerRow(nPartsPerRow_),
    lowestPartNum(lowestPartNum_),
    useInputParts(useInputParts_)
  {
    for (int j = 0; j < nCoordPerRank; j++)
      ids[j] = me * nCoordPerRank + j;

    for (int i = 0, j = 0; i < nPartsPerRow; i++)
      for (int k = 0; k < nCoordPerRank / nPartsPerRow; k++, j++)
        inputparts[j] = lowestPartNum + i + me*nPartsPerRow;

    for (int j = 0; j < nCoordPerRank; j++) {
      coords[0][j] = scalar_t(j);
      coords[1][j] = scalar_t(me);
      if (nCoordDim > 2) coords[2][j] = scalar_t(np);
    }
  }

  // Print the private data
  void printPrivate(std::string hi) {
    // print ids
    std::cout << hi << " private Rank " << me << " ids:    ";
    for (int j = 0; j < nCoordPerRank; j++) std::cout << ids[j] << " ";
    std::cout << std::endl;

    // print coords
    std::cout << hi << " private Rank " << me << " coords: ";
    for (int j = 0; j < nCoordPerRank; j++) {
      std::cout << "(";
      for (int k = 0; k < nCoordDim; k++)
        std::cout << coords[k][j] << ",";
      std::cout << ") ";
    }
    std::cout << std::endl;
 
    // print inputparts
    std::cout << hi << " private Rank " << me << " parts:  ";
    for (int j = 0; j < nCoordPerRank; j++) std::cout << inputparts[j] << " ";
    std::cout << std::endl;
  }
    
  // Print the data as received by the methods; 
  // should match output of printPrivate()
  void printMethods(std::string hi) {
    // print ids
    const gno_t *mids;
    getIDsView(mids);
    std::cout << hi << " methods Rank " << me << " ids:    ";
    for (size_t j = 0; j < getLocalNumIDs(); j++) std::cout << mids[j] << " ";
    std::cout << std::endl;

    // print coords
    const scalar_t **mcoords = new const scalar_t*[getNumEntriesPerID()];
    int *stride = new int[getNumEntriesPerID()];
    for (int k = 0; k < getNumEntriesPerID(); k++)
      getEntriesView(mcoords[k], stride[k], k);

    std::cout << hi << " methods Rank " << me << " coords: ";
    for (size_t j = 0; j < getLocalNumIDs(); j++) {
      std::cout << "(";
      for (int k = 0; k < getNumEntriesPerID(); k++)
        std::cout << mcoords[k][j*stride[k]] << ",";
      std::cout << ") ";
    }
    std::cout << std::endl;
    delete [] mcoords;
    delete [] stride;

    // print inputparts
    const part_t *minputparts;
    getPartsView(minputparts);
    std::cout << hi << " methods Rank " << me << " parts:  ";
    if (minputparts != NULL)
      for (size_t j = 0; j < getLocalNumIDs(); j++)
        std::cout << minputparts[j] << " ";
    else
      std::cout << "not provided";
    std::cout << std::endl;
  }

  // Methods for the VectorAdapter interface
  size_t getLocalNumIDs() const { return nCoordPerRank; }
  void getIDsView(const gno_t *&Ids) const { Ids = ids; }

  int getNumEntriesPerID() const { return nCoordDim; }
  void getEntriesView(const scalar_t *&Coords, int &Stride, int Idx) const
  {
    Coords = coords[Idx];
    Stride = 1;
  }

  void getPartsView(const part_t *&InputPart) const {
    if (useInputParts)
      InputPart = inputparts;
    else
      InputPart = NULL;
  }

private:
  int me;    // my rank
  int np;    // number of ranks
  int nPartsPerRow;   // number of parts created per row of the grid
  int lowestPartNum;  // lowest part number; not necessarily zero
  bool useInputParts; // use the input partition in the tests?  If not,
                      // Zoltan2 uses rank as default part number
  gno_t ids[nCoordPerRank];                  // IDs of this rank's grid points
  part_t inputparts[nCoordPerRank];          // Input part for each local point
  scalar_t coords[nCoordDim][nCoordPerRank]; // Coordinate for each local point
};


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  
  // TEST 1
  {
    // Create a very simple input adapter.
    VerySimpleVectorAdapter<zzuser_t> ia(*comm, 1, 0, true);
    ia.printPrivate("test1");
    ia.printMethods("test1");
  }

  // TEST 2
  {
    // Create a very simple input adapter.
    VerySimpleVectorAdapter<zzuser_t> ia(*comm, 1, 0, false);
    ia.printPrivate("test2");
    ia.printMethods("test2");
  }

  // TEST 3
  {
    // Create a very simple input adapter.
    VerySimpleVectorAdapter<zzuser_t> ia(*comm, 2, 4, true);
    ia.printPrivate("test3");
    ia.printMethods("test3");
  }

  // TEST 4
  {
    // Create a very simple input adapter.
    VerySimpleVectorAdapter<zzuser_t> ia(*comm, 3, 10, true);
    ia.printPrivate("test4");
    ia.printMethods("test4");
  }

  // How do we test a mapping problem and solution?
  // All returned processors must be in range [0,np-1].
  // All returned parts must be in range valid part numbers from input.
  // Exercise thrown error conditions.
  // Correctness of a particular mapping algorithm is beyond the scope of
  // this test.

  std::cout << "PASS" << std::endl;
}
