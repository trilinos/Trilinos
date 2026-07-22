// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//
// Test the MappingProblem and MappingSolution classes.
//

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommHelpers.hpp>

#include <Zoltan2_TestHelpers.hpp>
#include <Zoltan2_VectorAdapter.hpp>

#include <Zoltan2_MappingProblem.hpp>
#include <Zoltan2_MappingSolution.hpp>

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

  // Print the data as received by the methods.
  void print(std::string hi) {
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

  // Return values given to the constructor
  bool adapterUsesInputParts() { return useInputParts; }
  int adapterNPartsPerRow() { return nPartsPerRow; }
  int adapterLowestPartNum() { return lowestPartNum; }

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

//////////////////////////////////////////////////////////////////////////////
// Validate a mapping solution obtained without a partitioning solution
template <typename Adapter>
bool validMappingSolution(
  Zoltan2::MappingSolution<Adapter> &msoln,
  Adapter &ia,
  const Teuchos::Comm<int> &comm
)
{
  // Correctness of a particular mapping algorithm is beyond the scope of
  // this test.
  typedef typename Adapter::part_t part_t;

  bool aok = true;

  // All returned processors must be in range [0,np-1].

  if (ia.adapterUsesInputParts()) {
    // Adapter provided input parts
    const part_t *inputParts; 
    ia.getPartsView(inputParts);

    aok = validMappingSolution(msoln, ia, inputParts, comm);
  }
  else {
    // Adapter did not provide input parts; use default part = me for all IDs
    int me = comm.getRank();

    part_t *defaultParts = new part_t[ia.getLocalNumIDs()];
    for (size_t i = 0; i < ia.getLocalNumIDs(); i++)
      defaultParts[i] = me;

    aok = validMappingSolution(msoln, ia, defaultParts, comm);

    delete [] defaultParts;
  }

  return aok;
}

//////////////////////////////////////////////////////////////////////////////
// Validate a mapping solution obtained with a partitioning solution
template <typename Adapter>
bool validMappingSolution(
  Zoltan2::MappingSolution<Adapter> &msoln,
  Adapter &ia,
  Zoltan2::PartitioningSolution<Adapter> &psoln,
  const Teuchos::Comm<int> &comm
)
{
  typedef typename Adapter::part_t part_t;

  const part_t *assignedParts = psoln.getPartListView();

  return(validMappingSolution(msoln, ia, assignedParts, comm));
}

//////////////////////////////////////////////////////////////////////////////
// Validate a mapping solution.
// This test checks only whether the mapping solution is valid.  Correctness
// of a particular mapping algorithm is beyond the scope of this test.
template <typename Adapter>
bool validMappingSolution(
  Zoltan2::MappingSolution<Adapter> &msoln,
  Adapter &ia,
  const typename Adapter::part_t *useTheseParts, 
                                 // Parts to check for correct mapping to ranks;
                                 // may be from Adapter, 
                                 // from PartitioningSolution, or 
                                 // default to current rank
  const Teuchos::Comm<int> &comm
)
{
  typedef typename Adapter::part_t part_t;

  int np = comm.getSize();

  bool aok = true;

  // Min and max part numbers in useTheseParts
  part_t globalmin, localmin = std::numeric_limits<part_t>::max();
  part_t globalmax, localmax = 0;

  for (size_t i = 0; i < ia.getLocalNumIDs(); i++) {

    // All returned processors must be in range [0,np-1].
    int r = msoln.getRankForPart(useTheseParts[i]);
    if (r < 0 || r >= np) {
      aok = false;
      std::cout << __FILE__ << ":" << __LINE__ << " "
                << "Invalid rank " << r << " of " << np << " returned" 
                << std::endl;
    }

    // Find local max/min part number
    part_t p = useTheseParts[i];
    if (p > localmax) localmax = p;
    if (p < localmin) localmin = p;
  }

  // Find global max/min part number
  Teuchos::reduceAll<int, part_t>(comm, Teuchos::REDUCE_MAX, 1, 
                                  &localmax, &globalmax);
  Teuchos::reduceAll<int, part_t>(comm, Teuchos::REDUCE_MIN, 1, 
                                  &localmin, &globalmin);

  part_t *parts;
  part_t nParts;
  msoln.getMyPartsView(nParts, parts);

  for (part_t i = 0; i < nParts; i++) {

    // All returned parts must at least be in the range of part numbers 
    // from useTheseParts
    part_t p = parts[i];
    if ((p < globalmin) || (p > globalmax)) {
      aok = false;
      std::cout << __FILE__ << ":" << __LINE__ << " "
                << "Invalid part " << p << " of " << np << " returned" 
                << std::endl;
    }
  }

  // Test the error checking in mapping solution;
  // each test should throw an error.

  part_t ret;
  bool errorThrownCorrectly = false;
  part_t sillyPart = globalmax+10;
  try {
    ret = msoln.getRankForPart(sillyPart);
  }
  catch (std::exception &e) {
    errorThrownCorrectly = true;
  }
  if (errorThrownCorrectly == false) {
    aok = false;
    std::cout << __FILE__ << ":" << __LINE__ << " "
              << "Mapping Solution accepted a too-high part number "
              << sillyPart << " returned " << ret << std::endl;
  }

  errorThrownCorrectly = false;
  sillyPart = globalmin - 1;
  try {
    ret = msoln.getRankForPart(sillyPart);
  }
  catch (std::exception &e) {
    errorThrownCorrectly = true;
  }
  if (errorThrownCorrectly == false) {
    aok = false;
    std::cout << __FILE__ << ":" << __LINE__ << " "
              << "Mapping Solution accepted a too-low part number "
              << sillyPart << " returned " << ret << std::endl;
  }

  return aok;
}

//////////////////////////////////////////////////////////////////////////////

template <typename Adapter>
bool runTest(
  Adapter &ia,
  const RCP<const Teuchos::Comm<int> > &comm,
  std::string hi
)
{
  typedef typename Adapter::part_t part_t;
  typedef typename Adapter::scalar_t scalar_t;

  int me = comm->getRank();
  int np = comm->getSize();

  bool allgood = true;

  Teuchos::ParameterList params;
  params.set("mapping_algorithm", "block");

  // Test mapping using default machine
  if (me == 0)
    std::cout << "Testing Mapping using default machine" << std::endl;

  Zoltan2::MappingProblem<Adapter> mprob1(&ia, &params, comm);
  mprob1.solve();

  Zoltan2::MappingSolution<Adapter> *msoln1 = mprob1.getSolution();

  if (!validMappingSolution<Adapter>(*msoln1, ia, *comm)) {
    allgood = false;
    if (me == 0) 
      std::cout << hi << " FAILED: invalid mapping solution" << std::endl;
  }


  // Test mapping explicitly using default machine
  typedef Zoltan2::MachineRepresentation<scalar_t, part_t> machine_t;
  machine_t defMachine(*comm);

#ifdef KDD
  if (me == 0)
    std::cout << "Testing Mapping using explicit machine" << std::endl;

  Zoltan2::MappingProblem<Adapter, machine_t> mprob2(&ia, &params, comm,
                                                     NULL, &defMachine);
  mprob2.solve();

  Zoltan2::MappingSolution<Adapter> *msoln2 = mprob2.getSolution();
 
  if (!sameMappingSolution(*msoln1, *msoln2, *comm)) {
    allgood = false;
    if (me == 0) 
      std::cout << hi << " FAILED: solution with explicit machine "
                   "differs from default" << std::endl;
  }
#endif
    
  // Test mapping with a partitioning solution
  if (me == 0)
    std::cout << "Testing Mapping using a partitioning solution" << std::endl;

  RCP<const Zoltan2::Environment> env = rcp(new Zoltan2::Environment(comm));
  Zoltan2::PartitioningSolution<Adapter> psoln(env, comm, 0);

  ArrayRCP<part_t> partList(ia.getLocalNumIDs());
  for (size_t i = 0; i < ia.getLocalNumIDs(); i++)
    partList[i] = (me + 1) % np;  

  psoln.setParts(partList);

#ifdef HAVE_ZOLTAN2_MPI
  // Use an MPI_Comm, just to exercise that bit of code
  // In real life, no one should extract the MPI_Comm from the Teuchos::Comm;
  // he should use the Teuchos::Comm.  But for testing,
  // we need to exercise the MPI_Comm interface.
  MPI_Comm mpicomm =  Teuchos::getRawMpiComm(*comm);
  Zoltan2::MappingProblem<Adapter, machine_t> mprob3(&ia, &params, mpicomm,
                                                     NULL, &defMachine);
#else
  Zoltan2::MappingProblem<Adapter, machine_t> mprob3(&ia, &params, comm,
                                                     NULL, &defMachine);
#endif

  mprob3.solve();

  Zoltan2::MappingSolution<Adapter> *msoln3 = mprob3.getSolution();
 
  if (!validMappingSolution(*msoln3, ia, psoln, *comm)) {
    allgood = false;
    if (me == 0) 
      std::cout << hi << " FAILED: invalid mapping solution "
                   "from partitioning solution" << std::endl;
  }
  return allgood;
}

//////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int me = comm->getRank();
  bool allgood = true;

  typedef VerySimpleVectorAdapter<zzuser_t> vecAdapter_t;
  //typedef vecAdapter_t::part_t part_t;
  //typedef vecAdapter_t::scalar_t scalar_t;

  // TEST 1
  {
    int nPartsPerRow = 1;
    int firstPart = 0;
    bool useInputParts = true;

    vecAdapter_t ia(*comm, nPartsPerRow, firstPart, useInputParts);
    ia.print("test1");

    allgood = runTest(ia, comm, "test1");
  }

#ifdef KDD
  // TEST 2
  {
    // Create a very simple input adapter.
    int nPartsPerRow = 1;
    int firstPart = 0;
    bool useInputParts = false;
    vecAdapter_t ia(*comm, nPartsPerRow, firstPart, useInputParts);
    ia.print("test2");
  }

  // TEST 3
  {
    // Create a very simple input adapter.
    int nPartsPerRow = 2;
    int firstPart = 4;
    bool useInputParts = true;
    vecAdapter_t ia(*comm, nPartsPerRow, firstPart, useInputParts);
    ia.print("test3");
  }

  // TEST 4
  {
    // Create a very simple input adapter.
    int nPartsPerRow = 3;
    int firstPart = 10;
    bool useInputParts = true;
    vecAdapter_t ia(*comm, nPartsPerRow, firstPart, useInputParts);
    ia.print("test4");
  }
#endif

  if (allgood && (me == 0))
    std::cout << "PASS" << std::endl;
}
