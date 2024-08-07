// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Program that demonstrates how to emulate distributed directories that
// can sum-up entries using Tpetra classes (Map, Vector).
// Similar functionality (but, of course, without the sum) using a Zoltan_DD
// is also written.
// It would be fair to compare the runtimes of the Tpetra version with the
// Zoltan version, even without the sum operation.  And if/when we implement
// the sum operation, we could add a test of it to this program.

#include <Teuchos_RCP.hpp>
#include <Teuchos_ArrayView.hpp>

#include <unordered_map>

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>

#include <Zoltan2_TPLTraits.hpp>

#include <zoltan_dd_cpp.h>
#include <unordered_set>


static const size_t TOOMANY = 100;

//////////////////////////////////////////////////////////////////////////////
// Functions for debugging
template <typename T>
void printTpetraThing(T &thing, const std::string &msg, 
                 std::ostream &ostr = std::cout)
{
  ostr << msg << std::endl;
  Teuchos::FancyOStream fancy(Teuchos::rcpFromRef(ostr));
  thing.describe(fancy, Teuchos::VERB_EXTREME);
}

template <typename vector_t>
void printVector(vector_t &vec, const std::string &msg,
                 std::ostream &ostr = std::cout)
{
  if (vec.getGlobalLength() > TOOMANY) return;
  printTpetraThing<vector_t>(vec, msg, ostr);
}

template <typename map_t>
void printMap(map_t &map, const std::string &msg,
              std::ostream &ostr = std::cout)
{
  if (map.getGlobalNumElements() > TOOMANY) return;
  printTpetraThing<map_t>(map, msg, ostr);
}

//////////////////////////////////////////////////////////////////////////////
// Class to generate input IDs based on input values
template <typename id_t>
class IDs {
public:

  IDs(size_t nIds_, float fracShared_, id_t idBase_, int idStride_,
      Teuchos::RCP<const Teuchos::Comm<int> > &comm_) :
      nIds(nIds_),
      nShared(std::ceil(nIds * fracShared_)), 
      nUnique(nIds - nShared),
      idBase(idBase_), 
      idStride(idStride_),
      contiguous(idStride_ == 1),
      ids(nIds, 0), comm(comm_)
  {
    int me = comm->getRank();
    int np = comm->getSize();

    // Generate the uniquely owned IDs; 
    size_t cnt = 0;
    for (size_t i = 0; i < nUnique; i++) 
      ids[cnt++] = idBase + ((me * nUnique + i) * idStride);

    // Generate the shared IDs; they are the highest-numbered IDs
    size_t gnUnique = nUnique * np;
    std::srand(me);

    for (size_t i = 0; i < nShared; i++) {
      size_t r = rand() % nShared;
      ids[cnt++] = idBase + ((gnUnique + r) * idStride);
    }

    print();
  }

  void print()
  { 
    if (nIds > TOOMANY) return;

    int me = comm->getRank();

    std::cout << me << " nIds = " << nIds << "; nUnique = " << nUnique
              << "; nShared = " << nShared << std::endl;

    std::cout << me << " Unique: ";
    for (size_t i = 0; i < nUnique; i++) std::cout << ids[i] << " ";
    std::cout << std::endl;

    std::cout << me << " Shared: ";
    for (size_t i = 0; i < nShared; i++) std::cout << ids[nUnique+i] << " ";
    std::cout << std::endl;
  }

  bool TpetraDDTest();

  bool ZoltanDDTest();

private:

  typedef int scalar_t;    // use int since we are counting occurrences


  size_t nIds;             // Number of IDs per processor
  size_t nShared;          // Number of shared IDs per processor
  size_t nUnique;          // Number of unique IDs per processor

  id_t idBase;             // Smallest possible ID
  int idStride;            // Offset between IDs; 1 provides contiguous 
  bool contiguous;         // Flag indicating whether IDs are contiguous

  std::vector<id_t> ids;   // Ids generated on this proc

  Teuchos::RCP<const Teuchos::Comm<int> > comm;
};


//////////////////////////////////////////////////////////////////////////////
// Test of DDirectory-like functionality using Tpetra
// Basic steps:
// T1:  Create an overlapped map M1 and a vector V1 that uses it; 
//      V1[i] = number of local occurrencts of ID i
// T2:  Call createOneToOne to create a one-to-one map M2; 
//      create a vector M2 that uses it.
// T3:  Create an Import object between the overlapped map and 
//      the one-to-one map.
// T4:  Import with Tpetra::ADD from V1 to V2 to count the occurrences
// T5:  Import with Tpetra::REPLACE from V2 to V1 to return the 
//      result to the original processors


template <typename id_t>
bool IDs<id_t>::TpetraDDTest()
{
  typedef typename Tpetra::Map<int, id_t> map_t;
  typedef typename Teuchos::RCP<const map_t> rcpmap_t;

  typedef typename Tpetra::Vector<scalar_t, int, id_t> vector_t;
  typedef typename Teuchos::ArrayRCP<scalar_t> vectordata_t;

  // Step T1

  // Find unique IDs numbers and indexBase on this processor for 
  // constructing overlapping Tpetra::Map for IDs.

  id_t minId = std::numeric_limits<id_t>::max();

  std::unordered_map<id_t, size_t> uniqueIds;
  uniqueIds.reserve(nIds);  // Worst case

  for (size_t i = 0; i < nIds; i++) {
    id_t id = ids[i];
    if (id < minId) minId = id;
    if (uniqueIds.find(id) == uniqueIds.end())
      uniqueIds[id] = 1;
    else
      uniqueIds[id]++;
  }

  // Compute global indexBase

  id_t indexBase;
  Teuchos::reduceAll(*comm, Teuchos::REDUCE_MIN, 1, &minId, &indexBase);

  // Create list of locally unique IDs to use in Tpetra::Map creation

  size_t nUniqueIds = uniqueIds.size();
  Teuchos::Array<id_t> uniqueIdsList(nUniqueIds);

  size_t cnt = 0;
  for (auto it = uniqueIds.begin(); it != uniqueIds.end(); it++)
    uniqueIdsList[cnt++] = it->first;

  // Build Tpetra::Map for the given local ids

  size_t dummy = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  rcpmap_t idMap = Teuchos::rcp(new map_t(dummy, uniqueIdsList(), indexBase,
                                           comm), true);

  // Create Vector using this map.
  // This vector will store number of occurrences of each id across procs

  vector_t idVec(idMap, 0.);

  // Use the multivector for counting number of occurrences of each id

  {
    vectordata_t idData = idVec.getDataNonConst();
    for (auto it = uniqueIds.begin(); it != uniqueIds.end(); it++) {
      id_t idx = idMap->getLocalElement(it->first);
      idData[idx] = it->second;
    }
  }

  printMap(*idMap, "idMap ");

  printVector(idVec, "idVec before ");

  // Step T2

  // Create a one-to-one map and a vector that uses it

  rcpmap_t oto_idMap;
  if (contiguous) {
    // For contigous ids, can use Tpetra default Map
    id_t gnUnique = idMap->getMaxAllGlobalIndex() - indexBase + 1;
    oto_idMap = Teuchos::rcp(new map_t(gnUnique, indexBase, comm));
  }
  else {
    // Since ids may not be contiguous, cannot use default Tpetra Map
    oto_idMap = Tpetra::createOneToOne(idMap);
  }

  vector_t oto_idVec(oto_idMap);

  // Step T3

  // Create an exporter between the two maps

  typedef Tpetra::Export<int, id_t> export_t;

  export_t idExporter(idMap, oto_idMap);

  printTpetraThing(idExporter, "exporter ");

  // Step T4

  // Compute the number of occurrences of shared ids
  printMap(*oto_idMap, "oto_idMap ");

  printVector(oto_idVec, "oto_idVec BEFORE EXPORT");

  oto_idVec.doExport(idVec, idExporter, Tpetra::ADD);

  printVector(oto_idVec, "oto_idVec AFTER EXPORT");

  // Step T5

  // Send the result back to orig processors
  // so that we can identify shared IDs in orig procs

  idVec.doImport(oto_idVec, idExporter, Tpetra::REPLACE);

  printVector(idVec, "idVec after ");

  // Check the result
  size_t cntShared = 0;
  {
    auto idData = idVec.getDataNonConst();
    for (size_t i = 0; i < idVec.getLocalLength(); i++)
      if (idData[i] > 1) cntShared++;
  }

  std::cout << comm->getRank() << " cntShared = " << cntShared
                               << "; nShared = " << nShared << std::endl;

  return (cntShared == nShared);
}

//////////////////////////////////////////////////////////////////////////////
// Perform nearly same operations using Zoltan DD.
// Note that Zoltan cannot current do Tpetra::ADD; it will do Tpetra::REPLACE.
// (That is why we are doing this study!)
// Thus, the computed results will not be correct, but except for the addition,
// the mechanics are all the same.

template <typename id_t>
bool IDs<id_t>::ZoltanDDTest()
{

  if (nIds > size_t(std::numeric_limits<int>::max()))
    throw std::runtime_error("Problem too large for Zoltan_DD");

  int inIds = int(nIds);

  // Build a Zoltan directory for the IDs; don't care about local IDs

#ifdef HAVE_MPI
  MPI_Comm mpicomm = Teuchos::getRawMpiComm(*comm);
#else
  {  // Use siMPI from Zoltan here; make sure it is initialized
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
      int narg = 0;
      char **argv = NULL;
      MPI_Init(&narg, &argv);
    }
  }
  MPI_Comm mpicomm = MPI_COMM_WORLD;
#endif

  int nIdEnt = Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, id_t>::NUM_ID;

  Zoltan_DD zz(mpicomm, nIdEnt, 0, sizeof(scalar_t), inIds, 0);

  // TODO:  Can think about appropriate hash functions here.
  // To match Tpetra's default map, we could use something like DD_Hash_Fn1.

  // Allocate space for user data.  
  // In this case, the user data is occurrence counts that we would like to
  // add up in the directory.
  std::vector<scalar_t> user(nIds, 1.);
  
  // Depending on size of id_t, may need to copy IDs to array of ZOLTAN_ID_TYPE
  // TODO:  Add TPL_Traits that don't require ArrayView.  Doh.

  ZOLTAN_ID_PTR zgids = NULL;
  if (nIds) {
    Teuchos::ArrayView<id_t> av(&(ids[0]), nIds);
    Zoltan2::TPL_Traits<ZOLTAN_ID_PTR, id_t>::ASSIGN_ARRAY(&zgids, av);
  }

  // To do the summation, we'd need this function to be UpdateAdd or 
  // Update with an operator argument like Tpetra::ADD.
  // For now, we'll just do the update (which does a replace-like operation).
  zz.Update(zgids, NULL, (char*)(nIds ? &(user[0]) : NULL), NULL, inIds);

  // Retrieve the result for all local IDs.
  zz.Find(zgids, NULL, (char*)(nIds ? &(user[0]) : NULL), NULL, inIds, NULL);

  // The following step is needed only to test the results;
  // for general use, if user[i] > 1 in the summation, id[i] is shared.
  size_t cntShared = 0;
  std::unordered_set<id_t> alreadyCounted;
  for (size_t i = 0; i < nIds; i++) {
    if (user[i] > 1) {
      // Id is shared; have we already counted it locally?
      if (alreadyCounted.find(ids[i]) == alreadyCounted.end()) {
        alreadyCounted.insert(ids[i]);
        cntShared++;
      }
    }
  }

  if ((nIds * comm->getSize()) <= TOOMANY) zz.Print();

  return (cntShared == nShared);
}

//////////////////////////////////////////////////////////////////////////////
int main(int narg, char **arg)
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  // Test with contiguous IDs of default Tpetra GO type
  {
    size_t nIds = 5;        // Number of IDs per processor
    float fracShared = 0.2; // Fraction of IDs that should be shared
    size_t idBase = 0;      // Smallest possible ID
    int idStride = 1;       // Offset between IDs; 1 gives contiguous numbering

    typedef Tpetra::Map<>::global_ordinal_type gno_t;
    IDs<gno_t> myIds(nIds, fracShared, idBase, idStride, comm);

    myIds.TpetraDDTest();

    myIds.ZoltanDDTest();
  }

  // Test with non-contiguous IDs of default Tpetra GO starting at 20
  {
    size_t nIds = 5;        // Number of IDs per processor
    float fracShared = 0.2; // Fraction of IDs that should be shared
    size_t idBase = 20;      // Smallest possible ID
    int idStride = 3;       // Offset between IDs; 1 gives contiguous numbering

    typedef Tpetra::Map<>::global_ordinal_type gno_t;
    IDs<gno_t> myIds(nIds, fracShared, idBase, idStride, comm);

    myIds.TpetraDDTest();

    myIds.ZoltanDDTest();
  }

#ifdef HAVE_TPETRA_INT_INT
  // Test with contiguous integer IDs
  {
    size_t nIds = 5;        // Number of IDs per processor
    float fracShared = 0.2; // Fraction of IDs that should be shared
    size_t idBase = 0;      // Smallest possible ID
    int idStride = 1;       // Offset between IDs; 1 gives contiguous numbering

    IDs<int> myIds(nIds, fracShared, idBase, idStride, comm);

    myIds.TpetraDDTest();

    myIds.ZoltanDDTest();
  }

  // Test with non-contiguous integer IDs starting at 20
  {
    size_t nIds = 5;        // Number of IDs per processor
    float fracShared = 0.2; // Fraction of IDs that should be shared
    size_t idBase = 20;      // Smallest possible ID
    int idStride = 3;       // Offset between IDs; 1 gives contiguous numbering

    IDs<int> myIds(nIds, fracShared, idBase, idStride, comm);

    myIds.TpetraDDTest();

    myIds.ZoltanDDTest();
  }
#endif

#ifdef HAVE_TPETRA_INT_LONG_LONG
  // Test with non-contiguous long long IDs starting at 200
  {
    size_t nIds = 5;        // Number of IDs per processor
    float fracShared = 0.4; // Fraction of IDs that should be shared
    size_t idBase = 200;    // Smallest possible ID
    int idStride = 4;       // Offset between IDs; 1 gives contiguous numbering

    IDs<long long> myIds(nIds, fracShared, idBase, idStride, comm);

    myIds.TpetraDDTest();

    myIds.ZoltanDDTest();
  }
#endif
}
