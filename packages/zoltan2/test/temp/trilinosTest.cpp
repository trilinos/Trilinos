// @HEADER
// ***********************************************************************
//                Copyright message goes here. 
// ***********************************************************************
// @HEADER

#include <Tpetra_MultiVector.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_TimeMonitor.hpp>

using std::cout;
using std::cerr;
using std::endl;
using std::bad_alloc;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Comm;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::CommandLineProcessor;
using Teuchos::TimeMonitor;

typedef int lno_t;
typedef long gno_t;
typedef double scalar_t;

typedef RCP<Teuchos::Time> ttime_t;
typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t> tMVector_t;
typedef Tpetra::Map<lno_t, gno_t> tMap_t;

/*! \brief Create a mesh of approximately the desired size.
 *
 *  We want 3 dimensions close to equal in length.
 */
ArrayRCP<const ArrayView<const scalar_t> > getMeshCoordinates(
    const RCP<const Teuchos::Comm<int> > & comm,
    gno_t numGlobalCoords)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  double k = log(numGlobalCoords) / 3;
  double xdimf = exp(k) + 0.5;
  gno_t xdim = static_cast<int>(floor(xdimf));
  gno_t ydim = xdim;
  gno_t zdim = numGlobalCoords / (xdim*ydim);
  gno_t num=xdim*ydim*zdim;
  gno_t diff = numGlobalCoords - num;
  gno_t newdiff = 0;

  while (diff > 0){
    if (zdim > xdim && zdim > ydim){
      zdim++;
      newdiff = diff - (xdim*ydim);
      if (newdiff < 0)
        if (diff < -newdiff)
          zdim--;
    }
    else if (ydim > xdim && ydim > zdim){
      ydim++;
      newdiff = diff - (xdim*zdim);
      if (newdiff < 0)
        if (diff < -newdiff)
          ydim--;
    }
    else{
      xdim++;
      newdiff = diff - (ydim*zdim);
      if (newdiff < 0)
        if (diff < -newdiff)
          xdim--;
    }

    diff = newdiff;
  }

  num=xdim*ydim*zdim;
  diff = numGlobalCoords - num;
  if (diff < 0)
    diff /= -numGlobalCoords;
  else
    diff /= numGlobalCoords;

  if (rank == 0){
    if (diff > .01)
      cout << "Warning: Difference " << diff*100 << " percent" << endl;
    cout << "Mesh size: " << xdim << "x" << ydim << "x" <<
      zdim << ", " << num << " vertices." << endl;
  }

  // Divide coordinates.

  gno_t numLocalCoords = num / nprocs;
  gno_t leftOver = num % nprocs;
  gno_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = gno_t(rank) * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) + 
           ((gno_t(rank) - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  gno_t gid1 = gid0 + numLocalCoords;

  gno_t *ids = new gno_t [numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayRCP<gno_t> idArray(ids, 0, numLocalCoords, true);

  for (gno_t i=gid0; i < gid1; i++)
    *ids++ = i;   

  RCP<const tMap_t> idMap = rcp(
    new tMap_t(num, idArray.view(0, numLocalCoords), 0, comm));

  // Create a Tpetra::MultiVector of coordinates.

  scalar_t *x = new scalar_t [numLocalCoords*3]; 
  if (!x)
    throw bad_alloc();
  ArrayRCP<scalar_t> coordArray(x, 0, numLocalCoords*3, true);

  scalar_t *y = x + numLocalCoords;
  scalar_t *z = y + numLocalCoords;

  gno_t xStart = 0;
  gno_t yStart = 0;
  gno_t xyPlane = xdim*ydim;
  gno_t zStart = gid0 / xyPlane;
  gno_t rem = gid0 % xyPlane;
  if (rem > 0){
    yStart = rem / xdim;
    xStart = rem % xdim;
  }

  lno_t next = 0;
  for (scalar_t zval=zStart; next < numLocalCoords && zval < zdim; zval++){
    for (scalar_t yval=yStart; next < numLocalCoords && yval < ydim; yval++){
      for (scalar_t xval=xStart; next < numLocalCoords && xval < xdim; xval++){
        x[next] = xval;
        y[next] = yval;
        z[next] = zval;
        next++;
      }
      xStart = 0;
    }
    yStart = 0;
  }

  ArrayView<const scalar_t> xArray(x, numLocalCoords);
  ArrayView<const scalar_t> yArray(y, numLocalCoords);
  ArrayView<const scalar_t> zArray(z, numLocalCoords);
  ArrayRCP<ArrayView<const scalar_t> > coordinates =
    arcp(new ArrayView<const scalar_t> [3], 0, 3);
  coordinates[0] = xArray;
  coordinates[1] = yArray;
  coordinates[2] = zArray;

  ArrayRCP<const ArrayView<const scalar_t> > constCoords =
   coordinates.getConst();

  return constCoords;
}


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();

  ttime_t t1 = TimeMonitor::getNewTimer("contig map");
  ttime_t t2 = TimeMonitor::getNewTimer("contig mvector");
  ttime_t t3 = TimeMonitor::getNewTimer("non-contig map");
  ttime_t t4 = TimeMonitor::getNewTimer("non-contig mvector");

  if (rank==0)
    cout << "Number of processes: " << nprocs << endl;

  ///////////////////////////////////
  // Get global number of coords.
  ///////////////////////////////////
  double val = 1000;
  CommandLineProcessor commandLine(false, true);
  commandLine.setOption("size", &val, 
    "Approximate number of global coordinates.");

  commandLine.parse(argc, argv);

  gno_t numGlobalCoords = static_cast<gno_t>(val);

  // Create 3-d coords.
  ArrayRCP<const ArrayView<const scalar_t> > coordinates =
    getMeshCoordinates(comm, numGlobalCoords);

  ///////////////////////////////////
  // Divide coordinates and create a map (contiguous)
  ///////////////////////////////////

  gno_t numLocalCoords = numGlobalCoords / nprocs;
  gno_t leftOver = numGlobalCoords % nprocs;
  gno_t gid0 = 0;

  if (rank <= leftOver)
    gid0 = gno_t(rank) * (numLocalCoords+1);
  else
    gid0 = (leftOver * (numLocalCoords+1)) +
           ((gno_t(rank) - leftOver) * numLocalCoords);

  if (rank < leftOver)
    numLocalCoords++;

  gno_t gid1 = gid0 + numLocalCoords;

  gno_t *ids = new gno_t [numLocalCoords];
  if (!ids)
    throw bad_alloc();
  ArrayRCP<gno_t> idArray(ids, 0, numLocalCoords, true);

  for (gno_t i=gid0; i < gid1; i++)
    *ids++ = i;

  t1->start();

  RCP<const tMap_t> idMap = rcp(
    new tMap_t(numGlobalCoords, idArray.view(0, numLocalCoords), 0, comm));

  t1->stop();

  ///////////////////////////////////
  // Create multivector
  ///////////////////////////////////

  t2->start();

  RCP<tMVector_t> meshCoords = rcp(new tMVector_t(
    idMap, coordinates.view(0,3), 3));

  t2->stop();

  ///////////////////////////////////
  // New map, non-contiguous
  ///////////////////////////////////

  gno_t idx = 0;
  ids = idArray.getRawPtr();
  for (gno_t i=rank; i < numGlobalCoords; i+=nprocs, idx++)
    ids[idx] = i;

  t3->start();

  idMap = rcp(new tMap_t(
    numGlobalCoords, idArray.view(0, numLocalCoords), 0, comm));

  t3->stop();
  ///////////////////////////////////
  // New multivector
  ///////////////////////////////////

  t4->start();

  meshCoords = rcp(new tMVector_t(idMap, coordinates.view(0,3), 3));

  t4->stop();

  TimeMonitor::summarize(std::cout);
  
  return 0;
}

