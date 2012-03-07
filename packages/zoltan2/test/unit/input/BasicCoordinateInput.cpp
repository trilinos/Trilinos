// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Test for Zoltan2::BasicCoordinateInput 

#include <Zoltan2_BasicCoordinateInput.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_SerialDenseVector.hpp>

using Teuchos::RCP;
using Teuchos::Comm;
using Teuchos::DefaultComm;

typedef Teuchos::SerialDenseVector<lno_t, scalar_t> tvec_t;

typedef Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> userTypes_t;

int checkBasicCoordinate(
  Zoltan2::BasicCoordinateInput<userTypes_t> *ia, 
  int len, int glen, gno_t *ids,
  scalar_t *xyz,
  scalar_t *weights,
  int nCoords, int nWeights)
{
  int fail = 0;

  if (ia->getCoordinateDimension() != nCoords)
    fail = 100;

  if (!fail && ia->getNumberOfWeights() != nWeights)
    fail = 101;

  if (!fail && ia->getLocalNumberOfCoordinates() != len)
    fail = 102;

  if (!fail && ia->getGlobalNumberOfCoordinates() != glen)
    fail = 103;

  for (int x=0; !fail && x < nCoords; x++){
    const gno_t *idList;
    const scalar_t *vals;
    int stride;

    size_t nvals = ia->getCoordinates(x, idList, vals, stride);

    if (nvals != len*stride)
      fail = 104;

    scalar_t *coordVal = xyz + x;
    for (int i=0; !fail && i < len; i++, coordVal += 3){

      if (idList[i] != ids[i])
        fail = 105;

      if (!fail && vals[stride*i] != *coordVal)
        fail = 106;
    }
  }

  for (int w=0; !fail && w < nWeights; w++){
    const scalar_t *wgts;
    int stride;

    size_t nvals = ia->getCoordinateWeights(w, wgts, stride);

    if (nvals != len)
      fail = 108;

    scalar_t *weightVal = weights + len*w;
    for (int i=0; !fail && i < len; i++, weightVal++){
      if (wgts[stride*i] != *weightVal)
        fail = 110;
    }
  }

  return fail;
}
  

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0;

  // Create global Ids, x-, y- and z-coordinates, and also arrays of weights.

  lno_t numLocalIds = 10;
  gno_t *myIds = new gno_t [numLocalIds];
  gno_t base = rank * numLocalIds;
  
  int wdim = 2;
  scalar_t *weights = new scalar_t [numLocalIds*wdim];
  const scalar_t **weightPtrs = new const scalar_t * [wdim];

  scalar_t *x_values= new scalar_t [numLocalIds];
  scalar_t *y_values= new scalar_t [numLocalIds];
  scalar_t *z_values= new scalar_t [numLocalIds];
  scalar_t *xyz_values = new scalar_t [3*numLocalIds];

  for (lno_t i=0; i < numLocalIds; i++)   // global Ids
    myIds[i] = base+i;

  tvec_t weightVec(Teuchos::View, weights, numLocalIds*wdim); // random weights 
  weightVec.random(); 

  for (int w=0; w < wdim; w++){
    weightPtrs[w] = weights + numLocalIds*w;
  }

  tvec_t xVec(Teuchos::View, x_values, numLocalIds); // random x coordinates
  xVec.random();

  tvec_t yVec(Teuchos::View, y_values, numLocalIds); // random y coordinates
  yVec.random();

  tvec_t zVec(Teuchos::View, z_values, numLocalIds); // random z coordinates
  zVec.random();

  scalar_t *x = xyz_values;                // a stride-3 coordinate array
  scalar_t *y = xyz_values + 1;
  scalar_t *z = xyz_values + 2;

  for (int i=0, ii=0; i < numLocalIds; i++, ii += 3){
    x[ii] = x_values[i];
    y[ii] = y_values[i];
    z[ii] = z_values[i];
  }

  RCP<Zoltan2::BasicCoordinateInput<userTypes_t> > ia;

  {
    ////////////////////////////////////////////////////////////////
    // 3-dimensional coordinates with stride one and no weights,
    //   using simpler constructor
  
    int ncoords = 3;
    int nweights = 0;
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds, x_values, y_values, z_values));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 0", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, 
      numLocalIds*nprocs, 
      myIds, xyz_values, weights, ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 0", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 3-dimensional coordinates with stride one and no weights
  
    int ncoords = 3;
    int nweights = 0;

    std::vector<const scalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(x_values);
    values.push_back(y_values);
    values.push_back(z_values);
    valueStrides.push_back(1);
    valueStrides.push_back(1);
    valueStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds, values, valueStrides, weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 1", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, 
      numLocalIds*nprocs, myIds, xyz_values, weights, ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 1", fail);
  
    // Try using the default: no strides supplied means strides are one.

    std::vector<int> emptyStrides;
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds, values, emptyStrides, weightValues, emptyStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 2", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, 
      numLocalIds*nprocs, myIds, xyz_values, weights, ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 2", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 2-dimensional coordinates with stride three and two weights
  
    int ncoords = 2;
    int nweights = 2;

    std::vector<const scalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(xyz_values);
    values.push_back(xyz_values + 1);
    valueStrides.push_back(3);
    valueStrides.push_back(3);
  
    weightValues.push_back(weights);
    weightValues.push_back(weights + numLocalIds);
    weightStrides.push_back(1);
    weightStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds, values, valueStrides, weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 3", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, 
      numLocalIds*nprocs, myIds, xyz_values, weights, ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 3", fail);
  
    // Try using default weight strides

    std::vector<int> emptyStrides;
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds, values, valueStrides, weightValues, emptyStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, 
      numLocalIds*nprocs, myIds, xyz_values, weights, ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 4", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 1-dimensional coordinates with stride one and two weights

    int ncoords = 1;
    int nweights = 2;

    std::vector<const scalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(x_values);
    valueStrides.push_back(1);
  
    weightValues.push_back(weights);
    weightValues.push_back(weights + numLocalIds);
    weightStrides.push_back(1);
    weightStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds, values, valueStrides, weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, 
      numLocalIds*nprocs, myIds, xyz_values, weights, ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 4", fail);
  }

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return fail;
}

