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
using Teuchos::Array;

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
  int fail = 0;

  // Get some coordinates

  typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> mv_t;

  outputFlag_t flags;
  flags.set(OBJECT_COORDINATES);

  RCP<UserInputForTests> uinput;
  std::string fname = testDataFilePath+std::string("/simple.mtx");

  try{
    uinput = rcp(new UserInputForTests(fname, comm, flags));
  }
  catch(std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "input constructor", 1);

  RCP<mv_t> coords;

  try{
    coords = uinput->getCoordinates();
  }
  catch(std::exception &e){
    fail=1;
  }

  TEST_FAIL_AND_EXIT(*comm, !fail, "getting coordinates", 1);

  int numLocalIds = coords->getLocalLength();
  int numGlobalIds = coords->getGlobalLength();
  int coordDim = coords->getNumVectors();
  ArrayView<const gno_t> idList = coords->getMap()->getNodeElementList();

  // Create global Ids, x-, y- and z-coordinates, and also arrays of weights.

  Array<gno_t> myIds(numLocalIds);
  gno_t base = rank * numLocalIds;
  
  int wdim = 2;
  Array<scalar_t> weights(numLocalIds*wdim);

  scalar_t *x_values= coords->getDataNonConst(0).getRawPtr();
  scalar_t *y_values= x_values;  // fake 3 dimensions if needed
  scalar_t *z_values= x_values;

  if (coordDim > 1){
    y_values= coords->getDataNonConst(1).getRawPtr();
    if (coordDim > 2)
      z_values= coords->getDataNonConst(2).getRawPtr();
  }

  Array<scalar_t> xyz_values(3*numLocalIds);

  for (lno_t i=0; i < numLocalIds; i++)   // global Ids
    myIds[i] = base+i;

  scalar_t *x = xyz_values.getRawPtr();   // a stride-3 coordinate array
  scalar_t *y = x+1;
  scalar_t *z = y+1;

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
       numLocalIds, myIds.getRawPtr(), x_values, y_values, z_values));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 0", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
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
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 1", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 1", fail);
  
    // Try using the default: no strides supplied means strides are one.

    std::vector<int> emptyStrides;
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, emptyStrides, 
       weightValues, emptyStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 2", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 2", fail);
  }

  {
    ////////////////////////////////////////////////////////////////
    // 2-dimensional coordinates with stride three and two weights
  
    int ncoords = 2;
    int nweights = 2;

    std::vector<const scalar_t *> values, weightValues;
    std::vector<int> valueStrides, weightStrides;
  
    values.push_back(xyz_values.getRawPtr());
    values.push_back(xyz_values.getRawPtr() + 1);
    valueStrides.push_back(3);
    valueStrides.push_back(3);
  
    weightValues.push_back(weights.getRawPtr());
    weightValues.push_back(weights.getRawPtr() + numLocalIds);
    weightStrides.push_back(1);
    weightStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 3", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 3", fail);
  
    // Try using default weight strides

    std::vector<int> emptyStrides;
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, emptyStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
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
  
    weightValues.push_back(weights.getRawPtr());
    weightValues.push_back(weights.getRawPtr() + numLocalIds);
    weightStrides.push_back(1);
    weightStrides.push_back(1);
  
    try{
     ia = rcp(new Zoltan2::BasicCoordinateInput<userTypes_t>(
       numLocalIds, myIds.getRawPtr(), values, valueStrides, 
       weightValues, weightStrides));
    }
    catch (std::exception &e){
      fail = 1;
    }
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "constructor 4", fail);
  
    fail = checkBasicCoordinate(ia.getRawPtr(), numLocalIds, numGlobalIds,
      myIds.getRawPtr(), xyz_values.getRawPtr(), 
      weights.getRawPtr(), ncoords, nweights);
  
    TEST_FAIL_AND_RETURN_VALUE(*comm, fail==0, "check adapter 4", fail);
  }

  if (rank == 0)
    std::cout << "PASS" << std::endl;

  return fail;
}

