/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef RFGen_API_KLSolver_h
#define RFGen_API_KLSolver_h

#include "mpi.h"

#include "RFGen_Shards.h"

namespace RFGen
{

class API_KLSolver
{
 public:
  explicit API_KLSolver() {}

  virtual ~API_KLSolver() {}

  virtual unsigned getSpatialDim() const  = 0;

  virtual void computeLocalIntgDataSizes(
    int &localNumElem,
    int &localMaxIntgPts) = 0;

  virtual void computeLocalIntgData(
    shards::Array<double,shards::NaturalOrder,Cell,Point,Dim> &localIntgPtCoords,
    shards::Array<double,shards::NaturalOrder,Cell,Point> &localVolumeWeights) = 0;

  virtual MPI_Comm getParallelComm() const = 0;

  // app must copy eigenvectors to its native data structures
  virtual void setKLSolution(
    const int &numTerms,
    const shards::Array<double,shards::NaturalOrder,Eigen> &eigenValues,
    const shards::Array<double,shards::NaturalOrder,Eigen,Cell> &eigenValuesCell) = 0;
};

}

#endif // RFGen_API_KLSolver_h
