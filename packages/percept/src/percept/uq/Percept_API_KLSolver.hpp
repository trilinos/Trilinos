/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef percept_Percept_API_KLSolver_h
#define percept_Percept_API_KLSolver_h

#include <percept/rfgen/RFGen_API_KLSolver.h>

#include <stk_mesh/base/Field.hpp>

namespace stk { namespace mesh { 
    class BulkData; 
}}

namespace percept
{

class Percept_API_KLSolver : public RFGen::API_KLSolver
{
 public: 
  explicit 
  Percept_API_KLSolver(const stk::mesh::BulkData & mesh, 
                       const stk::mesh::Field<double> & phi,
                       std::vector<double> & lambda);

  ~Percept_API_KLSolver() {}

  unsigned getSpatialDim() const override;

  void computeLocalIntgDataSizes(
    int &localNumElem,
    int &localMaxIntgPts) override;

  void computeLocalIntgData(
    shards::Array<double,shards::NaturalOrder,RFGen::Cell,RFGen::Point,RFGen::Dim> &localIntgPtCoords,
    shards::Array<double,shards::NaturalOrder,RFGen::Cell,RFGen::Point> &localVolumeWeights) override;

  MPI_Comm getParallelComm() const override;

  void setKLSolution(
    const int &numTerms,
    const shards::Array<double,shards::NaturalOrder,RFGen::Eigen> &eigenValues,
    const shards::Array<double,shards::NaturalOrder,RFGen::Eigen,RFGen::Cell> &eigenVectors) override;

 private:
  const stk::mesh::BulkData & m_mesh;
  const stk::mesh::MetaData & m_meta;
  const stk::mesh::Field<double> & m_phi;
  std::vector<double> & m_lambda;
};

}

#endif // percept_Percept_API_KLSolver_h
