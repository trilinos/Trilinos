/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <percept/uq/Percept_API_KLSolver.hpp>
#include <percept/PerceptUtils.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <Shards_CellTopologyData.h>

namespace percept
{

Percept_API_KLSolver::Percept_API_KLSolver(
   const stk::mesh::BulkData & mesh,
   const stk::mesh::Field<double> & phi,
   std::vector<double> & lambda)
  : 
  RFGen::API_KLSolver(),
  m_mesh(mesh),
  m_meta(m_mesh.mesh_meta_data()),
  m_phi(phi),
  m_lambda(lambda)
{}

unsigned 
Percept_API_KLSolver::getSpatialDim() const
{
  return m_meta.spatial_dimension();
}

void 
Percept_API_KLSolver::computeLocalIntgDataSizes(
  int &localNumElem,
  int &localMaxIntgPts)
{
  std::vector<size_t> count;
  stk::mesh::Selector selector(m_meta.locally_owned_part());
   
  stk::mesh::count_entities( selector, m_mesh, count );

  localNumElem = count[stk::topology::ELEMENT_RANK];

  // one point integration to start
  localMaxIntgPts = 1;
}

void 
Percept_API_KLSolver::computeLocalIntgData(
  shards::Array<double,shards::NaturalOrder,RFGen::Cell,RFGen::Point,RFGen::Dim> &localIntgPtCoords,
  shards::Array<double,shards::NaturalOrder,RFGen::Cell,RFGen::Point> &localVolumeWeights)
{
  typedef stk::mesh::Field<double> VectorFieldType;
  const VectorFieldType * coordsField = static_cast<const VectorFieldType*>(m_meta.coordinate_field());

  double centroid[3];
  const int spatialDim = m_meta.spatial_dimension();

  int cell_index = 0;

  stk::mesh::Selector selector(m_meta.locally_owned_part());
  const stk::mesh::BucketVector & buckets = m_mesh.get_buckets(stk::topology::ELEMENT_RANK, selector);
  for (stk::mesh::BucketVector::const_iterator k = buckets.begin(); k != buckets.end(); ++k ) {
    stk::mesh::Bucket & bucket = **k;

    const int num_elements_in_bucket = bucket.size();
    const CellTopologyData * cell_topo_data = stk::mesh::get_cell_topology(bucket.topology()).getCellTopologyData();
    
    for (int iElement = 0; iElement < num_elements_in_bucket; iElement++) {
      stk::mesh::Entity element = bucket[iElement];
      
      computeCentroid(element, centroid, *coordsField);
      
      for (int d=0; d<spatialDim; d++) {
        localIntgPtCoords(cell_index,0,d) = centroid[d];
      }
      localVolumeWeights(cell_index,0) = volume(element, coordsField, cell_topo_data);

      cell_index++;
    }
  }
}

MPI_Comm 
Percept_API_KLSolver::getParallelComm() const
{
  return m_mesh.parallel();
}

void 
Percept_API_KLSolver::setKLSolution(
  const int &numTerms,
  const shards::Array<double,shards::NaturalOrder,RFGen::Eigen> &eigenValues,
  const shards::Array<double,shards::NaturalOrder,RFGen::Eigen,RFGen::Cell> &eigenVectors)
{
  for (int ev=0; ev<numTerms; ev++) {
    m_lambda[ev] = eigenValues(ev);
  }
  
  int cell_index = 0;

  stk::mesh::Selector selector(m_meta.locally_owned_part());
  const stk::mesh::BucketVector & buckets = m_mesh.get_buckets(stk::topology::ELEMENT_RANK, selector);
  for (stk::mesh::BucketVector::const_iterator k = buckets.begin(); k != buckets.end(); ++k ) {
    stk::mesh::Bucket & bucket = **k;

    const int num_elements_in_bucket = bucket.size();
    
    double * phi = stk::mesh::field_data(m_phi, bucket);

    for (int iElement = 0; iElement < num_elements_in_bucket; iElement++) {
      for (int ieigen=0; ieigen<numTerms; ieigen++) {
        phi[numTerms*iElement + ieigen] = eigenVectors(ieigen, cell_index);
      }
      cell_index++;
    }
  }     
}
  
}
