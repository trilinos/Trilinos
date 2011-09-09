#ifndef stk_encr_FieldFunction_hpp
#define stk_encr_FieldFunction_hpp

#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <utility>

#include <stk_percept/Percept.hpp>

#include <Intrepid_Basis.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/fem/CoordinateSystems.hpp>

#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/internal/Searcher.hpp>
#include "Teuchos_RCP.hpp"

#include <stk_percept/PerceptMesh.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_search/diag/IdentProc.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>

#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

//using namespace sierra;
using namespace Intrepid;

namespace stk
{
  namespace percept
  {
    
    //class Helper;
    /** Evaluate the function at this input point (or points) returning value(s) in output_field_values
     *
     *   In the following, the arrays are dimensioned using the notation (from Intrepid's doc):
     *
     *   [C]         - num. integration domains (cells/elements)
     *   [F]         - num. Intrepid "fields" (number of bases within an element == num. nodes typically)
     *   [P]         - num. integration (or interpolation) points within the element
     *   [D]         - spatial dimension
     *   [D1], [D2]  - spatial dimension
     *
     *   Locally, we introduce this notation:
     *
     *   [DOF]       - number of degrees-of-freedom per node of the interpolated stk Field.  For example, a vector field in 3D has [DOF] = 3
     *
     */
    class FieldFunction : public Function 
    {
    public:
      static bool m_parallelEval;

      enum SearchType { SIMPLE_SEARCH, STK_SEARCH };
      
      //========================================================================================================================
      // high-level interface
      FieldFunction(const char *name, mesh::FieldBase *field, PerceptMesh& mesh, 
                    int domain_dimension=3,
                    int codomain_dimension=1,
                    SearchType searchType = SIMPLE_SEARCH,
                    unsigned integration_order = 0);

      mesh::FieldBase *getField();
      void interpolateFrom(Function& function);

      //virtual void value(MDArray& in, MDArray& out, double time_value_optional=0.0);

      //========================================================================================================================
      // low-level interface
      FieldFunction(const char *name, mesh::FieldBase *field, mesh::BulkData *bulk, 
                    Dimensions domain_dimensions = Dimensions(),
                    Dimensions codomain_dimensions = Dimensions(),
                    SearchType searchType = SIMPLE_SEARCH,
                    unsigned integration_order = 0);

      FieldFunction(const char *name, mesh::FieldBase *field, PerceptMesh& eMesh,
                    Dimensions domain_dimensions, // = Dimensions(),
                    Dimensions codomain_dimensions, // = Dimensions(),
                    SearchType searchType = SIMPLE_SEARCH,
                    unsigned integration_order = 0);

      virtual ~FieldFunction();

      virtual void operator()(MDArray& in, MDArray& out, double time_value_optional=0.0);
      virtual void localEvaluation(MDArray& in, MDArray& out, double time_value_optional=0.0);

      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Entity& element, const MDArray& parametric_coords, double time_value_optional=0.0);
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0);

      template<class BucketOrEntity>
      void helper(MDArray& input_phy_points, MDArray& output_field_values,
                  const BucketOrEntity& bucket_or_element, const MDArray& parametric_coordinates, double time_value_optional);

      mesh::BulkData *getBulkData();

      //void setBulkData(mesh::BulkData *bulk) { m_bulkData = bulk; }
      bool getFoundOnLocalOwnedPart() { return m_found_on_local_owned_part; }

    private:
      mesh::FieldBase *m_my_field;
      mesh::BulkData *m_bulkData;
      const mesh::Entity *m_cachedElement;
      Searcher* m_searcher;
      //const CellTopologyData * m_cached_topo;

      //typedef Intrepid::Basis<double, MDArray > IntrepidBasisType;
      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;

      unsigned m_cached_topo_key;
      BasisTypeRCP m_cached_basis;

      SearchType m_searchType;
      bool m_found_on_local_owned_part;

    };

    /** Evaluate the function on this element at the parametric coordinates and return in output_field_values.
     *
     *  Dimensions of parametric_coordinates are required to be ([P],[D])
     *  Dimensions of output_field_values are required to be ([P],[DOF]), (or in future, ([C],[P],[DOF]) )
     *
     */
#define EXTRA_PRINT_FF_HELPER 0
    template<class BucketOrEntity>
    void FieldFunction::helper(MDArray& input_phy_points, MDArray& output_field_values,
                               const BucketOrEntity& bucket_or_element, const MDArray& parametric_coordinates, double time_value_optional) 
    {
      VERIFY_OP(parametric_coordinates.rank(), ==, 2, "FieldFunction::operator() parametric_coordinates bad rank");
      VERIFY_OP(output_field_values.rank(), ==, 2, "FieldFunction::operator() output_field_values bad rank");

      int numInterpPoints = parametric_coordinates.dimension(0);

      VERIFY_OP(output_field_values.dimension(0), ==, numInterpPoints, "FieldFunction::operator() output_field_values bad dim(0)");

      const CellTopologyData * const cell_topo_data = PerceptMesh::get_cell_topology(bucket_or_element);

      unsigned stride = 0;
      //double * fdata_bucket = PerceptMesh::field_data( m_my_field , bucket, &stride);
      // intentionally ignoring return value to get around compiler warning
      PerceptMesh::field_data( m_my_field , bucket_or_element, &stride);

      unsigned nDOF = stride;
#ifndef NDEBUG
      int nOutDim = m_codomain_dimensions.back(); // FIXME for tensor
      // FIXME 
      VERIFY_OP((int)nDOF, == , nOutDim,
                "FieldFunction::operator(): invalid dimensions nDOF, m_codomain_dimensions[0]= ");
#endif

      int numCells = PerceptMesh::size1(bucket_or_element); // FIXME for multiple cells

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 1" << std::endl;

      shards::CellTopology topo(cell_topo_data);
      int numNodes = topo.getNodeCount();
      int cellDim  = topo.getDimension();
      if (0)
        {
          MDArray cellWorkset(numCells, numNodes, cellDim);
          if (0) cellWorkset(0,0,0) = 0.0;
        }
      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 2" << std::endl;

      // map cell topology to a basis
      PerceptMesh::BasisTypeRCP basis;
      if (m_cached_topo_key != cell_topo_data->key)
        {
          basis = PerceptMesh::getBasis(topo);
          m_cached_basis = basis;
          m_cached_topo_key = cell_topo_data->key;
        }
      else
        {
          basis = m_cached_basis;
        }

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 3" << std::endl;
      VERIFY_OP(basis.get(), != , 0, "FieldFunction::operator() basis is null");
      

      int numBases = basis->getCardinality();
      if (numBases != numNodes)
        {
          throw std::runtime_error(" (numBases != numNodes) ");
        }

      // [P] = 1
      //int numInterpPoints = 1;    // FIXME - this is now set based on parametric_coordinates dimensioning

      // ([F],[P]), or ([F],[P],[D]) for GRAD
      MDArray basis_values(numBases, numInterpPoints); 

      // ([C],[F],[P]), or ([C],[F],[P],[D]) for GRAD
      MDArray transformed_basis_values(numCells, numBases, numInterpPoints); 

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 4" << std::endl;

      // FIXME - it appears that Intrepid only supports the evaluation of scalar-valued fields, so we have
      //   to copy the field one DOF at a time into a local array, evaluate, then copy back
      // ([C],[F])
      MDArray field_data_values(numCells, numBases);
      MDArray field_data_values_dof(numCells, numBases, nDOF);

      //const mesh::PairIterRelation elem_nodes = bucket_or_element.relations( mesh::Node );

      // ([P],[D])  [P] points in [D] dimensions
      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 4a, par= \n" << parametric_coordinates << std::endl;

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) parametric_coordinates = \n " << parametric_coordinates << std::endl;
      {
        EXCEPTWATCH;
        basis->getValues(basis_values, parametric_coordinates, Intrepid::OPERATOR_VALUE);
      }

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) basis_values = \n " << basis_values << std::endl;

      // this function just spreads (copies) the values of the basis to all elements in the workset (numCells)
      FunctionSpaceTools::HGRADtransformVALUE<double>(transformed_basis_values, basis_values);
      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) transformed_basis_values =  " << transformed_basis_values << std::endl;

      // ([C],[P]) - place for results of evaluation
      MDArray loc_output_field_values(numCells, numInterpPoints);

      PerceptMesh::fillCellNodes(bucket_or_element, m_my_field, field_data_values_dof);

      // gather
      for (unsigned iDOF = 0; iDOF < nDOF; iDOF++)
        {
          for (int iCell = 0; iCell < numCells; iCell++)
            {
              for (int iNode = 0; iNode < numNodes; iNode++)
                {
                  field_data_values(iCell, iNode) = field_data_values_dof(iCell, iNode, iDOF);
                  if (EXTRA_PRINT_FF_HELPER) std::cout << "tmp iNode= " << iNode << "iDOF= " << iDOF << " fd= " << field_data_values(iCell, iNode) << std::endl;
                  if (0 && iDOF==0 && iNode < 4)
                    std::cout << "tmp iNode= " << iNode << " coord= " 
                              << field_data_values_dof(iCell, iNode, 0) << " "
                              << field_data_values_dof(iCell, iNode, 1) << " "
                      //<< field_data_values_dof(iCell, iNode, 2) << " "
                              << std::endl;
                }
            }

          /// NOTE: this is needed since FunctionSpaceTools::evaluate method assumes the output array is initialized to 0
          if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) evaluate ... " << std::endl;
          loc_output_field_values.initialize(0.0);
          FunctionSpaceTools::evaluate<double>(loc_output_field_values, field_data_values, transformed_basis_values);
          if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) evaluate done " << std::endl;

          for (int iCell = 0; iCell < numCells; iCell++)
            {
              for (int iPoint = 0; iPoint < numInterpPoints; iPoint++)
                {
                  //output_field_values(iCell, iPoint, iDOF) = loc_output_field_values(iCell, iPoint);
                  //output_field_values(0, iDOF) = loc_output_field_values(iCell, iPoint);
                  output_field_values(iPoint, iDOF) = loc_output_field_values(iCell, iPoint);
                  if (EXTRA_PRINT_FF_HELPER) std::cout << "tmp iDOF= " << iDOF << " ofd= " << output_field_values(iPoint, iDOF) << std::endl;
                }
            }
          if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) scatter done " << std::endl;
        }
    }


  }//namespace percept
}//namespace stk
#endif
