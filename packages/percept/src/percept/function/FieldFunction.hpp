// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_FieldFunction_hpp
#define stk_encr_FieldFunction_hpp

#include <string>
#include <cmath>
#include <math.h>
#include <vector>
#include <utility>

#include <percept/Percept.hpp>

#include <Intrepid_Basis.hpp>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/CoordinateSystems.hpp>

#include <percept/function/Function.hpp>
#include <percept/function/internal/Searcher.hpp>
#include "Teuchos_RCP.hpp"

#include <percept/PerceptMesh.hpp>

#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>

#include <Shards_CellTopology.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

#include <Intrepid_HGRAD_HEX_C1_FEM.hpp>

#include <Intrepid_CellTools.hpp>
#include <Intrepid_FunctionSpaceTools.hpp>

#include <percept/element/intrepid/BasisTable.hpp>

//using namespace sierra;
using namespace Intrepid;

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
      FieldFunction(const char *name, stk::mesh::FieldBase *field, PerceptMesh& mesh,
                    int domain_dimension=3,
                    int codomain_dimension=1,
                    SearchType searchType = SIMPLE_SEARCH,
                    unsigned integration_order = 0);

      stk::mesh::FieldBase *get_field();
      void interpolateFrom(Function& function);

      //virtual void value(MDArray& in, MDArray& out, double time_value_optional=0.0);

      //========================================================================================================================
      // low-level interface
      FieldFunction(const char *name, stk::mesh::FieldBase *field, stk::mesh::BulkData *bulk,
                    Dimensions domain_dimensions = Dimensions(),
                    Dimensions codomain_dimensions = Dimensions(),
                    SearchType searchType = SIMPLE_SEARCH,
                    unsigned integration_order = 0);

      FieldFunction(const char *name, stk::mesh::FieldBase *field, PerceptMesh& eMesh,
                    Dimensions domain_dimensions, // = Dimensions(),
                    Dimensions codomain_dimensions, // = Dimensions(),
                    SearchType searchType = SIMPLE_SEARCH,
                    unsigned integration_order = 0);

      virtual ~FieldFunction();

      virtual Teuchos::RCP<Function > derivative(MDArrayString& deriv_spec)
      {
        m_deriv_spec = deriv_spec;
        Dimensions domain_dimensions = getDomainDimensions();
        Dimensions codomain_dimensions = getCodomainDimensions();
        //int meta_dimension = stk::mesh::MetaData::get_meta_data(stk::mesh::MetaData::get(*m_bulkData)).get_spatial_dimension();
        int num_grad = deriv_spec.dimension(0);
        //domain_dimensions.push_back(num_grad);
        codomain_dimensions.back() = num_grad;

        FieldFunction *deriv = new FieldFunction(this->getName().c_str(),
                                                 this->m_my_field,
                                                 this->get_bulk_data(),
                                                 domain_dimensions, codomain_dimensions,
                                                 this->m_searchType,
                                                 this->getIntegrationOrder());
        deriv->m_deriv_spec = deriv_spec;
        deriv->m_get_derivative = true;
        Teuchos::RCP<Function > deriv_rcp = Teuchos::rcp(deriv);
        return deriv_rcp;
      }

      virtual Teuchos::RCP<Function > gradient(int spatialDim=3)
      {
        int meta_dimension = stk::mesh::MetaData::get(*m_bulkData).spatial_dimension();
        VERIFY_OP_ON(meta_dimension, ==, spatialDim, "gradient: mismatch in spatial dimensions");
        std::string xyz[] = {"x", "y", "z"};
        MDArrayString mda(meta_dimension);
        for (int i = 0; i < meta_dimension; i++)
          {
            mda(i) = xyz[i];
          }
        return derivative(mda);
      }

      virtual void operator()(MDArray& in, MDArray& out, double time_value_optional=0.0);
      virtual void localEvaluation(MDArray& in, MDArray& out, double time_value_optional=0.0);

      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0);
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0);

      void setup_searcher(int D_);

      template<class BucketOrEntity>
      void helper(const stk::mesh::BulkData& bulk, MDArray& input_phy_points, MDArray& output_field_values,
                  const BucketOrEntity& bucket_or_element, const MDArray& parametric_coordinates, double time_value_optional);

      stk::mesh::BulkData *get_bulk_data();

      //void setBulkData(stk::mesh::BulkData *bulk) { m_bulkData = bulk; }
      bool getFoundOnLocalOwnedPart() { return m_found_on_local_owned_part; }

      const stk::mesh::Bucket& mybucket(const stk::mesh::BulkData& bulkData, const stk::mesh::Bucket& bucket_or_element) const { return bucket_or_element; }
      const stk::mesh::Bucket& mybucket(const stk::mesh::BulkData& bulkData, const stk::mesh::Entity& bucket_or_element) const { return bulkData.bucket(bucket_or_element); }
      stk::mesh::Bucket& mybucket( stk::mesh::BulkData& bulkData,  stk::mesh::Bucket& bucket_or_element)  { return bucket_or_element; }
      stk::mesh::Bucket& mybucket( stk::mesh::BulkData& bulkData,  stk::mesh::Entity& bucket_or_element)  { return bulkData.bucket(bucket_or_element); }


    private:
      stk::mesh::FieldBase *m_my_field;
      stk::mesh::BulkData *m_bulkData;
      stk::mesh::Entity m_cachedElement;
      Searcher* m_searcher;
      //const CellTopologyData * m_cached_topo;

      //typedef Intrepid::Basis<double, MDArray > IntrepidBasisType;
      typedef Intrepid::Basis<double, MDArray > BasisType;
      typedef Teuchos::RCP<BasisType>           BasisTypeRCP;

      unsigned m_cached_topo_key;
      BasisTypeRCP m_cached_basis;

      SearchType m_searchType;
      bool m_found_on_local_owned_part;

      bool m_get_derivative;
      MDArrayString m_deriv_spec;

      stk::mesh::FieldBase *m_coordinatesField;

    };

    /** Evaluate the function on this element at the parametric coordinates and return in output_field_values.
     *
     *  Dimensions of parametric_coordinates are required to be ([P],[D])
     *  Dimensions of output_field_values are required to be ([P],[DOF]), (or in future, ([C],[P],[DOF]) )
     *
     */
#define EXTRA_PRINT_FF_HELPER 0
    template<class BucketOrEntity>
    void FieldFunction::helper(const stk::mesh::BulkData& bulk, MDArray& input_phy_points, MDArray& output_field_values,
                               const BucketOrEntity& bucket_or_element, const MDArray& parametric_coordinates, double time_value_optional)
    {
      VERIFY_OP_ON(parametric_coordinates.rank(), ==, 2, "FieldFunction::operator() parametric_coordinates bad rank");
      VERIFY_OP_ON(output_field_values.rank(), <=, 3, "FieldFunction::operator() output_field_values bad rank");

      int numInterpPoints = parametric_coordinates.dimension(0);
      int spatialDim = stk::mesh::MetaData::get(*m_bulkData).spatial_dimension();
      int num_grad = getCodomainDimensions().back();

      int numCells = PerceptMesh::size1(bucket_or_element); 

      if (output_field_values.rank() == 2)
        {
          VERIFY_OP_ON(output_field_values.dimension(0), ==, numInterpPoints, "FieldFunction::operator() output_field_values bad dim(0)");
        }
      else if (output_field_values.rank() == 3)
        {
          VERIFY_OP_ON(output_field_values.dimension(1), ==, numInterpPoints, "FieldFunction::operator() output_field_values bad dim(0)");
          VERIFY_OP_ON(output_field_values.dimension(0), ==, numCells, "FieldFunction::operator() output_field_values bad dim(0)");
        }

      const CellTopologyData * const cell_topo_data = stk::mesh::get_cell_topology(mybucket(bulk, bucket_or_element).topology()).getCellTopologyData();

      unsigned stride = 0;
      //  NKC, huh?, no-op..
      //stk::mesh::field_data( *m_my_field , bucket_or_element);

      if (1)
        {
          stk::mesh::EntityRank er = mybucket(bulk, bucket_or_element).entity_rank();
          const stk::mesh::FieldBase::Restriction & r = stk::mesh::find_restriction(*m_my_field, er, stk::mesh::MetaData::get(*m_my_field).universal_part());
          static const stk::mesh::FieldBase::Restriction empty ;

          if (r == empty)
            {
              unsigned nfr = m_my_field->restrictions().size();
              for (unsigned ifr = 0; ifr < nfr; ifr++)
                {
                  const stk::mesh::FieldRestriction& fr = m_my_field->restrictions()[ifr];
                  unsigned field_dimension = fr.num_scalars_per_entity() ;
                  if (field_dimension > 0)
                    {
                      stride = field_dimension;
                    }
                }
            }
          else
            {
              stride = r.num_scalars_per_entity() ;
            }
        }

      unsigned nDOF = stride;
#ifndef NDEBUG
      int nOutDim = m_codomain_dimensions.back(); // FIXME for tensor

      VERIFY_OP((int)(nDOF * (m_get_derivative? spatialDim : 1)), == , nOutDim,
                "FieldFunction::operator(): invalid dimensions nDOF, m_codomain_dimensions[0]= ");
#endif

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 1" << std::endl;

      shards::CellTopology topo(cell_topo_data);
      int numNodes = topo.getNodeCount();
      //int cellDim  = topo.getDimension();
      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 2" << std::endl;

      // map cell topology to a basis
      BasisTable::BasisTypeRCP basis;
      if (m_cached_topo_key != cell_topo_data->key)
        {
          basis = BasisTable::getBasis(topo);
          m_cached_basis = basis;
          m_cached_topo_key = cell_topo_data->key;
        }
      else
        {
          basis = m_cached_basis;
        }

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 3" << std::endl;
      if (basis.get() == 0)
        {
          std::cout << "FieldFunction::operator() basis is null, topo= " << topo.getName() << std::endl;
        }
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

      MDArray Jac;
      MDArray JacInverse;
      if (m_get_derivative)
        {
          basis_values.resize(numBases, numInterpPoints, spatialDim);
          transformed_basis_values.resize(numCells, numBases, numInterpPoints, spatialDim);
          MDArray cellWorkset(numCells, numNodes, spatialDim);
          PerceptMesh::fillCellNodes(*get_bulk_data(), bucket_or_element,  m_coordinatesField, cellWorkset, spatialDim);
          Jac.resize(numCells, numInterpPoints, spatialDim, spatialDim);
          JacInverse.resize(numCells, numInterpPoints, spatialDim, spatialDim);
          CellTools<double>::setJacobian(Jac, parametric_coordinates, cellWorkset, topo);
          CellTools<double>::setJacobianInv(JacInverse, Jac);
        }

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) 4" << std::endl;

      // FIXME - it appears that Intrepid only supports the evaluation of scalar-valued fields, so we have
      //   to copy the field one DOF at a time into a local array, evaluate, then copy back
      // ([C],[F])
      MDArray field_data_values(numCells, numBases);
      MDArray field_data_values_dof(numCells, numBases, nDOF);

      // ([P],[D])  [P] points in [D] dimensions
      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) parametric_coordinates = \n " << parametric_coordinates << std::endl;

      {
        EXCEPTWATCH;
        basis->getValues(basis_values, parametric_coordinates, m_get_derivative ?  Intrepid::OPERATOR_GRAD : Intrepid::OPERATOR_VALUE );
      }

      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) basis_values = \n " << basis_values << std::endl;

      // this function just spreads (copies) the values of the basis to all elements in the workset (numCells)
      if (m_get_derivative)
        {
          FunctionSpaceTools::HGRADtransformGRAD<double>(transformed_basis_values, JacInverse, basis_values);
        }
      else
        {
          FunctionSpaceTools::HGRADtransformVALUE<double>(transformed_basis_values, basis_values);
        }
      if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) transformed_basis_values =  " << transformed_basis_values << std::endl;

      // ([C],[P]) - place for results of evaluation
      MDArray loc_output_field_values(numCells, numInterpPoints);
      if (m_get_derivative)
        {
          loc_output_field_values.resize(numCells, numInterpPoints, spatialDim);
        }

      PerceptMesh::fillCellNodes(*get_bulk_data(), bucket_or_element, m_my_field, field_data_values_dof);

      // gather
      for (unsigned iDOF = 0; iDOF < nDOF; iDOF++)
        {
          for (int iCell = 0; iCell < numCells; iCell++)
            {
              for (int iNode = 0; iNode < numNodes; iNode++)
                {
                  field_data_values(iCell, iNode) = field_data_values_dof(iCell, iNode, iDOF);
                  if (EXTRA_PRINT_FF_HELPER) std::cout << "tmp iNode= " << iNode << "iDOF= " << iDOF << " fd= " << field_data_values(iCell, iNode) << std::endl;
                }
            }

          /// NOTE: this is needed since FunctionSpaceTools::evaluate method assumes the output array is initialized to 0
          if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) evaluate ... " << std::endl;
          loc_output_field_values.initialize(0.0);
          FunctionSpaceTools::evaluate<double>(loc_output_field_values, field_data_values, transformed_basis_values);
          if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) evaluate done " << std::endl;

          //VERIFY_OP_ON(numCells, ==, 1, "numCells...");

          if (m_get_derivative)
            {
              for (int iDim = 0; iDim < num_grad; iDim++)
                {
                  int jDim = -1;
                  if (m_deriv_spec(iDim) == "x") jDim = 0;
                  else if (m_deriv_spec(iDim) == "y") jDim = 1;
                  else if (m_deriv_spec(iDim) == "z") jDim = 2;

                  for (int iCell = 0; iCell < numCells; iCell++)
                    {
                      for (int iPoint = 0; iPoint < numInterpPoints; iPoint++)
                        {
                          // FIXME
                          // output_field_values(iPoint, iDOF, iDim) = loc_output_field_values(iCell, iPoint, jDim);
                          if (output_field_values.rank() == 2)
                            output_field_values(iPoint, iDOF*num_grad+ iDim) = loc_output_field_values(iCell, iPoint, jDim);
                          else if (output_field_values.rank() == 3)
                            output_field_values(iCell, iPoint, iDOF*num_grad+ iDim) = loc_output_field_values(iCell, iPoint, jDim);

                          if (EXTRA_PRINT_FF_HELPER) std::cout << "tmp iDOF= " << iDOF << " iDim= " << iDim <<
                            " ofd= " << output_field_values(iPoint, iDOF, iDim) << std::endl;
                        }
                    }
                }
            }
          else
            {
              for (int iCell = 0; iCell < numCells; iCell++)
                {
                  for (int iPoint = 0; iPoint < numInterpPoints; iPoint++)
                    {
                      if (output_field_values.rank() == 2)
                        output_field_values(iPoint, iDOF) = loc_output_field_values(iCell, iPoint);
                      else if (output_field_values.rank() == 3)
                        output_field_values(iCell, iPoint, iDOF) = loc_output_field_values(iCell, iPoint);
                      if (EXTRA_PRINT_FF_HELPER) std::cout << "tmp iDOF= " << iDOF << " ofd= " << output_field_values(iPoint, iDOF) << std::endl;
                    }
                }
            }

          if (EXTRA_PRINT_FF_HELPER) std::cout << "FieldFunction::operator()(elem,...) scatter done " << std::endl;
        }
    }


  }//namespace percept

#endif
