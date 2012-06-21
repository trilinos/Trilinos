
#include <stk_percept/function/FieldFunction.hpp>

#include <stk_percept/function/internal/SimpleSearcher.hpp>
#include <stk_percept/function/internal/STKSearcher.hpp>
#include <stk_percept/function/internal/STKSearcherDef.hpp>
#include <stk_percept/function/internal/BuildBoundingBoxesDef.hpp>
#include <stk_percept/ParallelUtil.hpp>

namespace stk
{
  namespace percept
  {
    //========================================================================================================================
    // high-level interface
    FieldFunction::FieldFunction(const char *name, mesh::FieldBase *field, PerceptMesh& mesh,
                  int domain_dimension,
                  int codomain_dimension,
                  SearchType searchType,
                  unsigned integration_order) :

      Function(name, Dimensions(domain_dimension), Dimensions(codomain_dimension), integration_order),
      m_my_field(field), m_bulkData(mesh.get_bulk_data()),
      m_cachedElement(0), m_searcher(0),
      m_cached_topo_key(0), m_cached_basis(0), m_searchType(searchType), m_get_derivative(false)
                                //, m_parallelEval(true)
    {
    }

    bool FieldFunction::m_parallelEval=true;

    mesh::FieldBase *FieldFunction::get_field() {return m_my_field; }

    void FieldFunction::interpolateFrom(Function& function)
    {
      EXCEPTWATCH;
      mesh::fem::FEMMetaData& metaData = stk::mesh::fem::FEMMetaData::get(*m_my_field);
      mesh::BulkData& bulkData = *m_bulkData;

      PerceptMesh meshUtil(&metaData, &bulkData);

      /// for each node in the codomain, evaluate the function_to_interpolate's function, assign to the codomain field
      meshUtil.nodalOpLoop(function, m_my_field);

    }//interpolateFrom


    //========================================================================================================================
    // low-level interface

    typedef mesh::Field<double>                     ScalarFieldType ;
    typedef mesh::Field<double, mesh::Cartesian>    VectorFieldType ;


    FieldFunction::FieldFunction(const char *name, mesh::FieldBase *field, mesh::BulkData *bulk,
                                 Dimensions domain_dimensions,
                                 Dimensions codomain_dimensions,
                                 SearchType searchType,
                                 unsigned integration_order) :
      Function(name, domain_dimensions, codomain_dimensions, integration_order),
      m_my_field(field), m_bulkData(bulk), m_cachedElement(0), m_searcher(0),
      m_cached_topo_key(0), m_cached_basis(0), m_searchType(searchType), m_get_derivative(false)
                                //, m_parallelEval(true)
    {
      mesh::fem::FEMMetaData& metaData = stk::mesh::fem::FEMMetaData::get(*bulk);
      m_coordinatesField = metaData.get_field<VectorFieldType >("coordinates");
    }

    FieldFunction::FieldFunction(const char *name, mesh::FieldBase *field, PerceptMesh& eMesh,
                                 Dimensions domain_dimensions,
                                 Dimensions codomain_dimensions,
                                 SearchType searchType,
                                 unsigned integration_order) :
      Function(name, domain_dimensions, codomain_dimensions, integration_order),
      m_my_field(field), m_bulkData(eMesh.get_bulk_data()), m_cachedElement(0), m_searcher(0),
      m_cached_topo_key(0), m_cached_basis(0), m_searchType(searchType), m_get_derivative(false)
    {

      mesh::fem::FEMMetaData& metaData = stk::mesh::fem::FEMMetaData::get(*m_bulkData);
      m_coordinatesField = metaData.get_field<VectorFieldType >("coordinates");

    }

    FieldFunction::~FieldFunction()
    {
      if (m_searcher)
        delete m_searcher;
    }

    mesh::BulkData *FieldFunction::get_bulk_data() {return m_bulkData; }



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
     *  Dimensions of input_phy_points are required to be either ([D]) or ([P],[D])
     *  Dimensions of output_field_values are required to be ([DOF]) or ([P],[DOF]) respectively
     *
     */

    void FieldFunction::localEvaluation(MDArray& input_phy_points, MDArray& output_field_values, double time)
    {
      m_parallelEval = false;
      (*this)(input_phy_points, output_field_values, time);
      m_parallelEval = true;
    }
    void FieldFunction::operator()(MDArray& input_phy_points, MDArray& output_field_values, double time)
    {
      EXCEPTWATCH;
      argsAreValid(input_phy_points, output_field_values);

      m_found_on_local_owned_part = false;

      //// single point only (for now)
      unsigned found_it = 0;

      int nDim = last_dimension(input_phy_points);
      int spatialDim = nDim;
      static MDArray found_parametric_coordinates_one(1, spatialDim);

      MDArray output_field_values_local = output_field_values;

      //FIXME
#if 1
      if (!m_searcher)
        {
          switch (m_searchType)
            {
            case SIMPLE_SEARCH:
              m_searcher = new SimpleSearcher(m_bulkData);
              break;
            case STK_SEARCH:
              {
                //int spDim = last_dimension(input_phy_points);
                if (spatialDim == 3)
                  m_searcher = new STKSearcher<3>(m_bulkData);
                else
                  {
                    //m_searcher = new STKSearcher<2>(this);
                    throw std::runtime_error("STK_SEARCH not ready for 2D, use SIMPLE_SEARCH");
                  }
              }
              break;
            default:
              throw std::runtime_error("FieldFunction::operator() unknown search type");
              break;
            }
          //std::cout << "setupSearch..." << std::endl;
          m_searcher->setupSearch();
          //std::cout << "setupSearch...done" << std::endl;
        }
#endif

      int numInputPoints = 1;
      int rankInput = input_phy_points.rank();
      switch(rankInput)
        {
        case 1:
          numInputPoints = 1;
          break;
        case 2:
          numInputPoints = input_phy_points.dimension(0);
          break;
        case 3:
          numInputPoints = input_phy_points.dimension(1);
          break;
        }


      // FIXME for tensor valued fields
      int nDOF = last_dimension(output_field_values_local);

      static MDArray input_phy_points_one(1,3);
      static MDArray output_field_values_one(1,3);
      if (nDOF != output_field_values_one.dimension(0))
        {
          output_field_values_one.resize(1,nDOF);
        }
      if (spatialDim != input_phy_points_one.dimension(1))
        {
          input_phy_points_one.resize(1, spatialDim);
        }

      int iStart = 0;
      int iEnd = iStart + numInputPoints;
      int nCells = 1;
      if (rankInput == 3)
        {
          nCells = input_phy_points.dimension(0);
        }
      for (int iCell = 0; iCell < nCells; iCell++)
        {
          for (int iPoint = iStart; iPoint < iEnd; iPoint++)
            {
              for (int iDim = 0; iDim < spatialDim; iDim++)
                {
                  EXCEPTWATCH;
                  if (rankInput == 2)
                    {
                      EXCEPTWATCH;
                      input_phy_points_one(0, iDim) = input_phy_points(iPoint, iDim);
                    }
                  else if (rankInput == 1)
                    {
                      EXCEPTWATCH;
                      input_phy_points_one(0, iDim) = input_phy_points(iDim);
                    }
                  else if (rankInput == 3)
                    {
                      EXCEPTWATCH;
                      input_phy_points_one(0, iDim) = input_phy_points(iCell, iPoint, iDim);
                    }
                }

              const stk::mesh::Entity *found_element = 0;
              {
                EXCEPTWATCH;
                //if (m_searchType==STK_SEARCH) std::cout << "find" << std::endl;

                found_element = m_searcher->findElement(input_phy_points_one, found_parametric_coordinates_one, found_it, m_cachedElement);
                //if (m_searchType==STK_SEARCH)                std::cout << "find..done found_it=" << found_it << std::endl;
              }

              if (0 || EXTRA_PRINT)
                {
                  if (Util::getFlag(9828)) std::cout << "P[" << Util::get_rank() << "] FieldFunction::operator() found_it = " << found_it << std::endl;
                }

              // if found element on the local owned part, evaluate
              if (found_it)
                {
                  m_found_on_local_owned_part = true;
                  if (( EXTRA_PRINT) && m_searchType==STK_SEARCH)
                    std::cout << "FieldFunction::operator() found element # = " << found_element->identifier() << std::endl;

                  (*this)(input_phy_points_one, output_field_values_one, *found_element, found_parametric_coordinates_one);

                  for (int iDOF = 0; iDOF < nDOF; iDOF++)
                    {
                      if (output_field_values_local.rank() == 1)
                        {
                          output_field_values_local( iDOF) = output_field_values_one(0, iDOF);
                        }
                      else if (output_field_values_local.rank() == 2)
                        {
                          output_field_values_local(iPoint, iDOF) = output_field_values_one(0, iDOF);
                        }
                      else if (output_field_values_local.rank() == 3)
                        {
                          output_field_values_local(iCell, iPoint, iDOF) = output_field_values_one(0, iDOF);
                        }
                    }
                }
              else
                {
                  if (!m_parallelEval)
                    {
                      std::cout << "P[" << Util::get_rank() << "] FieldFunction::operator() found_it = " << found_it << " points= "
                                << input_phy_points_one
                                << std::endl;

                      throw std::runtime_error("FieldFunction::operator() in local eval mode and didn't find element - logic error");
                    }
                  double max_val = std::numeric_limits<double>::max();
#if 0
                  for (int iDOF = 0; iDOF < nDOF; iDOF++)
                    {
                      if (output_field_values_local.rank() == 1)
                        {
                          output_field_values_local( iDOF) = max_val;
                        }
                      else if (output_field_values_local.rank() == 2)
                        {
                          output_field_values_local(iPoint, iDOF) = max_val;
                        }
                      else if (output_field_values_local.rank() == 3)
                        {
                          output_field_values_local(iCell, iPoint, iDOF) = max_val;
                        }
                    }
#else
                  output_field_values_local.initialize(max_val);
#endif
                }

              // make sure it is found somewhere
              if (m_parallelEval)
                {
                  all_reduce( m_bulkData->parallel() , ReduceMax<1>( & found_it ) );
                }

              if (EXTRA_PRINT)       std::cout << "FieldFunction::operator() global found_it = " << found_it << std::endl;

              if (!found_it)
                {
                  throw std::runtime_error("FieldFunction::operator() couldn't find element");
                }

              if (m_parallelEval)
                {
                  //void stk_percept_global_lex_min(stk::ParallelMachine comm,  int n , T local_min[] , T global_min[] );

                  stk_percept_global_lex_min( m_bulkData->parallel(), output_field_values.size(), &output_field_values_local[0], &output_field_values[0]);
#if 0
                  if (Util::getFlag(9828)) std::cout <<  "P[" << Util::get_rank() << "] ofv.size= " << output_field_values.size() << std::endl;
                  if (Util::getFlag(9828)) std::cout <<  "P[" << Util::get_rank() << "] ofv_l.size= " << output_field_values_local.size() << std::endl;
                  if (1 && Util::getFlag(9828)) std::cout <<  "P[" << Util::get_rank() << "] ofv= "
                                                     << output_field_values[0] << " "
                                                     << output_field_values[1] << " "
                                                     << output_field_values[2] << " "
                                                     << " ofv_local= "
                                                     << output_field_values_local[0] << " "
                                                     << output_field_values_local[1] << " "
                                                     << output_field_values_local[2] << " "
                                                     << std::endl;
#endif
                }
              else
                {
                  output_field_values = output_field_values_local;
                }
              m_cachedElement = found_element;
            }
        }

    }

    /** Evaluate the function on this element at the parametric coordinates and return in output_field_values.
     *
     *  Dimensions of parametric_coordinates are required to be ([P],[D])
     *  Dimensions of output_field_values are required to be ([P],[DOF]), (or in future, ([C],[P],[DOF]) )
     *
     */

    void FieldFunction::operator()(MDArray& input_phy_points, MDArray& output_field_values,
                                   const stk::mesh::Entity& element, const MDArray& parametric_coordinates, double time_value_optional)
    {
      EXCEPTWATCH;
      helper(input_phy_points, output_field_values, element, parametric_coordinates, time_value_optional);
    }

    void FieldFunction::operator()(MDArray& input_phy_points, MDArray& output_field_values,
                                   const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates, double time_value_optional)
    {
      EXCEPTWATCH;
#ifndef NDEBUG
      int num_elements_in_bucket   = bucket.size();
      VERIFY_OP(input_phy_points.dimension(0), ==, num_elements_in_bucket, "FieldFunction::operator() mismatch in input_phy_points and num_elements_in_bucket");
      VERIFY_OP(output_field_values.dimension(0), ==, num_elements_in_bucket, "FieldFunction::operator() mismatch in input_phy_points and num_elements_in_bucket");
#endif
      helper(input_phy_points, output_field_values, bucket, parametric_coordinates, time_value_optional);
    }


  }//namespace percept
}//namespace stk
