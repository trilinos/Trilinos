// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/FieldTypes.hpp>

#include <percept/function/FieldFunction.hpp>

#include <percept/function/internal/SimpleSearcher.hpp>
#include <percept/function/internal/STKSearcher.hpp>
#include <percept/function/internal/STKSearcherDef.hpp>
#include <percept/ParallelUtil.hpp>
#include <percept/norm/Norm.hpp>

  namespace percept
  {
    //========================================================================================================================
    // high-level interface
    FieldFunction::FieldFunction(const char *name, stk::mesh::FieldBase *field, PerceptMesh& mesh,
                  int domain_dimension,
                  int codomain_dimension,
                  SearchType searchType,
                  unsigned integration_order) :

      Function(name, Dimensions(domain_dimension), Dimensions(codomain_dimension), integration_order),
      m_my_field(field), m_bulkData(mesh.get_bulk_data()),
      m_cachedElement(), m_searcher(0),
      m_cached_topo_key(0), m_cached_basis(0), m_searchType(searchType), m_get_derivative(false)
                                //, m_parallelEval(true)
    {
    }

    bool FieldFunction::m_parallelEval=true;

    stk::mesh::FieldBase *FieldFunction::get_field() {return m_my_field; }

    void FieldFunction::interpolateFrom(Function& function)
    {
      EXCEPTWATCH;
      stk::mesh::MetaData& metaData = stk::mesh::MetaData::get(*m_my_field);
      stk::mesh::BulkData& bulkData = *m_bulkData;

      PerceptMesh meshUtil(&metaData, &bulkData);

      /// for each node in the codomain, evaluate the function_to_interpolate's function, assign to the codomain field
      nodalOpLoop(*meshUtil.get_bulk_data(), function, m_my_field);

    }//interpolateFrom


    //========================================================================================================================
    // low-level interface

    FieldFunction::FieldFunction(const char *name, stk::mesh::FieldBase *field, stk::mesh::BulkData *bulk,
                                 Dimensions domain_dimensions,
                                 Dimensions codomain_dimensions,
                                 SearchType searchType,
                                 unsigned integration_order) :
      Function(name, domain_dimensions, codomain_dimensions, integration_order),
      m_my_field(field), m_bulkData(bulk), m_cachedElement(), m_searcher(0),
      m_cached_topo_key(0), m_cached_basis(0), m_searchType(searchType), m_get_derivative(false)
                                //, m_parallelEval(true)
    {
      const stk::mesh::MetaData& metaData = bulk->mesh_meta_data();
      m_coordinatesField = metaData.get_field<double>(stk::topology::NODE_RANK, "coordinates");
    }

    FieldFunction::FieldFunction(const char *name, stk::mesh::FieldBase *field, PerceptMesh& eMesh,
                                 Dimensions domain_dimensions,
                                 Dimensions codomain_dimensions,
                                 SearchType searchType,
                                 unsigned integration_order) :
      Function(name, domain_dimensions, codomain_dimensions, integration_order),
      m_my_field(field), m_bulkData(eMesh.get_bulk_data()), m_cachedElement(), m_searcher(0),
      m_cached_topo_key(0), m_cached_basis(0), m_searchType(searchType), m_get_derivative(false)
    {

      const stk::mesh::MetaData& metaData = m_bulkData->mesh_meta_data();
      m_coordinatesField = metaData.get_field<double>(stk::topology::NODE_RANK, "coordinates");

    }

    FieldFunction::~FieldFunction()
    {
      if (m_searcher)
        delete m_searcher;
    }

    stk::mesh::BulkData *FieldFunction::get_bulk_data() {return m_bulkData; }


    void FieldFunction::localEvaluation(MDArray& input_phy_points, MDArray& output_field_values, double time)
    {
      m_parallelEval = false;
      (*this)(input_phy_points, output_field_values, time);
      m_parallelEval = true;
    }

    void FieldFunction::setup_searcher(int D_)
    {
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
                if (D_ == 3)
                  m_searcher = new STKSearcher(m_bulkData);
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
    }

    /** Evaluate the function at this input point (or points) returning value(s) in output_field_values
     *
     *   In the following, the arrays are dimensioned using the notation (from Intrepid2's doc):
     *
     *   [C]         - num. integration domains (cells/elements)
     *   [F]         - num. Intrepid2 "fields" (number of bases within an element == num. nodes typically)
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
     *  [R] is used for the rank of MDArray's
     */
    void FieldFunction::operator()(MDArray& input_phy_points, MDArray& output_field_values, double time)
    {
      EXCEPTWATCH;
      argsAreValid(input_phy_points, output_field_values);
      if (NormBase::m_is_norm_being_evaluated)
        throw std::runtime_error("can't use FieldFunction as sub-function of other functions with Norm not set to use TURBO_BUCKET or TURBO_ELEMENT");

      m_found_on_local_owned_part = false;

      //// single point only (for now)
      unsigned found_it = 0;

      int D_ = last_dimension(input_phy_points);
      MDArray found_parametric_coordinates_one("found_parametric_coordinates_one", 1, 1, D_);
      setup_searcher(D_);

      MDArray output_field_values_local("output_field_values_local",output_field_values.layout());
      Kokkos::deep_copy(output_field_values_local,output_field_values);
      //int R_output = output_field_values.rank();

      int R_input = input_phy_points.rank();
      int P_ = (R_input == 1 ? 1 : input_phy_points.extent_int(R_input-2));

      // FIXME for tensor valued fields
      //int DOF_ = last_dimension(output_field_values_local);

      MDArray input_phy_points_one("input_phy_points_one",1,1,D_);
      //MDArray output_field_values_one("output_field_values_one",1,DOF_);

      int C_ = 1;
      if (R_input == 3)
        {
          C_ = input_phy_points.extent_int(0);
        }
      for (int iC = 0; iC < C_; iC++)
        {
          for (int iP = 0; iP < P_; iP++)
            {
              for (int iD = 0; iD < D_; iD++)
                {
                  switch(R_input)
                    {
                    case 1: input_phy_points_one(0, 0, iD) = input_phy_points(iD); break;
                    case 2: input_phy_points_one(0, 0, iD) = input_phy_points(iP, iD); break;
                    case 3: input_phy_points_one(0, 0, iD) = input_phy_points(iC, iP, iD); break;
                    default: VERIFY_1("bad rank");
                    }
                }

              stk::mesh::Entity found_element = stk::mesh::Entity();
              {
                EXCEPTWATCH;
                found_element = m_searcher->findElement(input_phy_points_one, found_parametric_coordinates_one, found_it, m_cachedElement);
              }

              // if found element on the local owned part, evaluate
              if (found_it)
                {
                  m_found_on_local_owned_part = true;
                  if (( EXTRA_PRINT) && m_searchType==STK_SEARCH)
                    std::cout << "FieldFunction::operator() found element # = " << m_bulkData->identifier(found_element) << std::endl;

                  //auto found_parametric_coordinates = Kokkos::subview(found_parametric_coordinates_one, 0, Kokkos::ALL(), Kokkos::ALL());

                  MDArray found_parametric_coordinates("found_parametric_coordinates", 1, D_);
                  Kokkos::deep_copy(found_parametric_coordinates,
                                    Kokkos::subview(found_parametric_coordinates_one, 0, Kokkos::ALL(), Kokkos::ALL()));

                  (*this)(input_phy_points, output_field_values_local, found_element, found_parametric_coordinates);

                  /*
                  for (int iDOF = 0; iDOF < DOF_; iDOF++)
                    {
                      switch (R_output)
                        {
                        case 1: output_field_values_local( iDOF)        = output_field_values_one(0, iDOF); break;
                        case 2: output_field_values_local(iP, iDOF)     = output_field_values_one(0, iDOF); break;
                        case 3: output_field_values_local(iC, iP, iDOF) = output_field_values_one(0, iDOF); break;
                        default: VERIFY_1("bad rank");
                        }
                    }
                  */
                }
              else
                {
                  if (!m_parallelEval)
                    {
                      std::cout << "P[" << Util::get_rank() << "] FieldFunction::operator() found_it = " << found_it << " points= "
                                << printContainer(input_phy_points_one)
                                << std::endl;

                      throw std::runtime_error("FieldFunction::operator() in local eval mode and didn't find element - logic error");
                    }
                  double max_val = std::numeric_limits<double>::max();
                  Kokkos::deep_copy(output_field_values_local,max_val);
                }

              // make sure it is found somewhere
              if (m_parallelEval)
                {
                  all_reduce( m_bulkData->parallel() , stk::ReduceMax<1>( & found_it ) );
                }

              if (EXTRA_PRINT) std::cout << "FieldFunction::operator() global found_it = " << found_it << std::endl;

              if (!found_it)
                {
                  throw std::runtime_error("FieldFunction::operator() couldn't find element");
                }

              if (m_parallelEval)
                {
                  percept_global_lex_min( m_bulkData->parallel(), output_field_values.size(), &output_field_values_local[0], &output_field_values[0]);
                }
              else
                {
                  output_field_values=output_field_values_local;
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
                                   const stk::mesh::Entity element, const MDArray& parametric_coordinates, double time_value_optional)
    {
      EXCEPTWATCH;
      helper(*m_bulkData, input_phy_points, output_field_values, element, parametric_coordinates, time_value_optional);
    }

    void FieldFunction::operator()(MDArray& input_phy_points, MDArray& output_field_values,
                                   const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates, double time_value_optional)
    {
      EXCEPTWATCH;
#ifndef NDEBUG
      int num_elements_in_bucket   = bucket.size();
      VERIFY_OP(input_phy_points.extent_int(0), ==, num_elements_in_bucket, "FieldFunction::operator() mismatch in input_phy_points and num_elements_in_bucket");
      VERIFY_OP(output_field_values.extent_int(0), ==, num_elements_in_bucket, "FieldFunction::operator() mismatch in input_phy_points and num_elements_in_bucket");
#endif
      helper(*m_bulkData, input_phy_points, output_field_values, bucket, parametric_coordinates, time_value_optional);
    }


  }//namespace percept
