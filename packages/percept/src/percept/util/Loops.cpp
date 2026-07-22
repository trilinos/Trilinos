// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <percept/util/Loops.hpp>

#include <percept/ExceptionWatch.hpp>
#include <percept/function/FieldFunction.hpp>
#include <percept/FieldTypes.hpp>

namespace percept {

  /** \brief Loop over all buckets and apply \param bucketOp passing in the argument \param field to \param bucketOp */
  void bucketOpLoop(stk::mesh::BulkData& bulkData,
                    BucketOp& bucketOp,
                    stk::mesh::FieldBase *field,
                    stk::mesh::Part *part)
  {
    EXCEPTWATCH;

    if (part) {
      stk::mesh::Selector selector(*part);
      bucketOpLoop(bulkData, bucketOp, field, &selector);
    }
    else {
      bucketOpLoop(bulkData, bucketOp, field, (stk::mesh::Selector *)0);
    }
  }

  void bucketOpLoop(stk::mesh::BulkData& bulkData,
                    BucketOp& bucketOp,
                    stk::mesh::FieldBase *field,
                    stk::mesh::Selector *selector,
                    bool is_surface_norm)
  {
    EXCEPTWATCH;

    stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK;
    if (is_surface_norm) rank = bulkData.mesh_meta_data().side_rank();

    const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( rank );

    for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k ) {
      if (!selector || (*selector)(**k)) { // this is where we do part selection
        const stk::mesh::Bucket & bucket = **k ;
        bool breakLoop = bucketOp(bucket, field, bulkData);

        if (breakLoop) return;
      }
    }
  }

  /** \brief Loop over all elements and apply \param elementOp passing in the argument \param field to \param elementOp */
  void elementOpLoop(stk::mesh::BulkData& bulkData,
                     ElementOp& elementOp,
                     stk::mesh::FieldBase *field,
                     stk::mesh::Part *part)
  {
    EXCEPTWATCH;

    if (part)
      {
        stk::mesh::Selector selector(*part);
        elementOpLoop(bulkData, elementOp, field, &selector);
      }
    else
      {
        elementOpLoop(bulkData, elementOp, field, (stk::mesh::Selector *)0);
      }
  }

  /** \brief Loop over all elements and apply \param elementOp passing in the argument \param field to \param elementOp */
  void elementOpLoop(stk::mesh::BulkData& bulkData,
                     ElementOp& elementOp,
                     stk::mesh::FieldBase *field,
                     stk::mesh::Selector *selector,
                     bool is_surface_norm)
  {
    EXCEPTWATCH;
    //checkState("elementOpLoop");
    elementOp.init_elementOp();

    // FIXME consider caching the coords_field in FieldFunction
    //CoordinatesFieldType *coords_field = metaData.get_field<CoordinatesFieldType::value_type >(stk::topology::NODE_RANK, "coordinates");
    stk::mesh::EntityRank rank = stk::topology::ELEMENT_RANK;
    if (is_surface_norm) rank = bulkData.mesh_meta_data().side_rank();

    const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( rank );
    for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (!selector || (*selector)(**k))  // this is where we do part selection
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_elements_in_bucket   = bucket.size();

            // FIXME for multiple points
            for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
              {
                stk::mesh::Entity element = bucket[iElement];

                bool breakLoop = elementOp(element, field, bulkData);
                //std::cout << "PerceptMesh::elementOpLoop breakLoop= " << breakLoop << std::endl;
                if (breakLoop)
                  {
                    elementOp.fini_elementOp();
                    return;
                  }

              }

          }
      }
    elementOp.fini_elementOp();
  }

  void nodalOpLoop(stk::mesh::BulkData& bulkData,
                   GenericFunction& nodalOp,
                   stk::mesh::FieldBase *field,
                   stk::mesh::Selector* selector)
  {
#if STK_PERCEPT_LITE
    VERIFY_MSG("not available in PerceptMeshLite");
#else
    EXCEPTWATCH;

    stk::mesh::FieldBase *coords_field = bulkData.mesh_meta_data().get_field<double>(stk::topology::NODE_RANK, "coordinates");

    // for each node in the codomain, evaluate the function_to_interpolate's function, assign to the codomain field

    const std::vector<stk::mesh::Bucket*> & buckets = bulkData.buckets( stk::topology::NODE_RANK );

    int num_nodes = 0;

    for ( std::vector<stk::mesh::Bucket*>::const_iterator k = buckets.begin() ; k != buckets.end() ; ++k )
      {
        if (!selector || (*selector)(**k))  // this is where we do part selection
          {
            stk::mesh::Bucket & bucket = **k ;
            const unsigned num_nodes_in_bucket   = bucket.size();

            //unsigned spatialDim = 0;
            //double * coord = stk::mesh::field_data( *coords_field , bucket.begin() );
            double * coord = static_cast<double*>(stk::mesh::field_data(*coords_field, bucket));
            //if (Util::getFlag(9829)) std::cout << "spatialDim= " << spatialDim << std::endl;

            unsigned stride = 0;
            double * output_nodal_field = 0;

            if (field) {
              stride = field->max_size();
              output_nodal_field = (double*)stk::mesh::field_data(*field, bucket);
            }
            else
              stride = (nodalOp.getCodomainDimensions().size() ? nodalOp.getCodomainDimensions()[0] : 0);

            //int inDim = nodalOp.getDomainDimensions()[0];
            //int outDim = nodalOp.getCodomainDimensions()[0];
            int inDim = bulkData.mesh_meta_data().spatial_dimension();

            num_nodes += num_nodes_in_bucket;
            // FIXME for multiple points
            for (unsigned inode = 0; inode < num_nodes_in_bucket; inode++)
              {
                MDArray pt("pt", inDim);  // FIXME for spatialDim
                for (int iSpace = 0; iSpace < inDim; iSpace++)
                  {
                    pt(iSpace) = coord[iSpace];
                  }
                MDArray out("out", stride);

                // an optional setting of the codomain from existing values (allows for +=, etc.)
                // if(set_output) {
                if (field)
                  {
                    for (unsigned jout = 0; jout < stride; jout++)
                      {
                        out(jout) = output_nodal_field[jout];
                        //if (Util::getFlag(9829)) std::cout << "bef jout= " << jout << " val= " << out(jout) << std::endl;
                      }
                  }

                //if (Util::getFlag(9829)) std::cout << "nodalOp= " << nodalOp << std::endl;
                {
                  FieldFunction::m_parallelEval=false;
                  nodalOp(pt, out);
                  FieldFunction::m_parallelEval=true;
                }

                if (field)
                  {
                    for (unsigned jout = 0; jout < stride; jout++)
                      {
                        //if (Util::getFlag(9829)) std::cout << "aft jout= " << jout << " val= " << out(jout) << std::endl;
                        output_nodal_field[jout] = out(jout);
                      }
                  }

                if (field) output_nodal_field += stride;  // FIXME
                coord += inDim;  // FIXME
              }

          }
      }

    if (0) std::cout << "P[" << stk::parallel_machine_rank(bulkData.parallel()) << "] num_nodes= "<< num_nodes << std::endl;

#endif
  }

} // namespace percept
