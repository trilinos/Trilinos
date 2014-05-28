/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <assert.h>

#include <stk_util/parallel/Parallel.hpp>
#include <init/Ionit_Initializer.h>
#include <Ioss_SubSystem.h>

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>

#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/fem/TopologyDimensions.hpp>

#include <stk_io/IossBridge.hpp>

/** \addtogroup stk_io_module
 * \{
 */

/**
 * Example code showing a basic, but complete, mesh to results output
 * coding including subsetting and periodic field input and output.
 * Includes handling of nodeblocks, element blocks, nodesets,
 * and sidesets.  Attribute fields and distribution factor
 * fields are also supported.
 *
 * This example can serve as the basis for adding binary IO support to
 * an application.  The code here uses the Ioss to/from stk::mesh
 * bridge functions in the stk::io namespace defined in IossBridge.hpp
 * include file.
 */
namespace stk {
  namespace percept {
    namespace io_util {

      /// Declare "coordinates" field and put it on the universal part. This
      /// example also defines all Ioss::Field::TRANSIENT fields that exist on the
      /// Ioss::Nodeblock as fields on the universal part.
      void process_read_nodeblocks_meta    (Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta, int& spatial_dim);

      /// Declare a part for each element block on the Ioss::Region
      /// 'region' unless the element block has the "omitted" property set
      /// to the value 1. The example then iterates each element block and
      /// defines any Ioss::Field::ATTRIBUTE and Ioss::Field::TRANSIENT fields that exist on the
      /// Ioss::ElementBlock as fields on the corresponding part.
      void process_read_elementblocks_meta (Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta);

      /// Declare a part for each Ioss::NodeSet on the Ioss::Region
      /// 'region' unless the nodeset has the "omitted" property set
      /// to the value 1. The example then iterates each nodeset and
      /// defines any "distribution factor" and Ioss::Field::TRANSIENT fields that
      /// exist on the Ioss::NodeSet as fields on the corresponding
      /// part.
      void process_read_nodesets_meta      (Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta);

      /// Declare a part for each Ioss::SideSet on the Ioss::Region
      /// 'region' unless the sideset has the "omitted" property set
      /// to the value 1. The example then iterates each sideset and
      /// defines any "distribution factor" and Ioss::Field::TRANSIENT fields that
      /// exist on the Ioss::SideSet as fields on the corresponding
      /// part.
      ///
      /// Each sideblock in the active sidesets is then processed by
      /// defining a part for each Ioss::SideBlock on the Ioss::SideSet
      /// unless the sideblock has the "omitted" property set to the value
      /// 1. The example then iterates each sideblock and defines any
      /// "distribution factor" and Ioss::Field::TRANSIENT fields that exist on the
      /// Ioss::SideBlock as fields on the corresponding part.
      void process_read_sidesets_meta      (Ioss::Region &region, stk::mesh::fem::FEMMetaData &meta);

      /// NOTE: This must be called after the process_read_elementblocks() call
      /// since there may be nodes that exist in the database that are
      /// not part of the analysis mesh due to subsetting of the element
      /// blocks.
      ///
      /// Populates  the "coordinates" field for all active nodes in the model.
      void process_read_nodeblocks_bulk    (Ioss::Region &region, stk::mesh::BulkData &bulk);

      /// NOTE: This should be the first function called of any of the
      /// "process_read_X" type functions that take an stk::mesh::BulkData
      /// argument, especially if the input Ioss::Region mesh is going to
      /// be subsetted (have element blocks omitted).
      ///
      /// This function iterates all non-omitted element blocks and
      /// declares each element (and the corresponding nodes) in the
      /// element block. If there are any Ioss::Field::ATTRIBUTE fields on the element
      /// block (for example, shell thickness or particle radius), then
      /// that field data is alse read and the corresponding
      /// stk::mesh::Field populated.
      void process_read_elementblocks_bulk (Ioss::Region &region, stk::mesh::BulkData &bulk);

      /// Iterates each non-omitted Ioss::NodeSet and then iterates each
      /// node in the Ioss::NodeSet.  If the node exists (that is, it is
      /// connected to a non-omitted Ioss::ElementBlock), then that node
      /// is associated with the part corresponding to this
      /// Ioss::NodeSet. If the "distribution_factor" field exists, then
      /// that data is also associated with the field.
      void process_read_nodesets_bulk      (Ioss::Region &region, stk::mesh::BulkData &bulk);

      /// Process each non-omitted Ioss::SideSet and the contained
      /// non-omitted Ioss::SideBlock and associate each element-side pair with
      /// the corresponding part if the underlying element is active.  If
      /// the "distribution_factor" field exists, then that data is also
      /// associated with the corresponding field.
      void process_read_sidesets_bulk      (Ioss::Region &region, stk::mesh::BulkData &bulk);

      /// A minimal example function showing how field data on the
      /// Ioss::Region entities can be periodically transferred to the
      /// corresponding field(s) on the stk::mesh entities. This would be
      /// used to bring in initial condition data or interpolation data or
      /// any other scenario in which data on the mesh file needs to be
      /// transferred to the stk::mesh fields.
      void process_read_input_request (Ioss::Region &region, stk::mesh::BulkData &bulk, int step);

      /// A minimal example function showing how stk::mesh field data can
      /// periodically be output to a results, history, heartbeat, or
      /// restart database.  The scheduling would be done either in this
      /// function or at a higher level and is not shown here. The
      /// function iterates all parts and if there is a corresponding Ioss
      /// part on the Ioss::Region, all fields defined to be output are
      /// iterated and their data output to the corresponding
      /// Ioss::Field. The function calls the
      /// stk::io::is_valid_part_field() function to determine whether the
      /// field should be output and then calls the
      /// stk::io::field_data_to_ioss() function to do the actual output
      /// of the field.
      void process_output_request(Ioss::Region &region, stk::mesh::BulkData &bulk, int step);

      /// This function shows the basic calls needed to perform definition
      /// and input of the mesh model and definition and periodic output
      /// of a results database. The function is given the mesh filename
      /// and the output filename and goes through all steps of
      /// associating the filename with an Ioss::DatabaseIO object of the
      /// correct type ("exodusII" in this example); creating an
      /// Ioss::Region and then defining an stk::mesh corresponding to
      /// this mesh.  The function also provides an example of how
      /// specific element blocks existing in the mesh database could be
      /// omitted from the analysis model.
      ///
      /// The example then shows how to define a results database
      /// corresponding to the analysis model and periodically output the
      /// results in an execute loop.
      ///
      /// A true application would have to provide additional
      /// functionality and robustness, but the example shows how the
      /// basic functionality can be provided by an application.
      ///
      /// Note that the paradigm illustrated here is different than the
      /// mesh input and output paradigm provided in the current
      /// framework.  In this case, the application is responsible for the
      /// majority of the IO behavior and the toolkit only provides some
      /// helper functions to bridge between the Ioss and the stk::mesh.
      /// It is hoped that this paradigm will result in more functionality
      /// for the application with less complication and overhead.
    }
  }
}
