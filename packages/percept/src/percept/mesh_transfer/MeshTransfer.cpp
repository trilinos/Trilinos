// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "MeshTransfer.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>

#include <stk_mesh/base/FieldBLAS.hpp>

#include <Ioss_SubSystem.h>

#include <percept/xfer/STKMeshTransferSetup.hpp>

#include <percept/PerceptMesh.hpp>

#include "RotationTranslation.hpp"

namespace percept
{

void MeshTransfer::process_options()
{
  clp.setDocString("mesh_transfer options");

  clp.setOption("src-file",    &src_mesh,        "source file (ExodusII) containing mesh and fields to transfer." );
  clp.setOption("dst-mesh",    &dst_mesh,        "destination mesh file (ExodusII). If destination fields are found \nin this file, the resulting transferred data is added to the existing destination field data." );

  clp.setOption("dst-entity",  &dst_entity,      "destination entity type (node, element)." );
  clp.setOption("dst-name" ,   &dst_field_name,  "destination field name." );

  clp.setOption("exp-factor",  &coarse_search_expansion_factor,  "expansion factor for coarse search. Setting to a value less than one will result in no expansion. (default is 1.5)" );

  clp.setOption("src-field",   &field_name,      "source field name" );

  clp.setOption("src-rznvec",  &thvec_name,      "name of normal component (out of plane) in rz-coordinates" );
  clp.setOption("src-rzpvec",  &rzvec_name,      "name of in-plance components (2) in rz coordinates");

  clp.setOption("xrot",        &xrot,            "rotation angle (degrees) for x-axis" );
  clp.setOption("yrot",        &yrot,            "rotation angle (degrees) for y-axis" );
  clp.setOption("zrot",        &zrot,            "rotation angle (degrees) for z-axis" );

  clp.setOption("xtrans",      &xtrans,          "translation x-component" );
  clp.setOption("ytrans",      &ytrans,          "translation y-component" );
  clp.setOption("ztrans",      &ztrans,          "translation z-component" );

  clp.setOption("target",      &target_mesh,        "target file (ExodusII) containing mesh and transferred fields." );
}

void MeshTransfer::run(int argc, char** argv)
{
  process_options();

  bool found_help = false;
  if (!stk::parallel_machine_rank(comm)) {
    for (int i = 0; i < argc; ++i) {
      const std::string s(argv[i]);
      if ( s == "-h" || s == "-help" || s == "--help") { 
	found_help = true;
      }
    }
  }
  if (!stk::parallel_machine_rank(comm) && found_help) {
    clp.printHelpMessage("mesh_transfer",std::cout);
    return;
  }
  
  clp.parse( argc, argv );
  
  if (!stk::parallel_machine_rank(comm) && src_mesh=="") { 
    std::cout << "Error: must provide --src-file option. For all options run with --help" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (!stk::parallel_machine_rank(comm) && dst_mesh=="") { 
    std::cout << "Error: must provide --dst-mesh option. For all options run with --help" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (!stk::parallel_machine_rank(comm) && target_mesh=="") { 
    std::cout << "Error: must provide --target option. For all options run with --help" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (!stk::parallel_machine_rank(comm) && src_mesh==dst_mesh) { 
    std::cout << "Error: --src-file and --dst-mesh cannot be the same file: " << src_mesh << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (!stk::parallel_machine_rank(comm) && field_name.size()==0 && thvec_name.size()==0 && rzvec_name.size()==0) {
    std::cout << "Error: must specify one of --fld, --thvec or --rzvec options" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  if (!stk::parallel_machine_rank(comm) && thvec_name.size()>0 && rzvec_name.size()>0) {
    std::cout << "Error: must specify exactly one of --thvec or --rzvec options" << std::endl;
    std::exit(EXIT_FAILURE);
  }

  // initialize to the generic type
  SrcFieldType srcFieldType = SRC_FIELD;

  if (thvec_name.size()>0) { 
    field_name = thvec_name;
    srcFieldType = SRC_RZN_FIELD;
  }
  else if (rzvec_name.size()>0) {
    field_name = rzvec_name;
    srcFieldType = SRC_RZP_FIELD;
  }

  std::shared_ptr<PerceptMesh> srcMesh(new percept::PerceptMesh(3,comm));
  std::shared_ptr<PerceptMesh> dstMesh(new percept::PerceptMesh(3,comm));
  srcMesh->open(src_mesh);
  dstMesh->open(dst_mesh);

  stk::mesh::FieldBase * fromField = srcMesh->get_field(field_name);

  if (!stk::parallel_machine_rank(comm) && NULL==fromField) {
    std::cout << "Error: unknown src field: " << field_name << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  stk::mesh::EntityRank fromRank = fromField->entity_rank();

  const int toDim = (srcFieldType == SRC_FIELD) ? 
      fromField->max_size()
    : dstMesh->get_spatial_dim();

  stk::mesh::EntityRank toRank = fromRank;
  if (dst_entity=="node")
    toRank = stk::topology::NODE_RANK;
  else if (dst_entity=="element")
    toRank = stk::topology::ELEMENT_RANK;
  else
    if (!stk::parallel_machine_rank(comm))
      std::cout << "Warning: unknown dst entity type: " << dst_entity << std::endl;

  if (dst_field_name.size()==0) dst_field_name = field_name;

  stk::mesh::FieldBase * toField = NULL;
  stk::mesh::FieldBase * toField_init = dstMesh->get_fem_meta_data()->get_field(toRank, dst_field_name);
  if (toField_init == NULL) {
    toField = dstMesh->add_field(dst_field_name, toRank, toDim);
    // Q: are fields init to zero?  
  }
  else {
    // save original field contents into a temp field,
    // which we add back in after the transfer
    std::string dst_field_name_init = dst_field_name + "_init";
    toField = toField_init;
    toField_init = dstMesh->add_field(dst_field_name_init, toRank, toDim, "universal_part",false); // false = no IO
  }
  
  srcMesh->commit();
  dstMesh->commit();

  std::shared_ptr<STKMeshTransfer> mesh_transfer =
    buildSTKMeshTransfer<STKMeshTransfer>(*(srcMesh->get_bulk_data()),
			 srcMesh->get_coordinates_field(),
			 fromField,
			 *(dstMesh->get_bulk_data()),
			 dstMesh->get_coordinates_field(),
			 toField,
             "transfer",
             srcFieldType,
             coarse_search_expansion_factor);

  if (!stk::parallel_machine_rank(comm))
    std::cout << "MeshTransfer: initializing transfer" << std::endl;

  initializeSTKMeshTransfer(&*mesh_transfer);

  const int num_time_steps = srcMesh->get_database_time_step_count();
  double current_time;
  
  stk::io::StkMeshIoBroker mesh_data(comm);
  mesh_data.set_bulk_data(*(dstMesh->get_bulk_data()));

  const size_t result_output_index = mesh_data.create_output_mesh(target_mesh, stk::io::WRITE_RESULTS);
  mesh_data.add_field(result_output_index, *toField);

  const stk::mesh::FieldVector &fields = mesh_data.meta_data().get_fields();
  for (size_t i=0; i < fields.size(); i++) {
    const Ioss::Field::RoleType* role = stk::io::get_field_role(*fields[i]);
    if ( role && *role == Ioss::Field::TRANSIENT )
      {
        mesh_data.add_field(result_output_index, *fields[i]);
      }
  }

  for (int ts=1; ts<=num_time_steps; ts++) { // Exodus is 1-based
    srcMesh->read_database_at_step(ts); 
    current_time = srcMesh->get_database_time_at_step(ts);

    if (toField_init) {
      stk::mesh::field_copy(*toField, *toField_init, stk::mesh::selectField(*toField_init));
    }

    // load any existing fields from dstMesh
    if (ts <= dstMesh->get_database_time_step_count()) {
      dstMesh->read_database_at_step(ts);
    }

    if (!stk::parallel_machine_rank(comm))
      std::cout << "MeshTransfer: performing transfer at time = "
                << current_time << std::endl;

    mesh_transfer->apply();

    // add in original field values if they exist
    if (toField_init) {
      stk::mesh::field_axpy(1.0, *toField_init, *toField, stk::mesh::selectField(*toField));
    }

    // apply rotations/translations
    applyRotation(dstMesh->get_coordinates_field(), xrot, yrot, zrot);
    if (toDim==dstMesh->get_spatial_dim()) {
      applyRotation(toField, xrot, yrot, zrot);
    }

    applyTranslation(dstMesh->get_coordinates_field(), xtrans, ytrans, ztrans);

    mesh_data.process_output_request(result_output_index, current_time);
  }

  srcMesh.reset();
  dstMesh.reset();
}

} //namespace percept
