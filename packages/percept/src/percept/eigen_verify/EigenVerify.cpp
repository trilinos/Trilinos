// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "EigenVerify.hpp"

#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_io/IossBridge.hpp>

#include <Ioss_SubSystem.h>

#include <percept/xfer/STKMeshTransferSetup.hpp>

#include <percept/function/FieldFunction.hpp>
#include <percept/norm/Norm.hpp>
#include <percept/FieldBLAS.hpp>
#include <percept/PerceptMesh.hpp>

namespace percept
{

void EigenVerify::process_options()
{
  clp.setDocString("eigen_verify options");

  clp.setOption("inmesh1",     &meshes_in[0],        "mesh file #1 (ExodusII)." );
  clp.setOption("inmesh2",     &meshes_in[1],        "mesh file #2 (ExodusII)." );

  clp.setOption("outmesh",  &mesh_out,            "out mesh file (ExodusII)." );
  clp.setOption("outfile",  &file_out,            "out results file (text)." );

  clp.setOption("fieldname", &field_name,        "name of nodal vector field to transfer." );
}

void EigenVerify::create_mesh_data(stk::io::StkMeshIoBroker * mesh_data,
				   const std::string &filename)
{
  mesh_data->property_add(Ioss::Property("FIELD_SUFFIX_SEPARATOR", ""));
  mesh_data->add_mesh_database(filename, "exodus", stk::io::READ_MESH);
  mesh_data->create_input_mesh();

  // This defines all fields found on the input mesh as stk fields
  mesh_data->add_all_mesh_fields_as_input_fields();
}

void EigenVerify::load_time_data(const int m)
{
  std::shared_ptr<Ioss::Region> ioss_region = mesh_data[m]->get_input_ioss_region();

  const int numTimeSteps =
    ioss_region->get_property("state_count").get_int();

  num_time_steps_all[m] = numTimeSteps;

  time_steps_all[m].resize(numTimeSteps);

  for (int i=0; i<numTimeSteps; ++i) {
    time_steps_all[m][i]=ioss_region->get_state_time(i+1);
  }
}

void EigenVerify::create_fields(const int num_time_steps)
{
  const std::string field_name_all      = field_name + ".all";
  const std::string xfer_field_name_all = field_name + ".xfer.all";

  for (int m=0; m<num_meshes; m++) {

    // allocate data for the eigenvectors (sending)
    fieldAll[m] = & (mesh_data[m]->meta_data().declare_field<double>(stk::topology::NODE_RANK, field_name_all));

    stk::mesh::put_field_on_mesh( *fieldAll[m],
			  mesh_data[m]->meta_data().universal_part(),
			  mesh_data[m]->meta_data().spatial_dimension(),
			  num_time_steps,
        nullptr);
    stk::io::set_field_output_type(*fieldAll[m], stk::io::FieldOutputType::VECTOR_3D);

    // allocate data for the eigenvectors (receiving)
    xferFieldAll[m] = & (mesh_data[m]->meta_data().declare_field<double>(stk::topology::NODE_RANK, xfer_field_name_all));

    stk::mesh::put_field_on_mesh( *xferFieldAll[m],
			  mesh_data[m]->meta_data().universal_part(),
			  mesh_data[m]->meta_data().spatial_dimension(),
			  num_time_steps,
        nullptr);

    // get eigenvector field
    inputField[m] = mesh_data[m]->meta_data().get_field<double>(stk::topology::NODE_RANK, field_name);

    if (NULL==inputField[m]) {
      std::cout << "Error: unknown displacement field: " << field_name 
		<< " on mesh " << meshes_in[m] << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  // create fields to store nodal errors
  std::string error_field_name = "error." + field_name;
  errorField = & ( mesh_data[1]->meta_data().declare_field<double>(stk::topology::NODE_RANK, error_field_name) );
  stk::mesh::put_field_on_mesh( *errorField, mesh_data[1]->meta_data().universal_part(), mesh_data[1]->meta_data().spatial_dimension(), nullptr);
  stk::io::set_field_output_type(*errorField, stk::io::FieldOutputType::VECTOR_3D);
}

void EigenVerify::load_field_data(const int m)
{
  for (int ts=0; ts<num_time_steps_all[m]; ts++) {

    mesh_data[m]->read_defined_input_fields(ts+1); // 1-based step

    field_copy_component(mesh_data[m]->bulk_data(), inputField[m], fieldAll[m], ts);
  }
}

int EigenVerify::get_num_time_steps()
{
  bool diff_num_time_steps = false;

  for (int m=1; m<num_meshes; m++) {
    if (num_time_steps_all[m] != num_time_steps_all[0])
      diff_num_time_steps = true;
  }

  if (diff_num_time_steps)
    std::exit(EXIT_FAILURE);

  return num_time_steps_all[0];
}

void compute_field_error(
  const stk::mesh::BulkData &bulkdata,
  const int m,
  stk::mesh::Field<double> * field,
  const int index_m,
  const double sign_factor,
  stk::mesh::Field<double> * xferField,
  stk::mesh::Field<double> * errorField)
{
  const stk::mesh::MetaData & meta = bulkdata.mesh_meta_data();

  const unsigned field_size = field->max_size();
  const unsigned nDim = meta.spatial_dimension();

  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::BucketVector const& entity_buckets =
    bulkdata.get_buckets( stk::topology::NODE_RANK, select_used );

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;

    const double * field_data = stk::mesh::field_data(*field, b);
    const double * xfer_field_data = stk::mesh::field_data(*xferField, b);
    double * error_data = stk::mesh::field_data(*errorField, b);

    const stk::mesh::Bucket::size_type length = b.size();

    unsigned offset_field, offset_xfer_field;

    for ( stk::mesh::Bucket::size_type k = 0 ; k < length ; ++k ) {

      offset_field = k*field_size + m*nDim;
      offset_xfer_field = k*field_size + index_m*nDim;

      for (unsigned idim=0; idim<nDim; idim++) {
	error_data[k*nDim+idim] = field_data[offset_field+idim]
	  - sign_factor * xfer_field_data[offset_xfer_field+idim];
      }
    }
  }
}

void
compute_permutations_and_signs(Kokkos::View<double**,Kokkos::HostSpace> & mac_values,
			       std::vector<int> & reorder_indices,
			       std::vector<double> & sign_factor)
{
  const int num_time_steps = sign_factor.size();

  for (int m=0; m<num_time_steps; m++) {
    int best_index = 0;
    double best_mac = fabs(mac_values(m,0));

    for (int n=1; n<num_time_steps; n++) {
      if (fabs(mac_values(m,n)) > best_mac) {
	best_index = n;
	best_mac = fabs(mac_values(m,n));
      }
    }
    reorder_indices[m] = best_index;
    sign_factor[m] = fabs(mac_values(m,best_index)) / mac_values(m,best_index);
  }
}

void
write_mac(Kokkos::View<double**,Kokkos::HostSpace> & mac_values,
	  std::vector<int> & reorder_indices,
	  std::vector<double> & sign_factor)
{
  const int num_time_steps = sign_factor.size();

  // DEBUG
  std::ofstream outmac("mac.dat");
  outmac.precision(3);
  for (int m=0; m<num_time_steps; m++) {
    for (int n=0; n<num_time_steps; n++) {
      outmac << std::fixed << m << " " << n << " " << mac_values(m,n) << std::endl;
    }
  }

  for (int m=0; m<num_time_steps; m++) {
    // DEBUG
    outmac << m << " "
	   << reorder_indices[m] << " "
	   << sign_factor[m] << std::endl;
  }
}

void EigenVerify::run(int argc, char** argv)
{
  process_options();
  clp.parse( argc, argv );

  for (int m=0; m<num_meshes; m++) {

    mesh_data.push_back(new stk::io::StkMeshIoBroker(comm));

    create_mesh_data(mesh_data[m], meshes_in[m]);

    load_time_data(m);
  }

  const int num_time_steps = get_num_time_steps();

  create_fields(num_time_steps);

  for (int m=0; m<num_meshes; m++) {

    mesh_data[m]->populate_bulk_data();

    load_field_data(m);
  }

  // init output mesh for fine mesh only
  const size_t result_output_index = mesh_data.back()->create_output_mesh(mesh_out, stk::io::WRITE_RESULTS);

  mesh_data.back()->add_field(result_output_index, *errorField);

  // build transfer: 0 => 1
  stk::mesh::Field<double> * coordinates_from = mesh_data[0]->meta_data().get_field<double>(stk::topology::NODE_RANK, "coordinates");
  stk::mesh::Field<double> * coordinates_to = mesh_data[1]->meta_data().get_field<double>(stk::topology::NODE_RANK, "coordinates");

  std::shared_ptr<STKMeshTransfer> mesh_transfer_01 =
    buildSTKMeshTransfer<STKMeshTransfer>(mesh_data[0]->bulk_data(),
		      coordinates_from,
		      fieldAll[0],
		      mesh_data[1]->bulk_data(),
		      coordinates_to,
		      xferFieldAll[1],
		      "transfer_01");

  initializeSTKMeshTransfer(&*mesh_transfer_01);

  mesh_transfer_01->apply();

  // fix up signs and ordering
  Kokkos::View<double**, Kokkos::HostSpace> mac_values("mac_values",num_time_steps, num_time_steps);

  field_compute_MAC(mesh_data[1]->bulk_data(), mesh_data[1]->meta_data(),
		    fieldAll[1], xferFieldAll[1], mac_values);

  std::vector<int> reorder_indices(num_time_steps);
  std::vector<double> sign_factor(num_time_steps);

  compute_permutations_and_signs(mac_values, reorder_indices, sign_factor);

  write_mac(mac_values, reorder_indices, sign_factor);

  const std::vector<double> & output_times = time_steps_all.back();

  std::ofstream outfile(file_out.c_str());

  for (int m=0; m<num_time_steps; m++) {
    // copy scaled/permuted evecs to error vector
    const int index_m = reorder_indices[m];
    const double sign_factor_m = sign_factor[m];

    compute_field_error(mesh_data[1]->bulk_data(),
			m,
			fieldAll[1],
			index_m,
			sign_factor_m,
			xferFieldAll[1],
			errorField);

    // compute L2 norm of error vector
    FieldFunction errorFunc("error", errorField, &(mesh_data[1]->bulk_data()),
			    Dimensions(3), Dimensions(3));

    stk::mesh::Selector element_selector = PerceptMesh::get_selector_of_rank(mesh_data[1]->meta_data(), stk::topology::ELEMENT_RANK);

    Norm<2> norm(mesh_data[1]->bulk_data(), &element_selector);
    ConstantFunction result(0.0, "result");

    // TODO add as global var to Exo output
    norm(errorFunc, result);

    outfile << output_times[m] << " " << result.getValue() << std::endl;

    // output error vector (for current grid time step)
    mesh_data.back()->process_output_request(result_output_index, output_times[m]);
  }
  delete mesh_data[0];
  delete mesh_data[1];
}

} //namespace percept
