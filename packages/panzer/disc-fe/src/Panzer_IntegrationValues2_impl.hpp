#ifndef PANZER_INTEGRATION_VALUES_2_IMPL_HPP
#define PANZER_INTEGRATION_VALUES_2_IMPL_HPP

/* Normally, we would fold this definition into the .cpp file since we
   use ETI to compile this object. However, for unit testing the stand
   alone device function swapQuadraturePoints, we need the definition
   at the calling point to avoid requiring relocatable device code for
   cuda. So for the IntegrationValues2 object, if there is a
   standalone device member function, put it in this file (for unit
   testing).
 */
namespace panzer {
  
template <typename Scalar>
void IntegrationValues2<Scalar>::
swapQuadraturePoints(int cell,
                     int a,
                     int b) const
{
  const int new_cell_point = a;
  const int old_cell_point = b;

  const int cell_dim = ref_ip_coordinates.extent(2);

#ifdef PANZER_DEBUG
  KOKKOS_ASSERT(cell < ip_coordinates.extent_int(0));
  KOKKOS_ASSERT(a < ip_coordinates.extent_int(1));
  KOKKOS_ASSERT(b < ip_coordinates.extent_int(1));
  KOKKOS_ASSERT(cell >= 0);
  KOKKOS_ASSERT(a >= 0);
  KOKKOS_ASSERT(b >= 0);
#endif

  // Using scratch_for_compute_side_measure instead of hold. This
  // works around UVM issues of DFAD temporaries in device code.
  // Scalar hold;

  scratch_for_compute_side_measure(0) = weighted_measure(cell,new_cell_point);
  weighted_measure(cell,new_cell_point) = weighted_measure(cell,old_cell_point);
  weighted_measure(cell,old_cell_point) = scratch_for_compute_side_measure(0);

  scratch_for_compute_side_measure(0) = jac_det(cell,new_cell_point);
  jac_det(cell,new_cell_point) = jac_det(cell,old_cell_point);
  jac_det(cell,old_cell_point) = scratch_for_compute_side_measure(0);

  for(int dim=0;dim<cell_dim;++dim){

    scratch_for_compute_side_measure(0) = ref_ip_coordinates(cell,new_cell_point,dim);
    ref_ip_coordinates(cell,new_cell_point,dim) = ref_ip_coordinates(cell,old_cell_point,dim);
    ref_ip_coordinates(cell,old_cell_point,dim) = scratch_for_compute_side_measure(0);

    scratch_for_compute_side_measure(0) = ip_coordinates(cell,new_cell_point,dim);
    ip_coordinates(cell,new_cell_point,dim) = ip_coordinates(cell,old_cell_point,dim);
    ip_coordinates(cell,old_cell_point,dim) = scratch_for_compute_side_measure(0);

    scratch_for_compute_side_measure(0) = surface_normals(cell,new_cell_point,dim);
    surface_normals(cell,new_cell_point,dim) = surface_normals(cell,old_cell_point,dim);
    surface_normals(cell,old_cell_point,dim) = scratch_for_compute_side_measure(0);

    for(int dim2=0;dim2<cell_dim;++dim2){

      scratch_for_compute_side_measure(0) = jac(cell,new_cell_point,dim,dim2);
      jac(cell,new_cell_point,dim,dim2) = jac(cell,old_cell_point,dim,dim2);
      jac(cell,old_cell_point,dim,dim2) = scratch_for_compute_side_measure(0);

      scratch_for_compute_side_measure(0) = jac_inv(cell,new_cell_point,dim,dim2);
      jac_inv(cell,new_cell_point,dim,dim2) = jac_inv(cell,old_cell_point,dim,dim2);
      jac_inv(cell,old_cell_point,dim,dim2) = scratch_for_compute_side_measure(0);
    }
  }

  // Rotation matrices are always in 3D
  for(int dim=0; dim<3; ++dim){
    for(int dim2=0; dim2<3; ++dim2){
      scratch_for_compute_side_measure(0) = surface_rotation_matrices(cell,new_cell_point,dim,dim2);
      surface_rotation_matrices(cell,new_cell_point,dim,dim2) = surface_rotation_matrices(cell,old_cell_point,dim,dim2);
      surface_rotation_matrices(cell,old_cell_point,dim,dim2) = scratch_for_compute_side_measure(0);
    }
  }
}

}

#endif
