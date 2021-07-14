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

template<typename Scalar,
         typename T0,typename T1,typename T2,typename T3,
         typename T4,typename T5,typename T6,typename T7>
void
swapQuadraturePoints(int cell,
                     int a,
                     int b,
                     T0& ref_ip_coordinates,
                     T1& ip_coordinates,
                     T2& weighted_measure,
                     T3& jac,
                     T4& jac_det,
                     T5& jac_inv,
                     T6& surface_normals,
                     T7& surface_rotation_matrices)
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

  // If this is a DFAD type, we will need to fix allocation to size the derivative array. 
  Scalar hold;

  hold = weighted_measure(cell,new_cell_point);
  weighted_measure(cell,new_cell_point) = weighted_measure(cell,old_cell_point);
  weighted_measure(cell,old_cell_point) = hold;

  hold = jac_det(cell,new_cell_point);
  jac_det(cell,new_cell_point) = jac_det(cell,old_cell_point);
  jac_det(cell,old_cell_point) = hold;

  for(int dim=0;dim<cell_dim;++dim){

    hold = ref_ip_coordinates(cell,new_cell_point,dim);
    ref_ip_coordinates(cell,new_cell_point,dim) = ref_ip_coordinates(cell,old_cell_point,dim);
    ref_ip_coordinates(cell,old_cell_point,dim) = hold;

    hold = ip_coordinates(cell,new_cell_point,dim);
    ip_coordinates(cell,new_cell_point,dim) = ip_coordinates(cell,old_cell_point,dim);
    ip_coordinates(cell,old_cell_point,dim) = hold;

    hold = surface_normals(cell,new_cell_point,dim);
    surface_normals(cell,new_cell_point,dim) = surface_normals(cell,old_cell_point,dim);
    surface_normals(cell,old_cell_point,dim) = hold;

    for(int dim2=0;dim2<cell_dim;++dim2){

      hold = jac(cell,new_cell_point,dim,dim2);
      jac(cell,new_cell_point,dim,dim2) = jac(cell,old_cell_point,dim,dim2);
      jac(cell,old_cell_point,dim,dim2) = hold;

      hold = jac_inv(cell,new_cell_point,dim,dim2);
      jac_inv(cell,new_cell_point,dim,dim2) = jac_inv(cell,old_cell_point,dim,dim2);
      jac_inv(cell,old_cell_point,dim,dim2) = hold;
    }
  }

  // Rotation matrices are always in 3D
  for(int dim=0; dim<3; ++dim){
    for(int dim2=0; dim2<3; ++dim2){
      hold = surface_rotation_matrices(cell,new_cell_point,dim,dim2);
      surface_rotation_matrices(cell,new_cell_point,dim,dim2) = surface_rotation_matrices(cell,old_cell_point,dim,dim2);
      surface_rotation_matrices(cell,old_cell_point,dim,dim2) = hold;
    }
  }
}

}

#endif
