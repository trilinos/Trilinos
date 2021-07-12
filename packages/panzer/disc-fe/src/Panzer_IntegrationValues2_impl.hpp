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

template <typename Scalar>
void IntegrationValues2<Scalar>::
convertNormalToRotationMatrix(const Scalar normal[3], Scalar transverse[3], Scalar binormal[3])
{
  using T = Scalar;
  
  const T n  = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

  // If this fails then the geometry for this cell is probably undefined
  if(n > 0.){
    // Make sure transverse is not parallel to normal within some margin of error
    transverse[0]=0.;transverse[1]=1.;transverse[2]=0.;
    if(std::fabs(normal[0]*transverse[0]+normal[1]*transverse[1])>0.9){
      transverse[0]=1.;transverse[1]=0.;
    }

    const T nt = normal[0]*transverse[0]+normal[1]*transverse[1]+normal[2]*transverse[2];

    // Note normal has unit length
    const T mult = nt/(n*n); // = nt

    // Remove normal projection from transverse
    for(int dim=0;dim<3;++dim){
      transverse[dim] = transverse[dim] - mult * normal[dim];
    }

    const T t = sqrt(transverse[0]*transverse[0]+transverse[1]*transverse[1]+transverse[2]*transverse[2]);
    KOKKOS_ASSERT(t != 0.);
    for(int dim=0;dim<3;++dim){
      transverse[dim] /= t;
    }

    // We assume a right handed system such that b = n \times t
    binormal[0] = (normal[1] * transverse[2] - normal[2] * transverse[1]);
    binormal[1] = (normal[2] * transverse[0] - normal[0] * transverse[2]);
    binormal[2] = (normal[0] * transverse[1] - normal[1] * transverse[0]);

    // Normalize binormal
    const T b = sqrt(binormal[0]*binormal[0]+binormal[1]*binormal[1]+binormal[2]*binormal[2]);
    for(int dim=0;dim<3;++dim){
      binormal[dim] /= b;
    }
  } else {
    transverse[0] = 0.;
    transverse[1] = 0.;
    transverse[2] = 0.;
    binormal[0] = 0.;
    binormal[1] = 0.;
    binormal[2] = 0.;
  }
}

}

#endif
