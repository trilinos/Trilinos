// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CONVERT_NORMAL_TO_ROTATION_MATRIX_HPP
#define PANZER_CONVERT_NORMAL_TO_ROTATION_MATRIX_HPP

namespace panzer {

/** \brief Builds an orthonormal (transverse, binormal) pair completing a
  * right-handed rotation matrix for the given unit normal.
  *
  * Given a surface normal, constructs a transverse vector orthogonal to
  * it and a binormal vector such that (transverse, binormal, normal)
  * form a right-handed orthonormal basis, i.e. binormal = normal x
  * transverse. If the normal has zero length, transverse and binormal
  * are both set to zero.
  *
  * \tparam Scalar the scalar type of the vector components.
  * \param[in] normal the unit surface normal, as a 3-component array.
  * \param[out] transverse the computed transverse vector, as a 3-component array.
  * \param[out] binormal the computed binormal vector, as a 3-component array.
  */
template <typename Scalar>
KOKKOS_INLINE_FUNCTION
void
convertNormalToRotationMatrix(const Scalar normal[3], Scalar transverse[3], Scalar binormal[3])
{
  using T = Scalar;
  
  const T n  = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);

  // If this fails then the geometry for this cell is probably undefined
  if(n > 0.){
    // Make sure transverse is not parallel to normal within some margin of error
    transverse[0]=0.;transverse[1]=1.;transverse[2]=0.;
    if(Kokkos::fabs(normal[0]*transverse[0]+normal[1]*transverse[1])>0.9){
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
