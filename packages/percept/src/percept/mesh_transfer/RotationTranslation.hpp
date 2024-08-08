// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.


#ifndef RotationTranslation_h
#define RotationTranslation_h

void applyRotation(stk::mesh::FieldBase *vectorField,
		   double&xrot, double&yrot, double&zrot) // angles in degrees
{
  stk::mesh::BulkData& bulkdata = vectorField->get_mesh();
  const stk::mesh::MetaData & meta = bulkdata.mesh_meta_data();
  stk::mesh::EntityRank rank = vectorField->entity_rank();

  const unsigned nDim = meta.spatial_dimension();

  const double xtheta = M_PI*xrot/180.0;
  const double ytheta = M_PI*yrot/180.0;
  const double ztheta = M_PI*zrot/180.0;

  double xrot_mat[3][3] = {{1,0,           0},
			   {0,cos(xtheta),-sin(xtheta)},
			   {0,sin(xtheta), cos(xtheta)}};

  double yrot_mat[3][3] = {{ cos(ytheta), 0, sin(ytheta)},
			   { 0,           1, 0},
			   {-sin(ytheta), 0, cos(ytheta)}};

  double zrot_mat[3][3] = {{cos(ztheta),-sin(ztheta), 0},
			   {sin(ztheta), cos(ztheta), 0},
			   {0,           0,           1}};

  double yxrot_mat[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for (unsigned i=0; i<3; i++) {
    for (unsigned j=0; j<3; j++) {
      for (unsigned k=0; k<3; k++) {
	yxrot_mat[i][j] += yrot_mat[i][k]*xrot_mat[k][j];
      }
    }
  }
  double zyxrot_mat[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
  for (unsigned i=0; i<3; i++) {
    for (unsigned j=0; j<3; j++) {
      for (unsigned k=0; k<3; k++) {
	zyxrot_mat[i][j] += zrot_mat[i][k]*yxrot_mat[k][j];
      }
    }
  }

  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::BucketVector const& entity_buckets =
    bulkdata.get_buckets( rank, select_used );

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;

    double * v = static_cast<double*>(stk::mesh::field_data(*vectorField, b));

    if (nDim == 3) {
      double tmp[3];
      const stk::mesh::Bucket::size_type length = b.size();

      for (unsigned ie=0; ie<length; ie++) {
        tmp[0] = zyxrot_mat[0][0]*v[0]+zyxrot_mat[0][1]*v[1]+zyxrot_mat[0][2]*v[2];
        tmp[1] = zyxrot_mat[1][0]*v[0]+zyxrot_mat[1][1]*v[1]+zyxrot_mat[1][2]*v[2];
        tmp[2] = zyxrot_mat[2][0]*v[0]+zyxrot_mat[2][1]*v[1]+zyxrot_mat[2][2]*v[2];

        v[0] = tmp[0];
        v[1] = tmp[1];
        v[2] = tmp[2];

        v += nDim;
      }
    }
    else if (nDim == 2) {
      double tmp[2];
      const stk::mesh::Bucket::size_type length = b.size();

      for (unsigned ie=0; ie<length; ie++) {
        tmp[0] = zyxrot_mat[0][0]*v[0]+zyxrot_mat[0][1]*v[1];
        tmp[1] = zyxrot_mat[1][0]*v[0]+zyxrot_mat[1][1]*v[1];

        v[0] = tmp[0];
        v[1] = tmp[1];

        v += nDim;
      }
    }
  }
}

void applyTranslation(stk::mesh::FieldBase *vectorField,
		      double&xtrans, double&ytrans, double&ztrans)
{
  stk::mesh::BulkData& bulkdata = vectorField->get_mesh();
  const stk::mesh::MetaData & meta = bulkdata.mesh_meta_data();
  stk::mesh::EntityRank rank = vectorField->entity_rank();

  stk::mesh::Selector select_used =
    meta.locally_owned_part() |
    meta.globally_shared_part() ;

  stk::mesh::BucketVector const& entity_buckets =
    bulkdata.get_buckets( rank, select_used );

  const unsigned nDim = meta.spatial_dimension();

  for ( stk::mesh::BucketVector::const_iterator ib = entity_buckets.begin();
        ib != entity_buckets.end() ; ++ib ) {

    stk::mesh::Bucket & b = **ib ;

    double * v = static_cast<double*>(stk::mesh::field_data(*vectorField, b));

    const stk::mesh::Bucket::size_type length = b.size();
    for (unsigned ie=0; ie<length; ie++) {
      v[0] += xtrans;
      v[1] += ytrans;
      if (nDim == 3) {
        v[2] += ztrans;
      }
      
      v += nDim;
    }
  }
}


#endif
