// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_ContourSubElement.hpp>
#include <Akri_DiagWriter.hpp>
#include <Akri_LevelSet.hpp>
#include <Akri_MathUtil.hpp>
#include <Akri_MeshHelpers.hpp>
#include <Akri_Utility.hpp>

#include <stk_util/parallel/ParallelComm.hpp>

#include <cmath>

#include <Akri_MasterElementDeterminer.hpp>

namespace krino{

bool ContourSubElement::is_more(const stk::math::Vector3d & x, const stk::math::Vector3d & y)
{
  // note that in the case of x==y, this will return false
  if (utility::is_more(x[0],y[0]) || (!utility::is_less(x[0],y[0]) &&
      (utility::is_more(x[1],y[1]) || (!utility::is_less(x[1],y[1]) &&
       (utility::is_more(x[2],y[2]))))))
    {
      return true;
    }
  return false;
}

ContourSubElement::ContourSubElement( const ContourElement *in_owner,
			const int in_subelement_depth,
			const int subelement_sign )
    : my_owner( in_owner ),
      my_subelement_depth( in_subelement_depth ),
      my_sign( subelement_sign )
{
}

int
ContourSubElement::build_facets( FacetedSurfaceBase & facets )
{
  int start_size = facets.size();

  if ( !my_subelements.empty() )
  {
    for ( auto && subelem : my_subelements )
    {
      subelem->build_facets( facets );
    }
  }
  else
  {
    const int numSides = get_num_sides();
    const int * sideIds = get_side_ids();
    for ( int iside = 0; iside < numSides; iside++ )
    {
      // side_ids == -2 indicates that sides is on interface
      if ( sideIds[iside] == -2 )
      {
        side_facets( facets, iside );
      }
    }
  }

  int added_facets = facets.size() - start_size;
  return added_facets;
}

void
ContourSubElement::dump_structure() const
{

  if ( !my_subelements.empty() )
  {
    for ( int i = 0; i < my_subelement_depth; i++ )
    {
      krinolog << "  ";
    }
    krinolog << "my subelement type = " << topology().name() << ", my # of subelements = " << my_subelements.size() << std::endl;
    for ( auto && subelem : my_subelements )
    {
      subelem->dump_structure();
    }
  }
  else
  {
    for ( int i = 0; i < my_subelement_depth; i++ )
    {
      krinolog << "  ";
    }
    krinolog << "my subelement type = " << topology().name() << std::endl;
  }
}

void
ContourSubElement::dump_details() const
{
  if ( !my_subelements.empty() )
  {
    krinolog << "my subelement type = " << topology().name() << ", my # of subelements = " << my_subelements.size() << std::endl;
    for ( auto && subelem : my_subelements )
    {
      subelem->dump_details();
    }
  }
  else
  {
    krinolog << "--------------------begin subelement definition--------------------" << std::endl;
    krinolog << *this;
    krinolog << "---------------------end subelement definition---------------------" << std::endl;
  }
}

int
ContourSubElement::side_facets( FacetedSurfaceBase & facets, int side ) const
{
  const std::string & owner_type = my_owner->dist_topology().name();
  const std::string & sub_type = topology().name();
  ThrowRuntimeError("Subelement decomposition for subelement of type '" << sub_type
      << "' which was generated from owning element of type '" << owner_type
      << "' is missing the capability to generate conformal facets.");
  return -1;
}

double
ContourSubElement::side_area( int side ) const
{
  const std::string & owner_type = my_owner->dist_topology().name();
  const std::string & sub_type = topology().name();
  ThrowRuntimeError("Subelement decomposition for subelement of type '" << sub_type
      << "' which was generated from owning element of type '" << owner_type
      << "' is missing the capability to compute its side area.");
  return 0.;
}

std::ostream &
ContourSubElement::put( std::ostream& os ) const
{
  os << "Subelement description:" << std::endl;
  os << "  type = " << topology().name()
     << ", relative volume = " << relative_volume()
     << ", parametric_quality = " << parametric_quality()
     << ", physical_quality = " << physical_quality() << std::endl;
  const int numNodes = get_num_nodes();
  const stk::math::Vector3d * nodeCoords = get_coordinates_at_nodes();
  const double * nodeDist = get_distance_at_nodes();
  const int * sideIds = get_side_ids();
  for ( int i = 0; i < numNodes; i++ )
    {
      stk::math::Vector3d x = my_owner->coordinates( nodeCoords[i] );
      os << "  coords[" << i << "] = ("
      << nodeCoords[i][0] << ","
      << nodeCoords[i][1] << ","
      << nodeCoords[i][2] << ")"
      << ",  x = ("
      << x[0] << ","
      << x[1] << ","
      << x[2] << ")"
      << ", dist = " << nodeDist[i]
      << ", sign = " << -1 + 2*LevelSet::sign_change(nodeDist[i],-1.) << std::endl;
    }
  for ( int i = 0; i < get_num_sides(); i++ )
    {
      os << "  side_ids[" << i << "] = " << sideIds[i]
	 << ", side_quality[" << i << "] = " << side_quality(i) << std::endl;
    }
  // matlab visualization
  os << "  matlabvertices = [";
  for ( int i = 0; i < numNodes; i++ )
    {
      os << nodeCoords[i][0] << " "
	 << nodeCoords[i][1] << " "
	 << nodeCoords[i][2] << "; ";
    }
  os << "];" << std::endl;
  os << "  physical space matlabvertices = [";
  for ( int i = 0; i < numNodes; i++ )
    {
      stk::math::Vector3d x = my_owner->coordinates( nodeCoords[i] );
      os << x[0] << " "
	 << x[1] << " "
	 << x[2] << "; ";
    }
  os << "];" << std::endl;

  return os ;
}

double
ContourSubElement::relative_volume() const
{
  if (topology() == stk::topology::TETRAHEDRON_4)
    return compute_tet_volume(get_coordinates_at_nodes());
  if (topology() == stk::topology::TRIANGLE_3_2D)
    return compute_tri_volume(get_coordinates_at_nodes());

  // This is a relative volume compared to the owner volume.
  // Actually this is a relative volume if the "parametric" volume of the element is unity.
  // Otherwise, it is off by a factor.
  const int nelem = 1;
  const int dim   = spatial_dim();
  const auto & masterElement = get_master_element();
  const int numNodes = get_num_nodes();
  const stk::math::Vector3d * nodeCoords = get_coordinates_at_nodes();
  const int nint  = masterElement.num_intg_pts();
  std::vector<double> coords(numNodes * dim, 0.);
  std::vector<double> det_J(nint, 0.);
  double error = 0.;

  // integration weights
  const double * intg_weights = masterElement.intg_weights();

  // load coords
  int count = 0;
  for ( int i = 0; i < numNodes; i++ )
    {
      for ( int j = 0; j < dim; j++ )
	{
	  coords[count++] = nodeCoords[i][j];
	}
    }

  // determinant at integration points
  masterElement.determinant( dim, nelem, coords.data(), det_J.data(), &error );

  double elem_volume = 0.;
  for ( int ip = 0; ip < nint; ip++ )
    {
      elem_volume += det_J[ip] * intg_weights[ip];
    }

  return elem_volume;
}

bool
ContourSubElement::have_interface_sides() const
{
  if ( !my_subelements.empty())
  {
    for ( auto && subelem : my_subelements )
    {
      if (subelem->have_interface_sides())
      {
        return true;
      }
    }
  }
  else
  {
    const int numSides = get_num_sides();
    const int * sideIds = get_side_ids();

    for ( int iside = 0; iside < numSides; iside++ )
    {
      // side_ids == -2 indicates that sides is on interface
      if ( sideIds[iside] == -2 )
      {
        return true;
      }
    }
  }
  return false;
}

ContourSubElement::~ContourSubElement()
{
  for ( auto && subelem : my_subelements )
    delete subelem;
}

double ContourSubElement::compute_area_of_interface() const
{
  double surfaceArea = 0.;
  if ( !my_subelements.empty() )
  {
    for ( auto && subelem : my_subelements )
      surfaceArea += subelem->compute_area_of_interface( );
    return surfaceArea;
  }
  else
  {
    const int numSides = get_num_sides();
    const int * sideIds = get_side_ids();
    for ( int iside = 0; iside < numSides; iside++ )
    {
      // side_ids == -2 indicates that sides is on interface
      if ( sideIds[iside] == -2 )
      {
        surfaceArea += side_area( iside );
      }
    }
  }
  return surfaceArea;
}

double ContourSubElement::compute_relative_signed_volume(const int signOfDomain) const
{
  STK_ThrowAssert(signOfDomain == -1 || signOfDomain == 1);
  if ( !my_subelements.empty() )
  {
    double signedVolume = 0.;
    for ( auto && subelem : my_subelements )
      signedVolume += subelem->compute_relative_signed_volume(signOfDomain);
    return signedVolume;
  }

  if (my_sign == signOfDomain)
    return relative_volume();
  return 0.;
}


int
ContourSubElement::num_intg_pts(const int intg_pt_sign)
{
  int num_pts = 0;

  if ( !my_subelements.empty() )
    {
      for ( auto && subelem : my_subelements )
	{
	  num_pts += subelem->num_intg_pts(intg_pt_sign);
	}
    }
  else
    {
      if ( 0 == intg_pt_sign ) // interface points
	{
          const auto & sideMasterElement = get_side_master_element();
          const int numSides = get_num_sides();
          const int * sideIds = get_side_ids();
	  for ( int iside = 0; iside < numSides; iside++ )
	    {
	      // side_ids == -2 indicates that side is on interface
	      if ( sideIds[iside] == -2 )
		{
		  num_pts += sideMasterElement.num_intg_pts();
		}
	    }
	}
      else // volume points
	{
	  if ( intg_pt_sign != my_sign )
	    {
	      return(0);
	    }
	  else
	    {
	      return ( get_master_element().num_intg_pts() );
	    }
	}
    }

  return num_pts;
}

int
ContourSubElement::gather_intg_pts( const int intg_pt_sign,
			     sierra::ArrayContainer<double,DIM,NINT> & intg_pt_locations,
			     sierra::ArrayContainer<double,NINT> & intg_weights,
			     sierra::ArrayContainer<double,NINT> & determinant,
			     int index )
{
  if ( !my_subelements.empty() )
    {
      for ( auto && subelem : my_subelements )
	{
	  index = subelem->gather_intg_pts( intg_pt_sign,
              intg_pt_locations,
              intg_weights,
              determinant,
              index );
	}
    }
  else
    {
      const int numNodes = get_num_nodes();
      const stk::math::Vector3d * nodeCoords = get_coordinates_at_nodes();

      if ( 0 == intg_pt_sign ) // interface points
	{
          const auto & sideMasterElement = get_side_master_element();
          const int numSides = get_num_sides();
          const int * sideIds = get_side_ids();

	  for ( int iside = 0; iside < numSides; iside++ )
	    {
	      // side_ids == -2 indicates that side is on interface
	      if ( sideIds[iside] == -2 )
		{
		  const int nelem = 1;
		  const int dim   = spatial_dim();
		  const stk::topology Top = topology();
		  const stk::topology sideTop = Top.side_topology(iside);
		  const int side_num_intg_pts = sideMasterElement.num_intg_pts();
		  const int side_num_nodes = sideTop.num_nodes();
		  const unsigned * const lnn = get_side_node_ordinals(Top, iside );
		  double error;

		  // temp arrays
		  sierra::ArrayContainer<double,DIM,NPE_COORD> coords(dim,side_num_nodes);
		  sierra::ArrayContainer<double,NINT> det_J(side_num_intg_pts);

		  // load coords
		  for ( int i = 0; i < side_num_nodes; i++ )
		    {
		      stk::math::Vector3d coordinates = my_owner->coordinates( nodeCoords[lnn[i]] );
		      for ( int d = 0; d < dim; d++ ) coords(d,i) = coordinates[d];
		    }

		  // determinant at integration points
		  sideMasterElement.determinant( dim, nelem, coords.ptr(), det_J.ptr(), &error );

		  // integration weights
		  const double * intg_wts_ptr = sideMasterElement.intg_weights();
		  const sierra::Array<const double,NINT> intg_wts(intg_wts_ptr,side_num_intg_pts);

		  // basis fns at integration point locations
		  const double * bf_ptr = sideMasterElement.shape_fcn();
		  const sierra::Array<const double,NPE_VAR,NINT> bf(bf_ptr,side_num_nodes,side_num_intg_pts);

		  for ( int ip = 0; ip < side_num_intg_pts; ++ip)
		    {
		      determinant(index) = det_J(ip);
		      intg_weights(index) = intg_wts(ip);

		      stk::math::Vector3d xi(stk::math::Vector3d::ZERO);
		      for ( int i = 0; i < side_num_nodes; i++ )
			xi += bf(i,ip) * nodeCoords[lnn[i]];

		      for ( int d = 0; d < dim; ++d )
			intg_pt_locations(d,index) = xi[d];

		      index++;
		    }
		}
	    }
	}
      else // volume points
	{
	  STK_ThrowAssert(-1 == intg_pt_sign || 1 == intg_pt_sign);

	  const auto & masterElement = get_master_element();

	  if ( intg_pt_sign != my_sign )
	    {
	      return(index);
	    }

	  const int nelem            = 1;
	  const int dim              = spatial_dim();
	  const int vol_num_intg_pts = masterElement.num_intg_pts();

	  // temp arrays
	  sierra::ArrayContainer<double,DIM,NPE_COORD> coords(dim,numNodes);
	  sierra::ArrayContainer<double,NINT> det_J(vol_num_intg_pts);

	  double error;

	  // integration weights
	  const double * intg_wts_ptr = masterElement.intg_weights();
	  const sierra::Array<const double,NINT> intg_wts(intg_wts_ptr,vol_num_intg_pts);

	  // load coords
	  for ( int i = 0; i < numNodes; i++ )
	    {
	      for ( int d = 0; d < dim; d++ ) coords(d,i) = nodeCoords[i][d];
	    }

	  // determinant at integration points
	  masterElement.determinant( dim, nelem, coords.ptr(), det_J.ptr(), &error );

	  // basis fns at integration point locations
	  const double * bf_ptr = masterElement.shape_fcn();
	  const sierra::Array<const double,NPE_VAR,NINT> bf(bf_ptr,numNodes,vol_num_intg_pts);

	  for ( int ip = 0; ip < vol_num_intg_pts; ++ip)
	    {
	      determinant(index) = det_J(ip);
	      intg_weights(index) = intg_wts(ip);

	      stk::math::Vector3d xi(stk::math::Vector3d::ZERO);
	      for ( int i = 0; i < numNodes; i++ )
		xi += bf(i,ip) * nodeCoords[i];

	      for ( int d = 0; d < dim; d++ )
		intg_pt_locations(d,index) = xi[d];

	      index++;
	    }
	}
    }
  return(index);
}

double
ContourSubElement::parametric_quality() const
{
  const int nelem = 1;
  const auto & masterElement = get_master_element();
  const int numNodes = get_num_nodes();
  const stk::math::Vector3d * nodeCoords = get_coordinates_at_nodes();
  const int nint  = masterElement.num_intg_pts();
  const int dim   = spatial_dim();
  std::vector<double> coords(numNodes * dim, 0.);
  std::vector<double> det_J(nint, 0.);
  double error = 0.;
  double sub_quality = 0.;

  // load coords
  int count = 0;
  for ( int i = 0; i < numNodes; i++ )
    {
      for ( int j = 0; j < dim; j++ )
	{
	  coords[count++] = nodeCoords[i][j];
	}
    }

  // determinant at integration points
  masterElement.determinant( dim, nelem, coords.data(), det_J.data(), &error );

  double min_det_J = 0., sum_det_J = 0.;
  for ( int ip = 0; ip < nint; ip++ )
    {
      if ( ip == 0 || det_J[ip] < min_det_J )
	min_det_J = det_J[ip];
      sum_det_J += std::fabs( det_J[ip] );
    }

  if ( sum_det_J < std::pow(std::numeric_limits<double>::epsilon(),1./dim) )
    sub_quality = 1.; // element too small to consider
  else
    sub_quality = min_det_J * nint / sum_det_J;

  return sub_quality;
}

double
ContourSubElement::physical_quality() const
{
  const int nelem = 1;
  const auto & masterElement = get_master_element();
  const int numNodes = get_num_nodes();
  const stk::math::Vector3d * nodeCoords = get_coordinates_at_nodes();
  const int nint  = masterElement.num_intg_pts();
  const int dim   = spatial_dim();
  std::vector<double> coords(numNodes * dim, 0.);
  std::vector<double> det_J(nint, 0.);
  double error = 0.;
  double sub_quality = 0.;

  // load coords
  int count = 0;
  for ( int i = 0; i < numNodes; i++ )
    {
      const stk::math::Vector3d phys_coords = my_owner->coordinates(nodeCoords[i]);
      for ( int j = 0; j < dim; j++ )
	{
	  coords[count++] = phys_coords[j];
	}
    }

  // determinant at integration points
  masterElement.determinant( dim, nelem, coords.data(), det_J.data(), &error );

  double min_det_J = 0., sum_det_J = 0.;
  for ( int ip = 0; ip < nint; ip++ )
    {
      if ( ip == 0 || det_J[ip] < min_det_J )
	min_det_J = det_J[ip];
      sum_det_J += std::fabs( det_J[ip] );
    }

  if ( sum_det_J < std::pow(std::numeric_limits<double>::epsilon(),1./dim) * std::pow(my_owner->length_scale(),1.*dim) )
    sub_quality = 1.; // element too small to consider
  else
    sub_quality = min_det_J * nint / sum_det_J;

  return sub_quality;
}

double
ContourSubElement::side_quality(const int side) const
{
  const int nelem = 1;
  const auto & sideMasterElement = get_side_master_element();
  const stk::math::Vector3d * nodeCoords = get_coordinates_at_nodes();
  const int nint  = sideMasterElement.num_intg_pts();
  const stk::topology Top = topology();
  const stk::topology sideTop = Top.side_topology(side);
  const int side_num_nodes = sideTop.num_nodes();
  const unsigned * const lnn = get_side_node_ordinals(Top, side);
  const int dim = spatial_dim();
  std::vector<double> coords(side_num_nodes * dim, 0.);
  std::vector<double> det_J(nint, 0.);
  double error = 0.;
  double quality = 0.;

  // load coords on side
  int count = 0;
  for ( int i = 0; i < side_num_nodes; i++ )
    {
      for ( int j = 0; j < dim; j++ )
	{
	  coords[count++] = nodeCoords[lnn[i]][j];
	}
    }

  // determinant at integration points
  sideMasterElement.determinant( dim, nelem, coords.data(), det_J.data(), &error );

  double min_det_J = 0., sum_det_J = 0.;
  for ( int ip = 0; ip < nint; ip++ )
    {
      if ( ip == 0 || det_J[ip] < min_det_J )
	min_det_J = det_J[ip];
      sum_det_J += std::fabs( det_J[ip] );
    }

  if ( sum_det_J < std::pow(std::numeric_limits<double>::epsilon(),1./dim) )
    quality = 1.; // element too small to consider
  else
    quality = min_det_J * nint / sum_det_J;

  return quality;
}

double
ContourSubElement::find_quadratic_crossing( double d0,
				     double d1,
				     double d2 )
{
  const double epsilon = std::numeric_limits<double>::epsilon()*std::sqrt(d0*d0 + d1*d1 + d2*d2);
  if ( std::fabs(d0) < epsilon ) return 0.0;
  if ( std::fabs(d1) < epsilon ) return 1.0;
  if ( std::fabs(d2) < epsilon ) return 0.5;

  STK_ThrowAssert(d0*d1 < 0.0 && (d0*d2 < 0.0 || d1*d2 < 0.0)); // Insist on one and only one crossing

  const double a = 2.0*(d0 - 2.0*d2 + d1);
  const double b = -3.0*d0 - d1 + 4.0*d2;
  const double c = d0;
  const int sign_b = ( b < 0.0 ) ? -1 : 1;
  const double q = -0.5*(b + sign_b*std::sqrt(b*b-4.0*a*c));

  const int sign_a = ( a < 0.0 ) ? -1 : 1;

  if (q*sign_a > 0.0 && q*sign_a < a*sign_a)
    {
      STK_ThrowAssert(!(c*(( q < 0.0 ) ? -1 : 1) > 0.0 && c*(( q < 0.0 ) ? -1 : 1) < q*(( q < 0.0 ) ? -1 : 1))); // Insist on only one crossing
      return (q/a);
    }
  else
    {
      STK_ThrowAssert(c*(( q < 0.0 ) ? -1 : 1) > 0.0 && c*(( q < 0.0 ) ? -1 : 1) < q*(( q < 0.0 ) ? -1 : 1));
      return (c/q);
    }
}

ContourSubElement_Quad_4::ContourSubElement_Quad_4(
  const std::array<stk::math::Vector3d,4> & coords,
  const std::array<int,4> &  side_ids,
  const ContourElement * in_owner )
    : ContourSubElementWithTopology<stk::topology::QUAD_4_2D>( coords,
		  side_ids,
		  in_owner,
		  0, /* in_subelement_depth*/
		  0  /* subelement_sign=0 for now, correct this below if this element is entirely on one side */ )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_linear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    non_conformal_decomposition();
  }
}

int
ContourSubElement_Quad_4::non_conformal_decomposition()
{
  int success = true; // optimism

  // create 4, 3-noded, adaptive triangles
  my_subelements.reserve(4);
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,3> sub_coords;
  std::array<int,3> sub_ids;

  stk::math::Vector3d center = 0.25*(myCoords[0]+myCoords[1]+myCoords[2]+myCoords[3]);

  // triangle #1
  sub_coords[0] = myCoords[0];
  sub_coords[1] = myCoords[1];
  sub_coords[2] = center;
  sub_ids[0] = mySideIds[0];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // triangle #2
  sub_coords[0] = myCoords[1];
  sub_coords[1] = myCoords[2];
  sub_coords[2] = center;
  sub_ids[0] = mySideIds[1];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // triangle #3
  sub_coords[0] = myCoords[2];
  sub_coords[1] = myCoords[3];
  sub_ids[0] = mySideIds[2];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // triangle #4
  sub_coords[0] = myCoords[3];
  sub_coords[1] = myCoords[0];
  sub_coords[2] = center;
  sub_ids[0] = mySideIds[3];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  return success;
}

ContourSubElement_Quad_9::ContourSubElement_Quad_9(
  const std::array<stk::math::Vector3d,9> & coords,
  const std::array<int,4> &  sideIds,
  const ContourElement *in_owner )
    : ContourSubElementWithTopology<stk::topology::QUAD_9_2D>( coords,
		  sideIds,
		  in_owner,
		  0, /* in_subelement_depth*/
		  0  /* subelement_sign=0 for now, correct this below if this element is entirely on one side */ )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_conservative_nonlinear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    non_conformal_decomposition();
  }
}

int
ContourSubElement_Quad_9::non_conformal_decomposition()
{
  int success = true; // optimism

  // create 4, 3-noded, adaptive triangles
  my_subelements.reserve(4);
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,3> sub_coords;
  std::array<int,3> sub_ids;

  // triangle #1
  sub_coords[0] = myCoords[0];
  sub_coords[1] = myCoords[1];
  sub_coords[2] = myCoords[8];
  sub_ids[0] = mySideIds[0];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // triangle #2
  sub_coords[0] = myCoords[1];
  sub_coords[1] = myCoords[2];
  sub_coords[2] = myCoords[8];
  sub_ids[0] = mySideIds[1];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // triangle #3
  sub_coords[0] = myCoords[2];
  sub_coords[1] = myCoords[3];
  sub_coords[2] = myCoords[8];
  sub_ids[0] = mySideIds[2];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // triangle #4
  sub_coords[0] = myCoords[3];
  sub_coords[1] = myCoords[0];
  sub_coords[2] = myCoords[8];
  sub_ids[0] = mySideIds[3];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  return success;
}

ContourSubElement_Tri_3::ContourSubElement_Tri_3(
  const std::array<stk::math::Vector3d,3> & coords,
  const std::array<int,3> &  sideIds,
  const ContourElement * in_owner,
  const int in_subelement_depth,
  const int subelement_sign )
    : ContourSubElementWithTopology<stk::topology::TRI_3_2D>( coords,
		  sideIds,
		  in_owner,
		  in_subelement_depth,
		  subelement_sign )
{
  // if this is a conformal element, return quickly
  if ( subelement_sign != 0 )
    {
      return;
    }

  // snap to mesh
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  for (int n = 0; n < 3; ++n)
    {
      if (std::fabs(myDist[n]) < snapTol )
	{
	  myDist[n] = 0.0;
	}
    }


  my_sign = ContourElement::compute_linear_distance_sign(snapTol, myDist.size(), myDist.data());

  if ( my_sign == 0 )
    {
      // attempt conformal decomposition
      int success = conformal_decomposition();
      STK_ThrowErrorMsgIf(!success, " Conformal decomposition failed.\n");
    }
}

int
ContourSubElement_Tri_3::conformal_decomposition()
{

  // create 4 conforming triangular subelements

  // create 4, 3-noded tris
  my_subelements.clear();
  my_subelements.reserve(4);
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,3> sub_coords;
  std::array<int,3> sub_ids;

  // For any edge with a crossing, we will move the
  // mid side node for that egdge to the crossing
  // we will keep the modified locations of the nodes
  // in a local vector of nodes (lcoords).
  // We will also create local vectors for the distance and side_ids
  // so that we can reorient the tri as discussed below.
  std::array<stk::math::Vector3d,6> lcoords = {myCoords[0], myCoords[1], myCoords[2], stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO};
  std::array<double,3> ldist = { myDist[0], myDist[1], myDist[2] };
  std::array<int,6> is_on_surf = {0,0,0,0,0,0};
  int sub_sign;
  std::array<int,6> edge_node_ids;

  // find edge crossings
  edge_node_ids[0] = 0;
  edge_node_ids[1] = 1;
  edge_node_ids[2] = 2;
  edge_node_ids[3] = process_edge( 0, 1, 3, is_on_surf, lcoords, ldist );
  edge_node_ids[4] = process_edge( 1, 2, 4, is_on_surf, lcoords, ldist );
  edge_node_ids[5] = process_edge( 2, 0, 5, is_on_surf, lcoords, ldist );

  const int zero_sign = LevelSet::sign(0.0);
  std::array<bool,4> sub_degenerate;

  sub_degenerate[0] = is_degenerate(edge_node_ids,0,3,5);
  sub_degenerate[1] = is_degenerate(edge_node_ids,3,1,4);
  sub_degenerate[2] = is_degenerate(edge_node_ids,5,4,2);
  sub_degenerate[3] = is_degenerate(edge_node_ids,3,4,5);

  // tri #1
  if (!sub_degenerate[0])
    {
      sub_coords[0] = lcoords[0];
      sub_coords[1] = lcoords[3];
      sub_coords[2] = lcoords[5];
      sub_sign = LevelSet::sign(myDist[0]);
      sub_ids[0] = mySideIds[0];
      sub_ids[1] = ((is_on_surf[3] && is_on_surf[5]) && (zero_sign != sub_sign || sub_degenerate[3])) ? -2 : -1;
      sub_ids[2] = mySideIds[2];
      sub = new ContourSubElement_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tri #2
  if (!sub_degenerate[1])
    {
      sub_coords[0] = lcoords[3];
      sub_coords[1] = lcoords[1];
      sub_coords[2] = lcoords[4];
      sub_sign = LevelSet::sign(myDist[1]);
      sub_ids[0] = mySideIds[0];
      sub_ids[1] = mySideIds[1];
      sub_ids[2] = ((is_on_surf[3] && is_on_surf[4]) && (zero_sign != sub_sign || sub_degenerate[3])) ? -2 : -1;
      sub = new ContourSubElement_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tri #3
  if (!sub_degenerate[2])
    {
      sub_coords[0] = lcoords[5];
      sub_coords[1] = lcoords[4];
      sub_coords[2] = lcoords[2];
      sub_sign = LevelSet::sign(myDist[2]);
      sub_ids[0] = ((is_on_surf[5] && is_on_surf[4]) && (zero_sign != sub_sign || sub_degenerate[3])) ? -2 : -1;
      sub_ids[1] = mySideIds[1];
      sub_ids[2] = mySideIds[2];
      sub = new ContourSubElement_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tri #4
  if (!sub_degenerate[3])
    {
      sub_coords[0] = lcoords[3];
      sub_coords[1] = lcoords[4];
      sub_coords[2] = lcoords[5];
      sub_sign = LevelSet::sign( (is_on_surf[3] ? 0.0 : myDist[0]+myDist[1]) +
				 (is_on_surf[4] ? 0.0 : myDist[1]+myDist[2]) +
				 (is_on_surf[5] ? 0.0 : myDist[2]+myDist[0]) );
      sub_ids[0] = ((is_on_surf[3] && is_on_surf[4]) && (zero_sign != sub_sign || sub_degenerate[1])) ? -2 : -1;
      sub_ids[1] = ((is_on_surf[4] && is_on_surf[5]) && (zero_sign != sub_sign || sub_degenerate[2])) ? -2 : -1;
      sub_ids[2] = ((is_on_surf[5] && is_on_surf[3]) && (zero_sign != sub_sign || sub_degenerate[0])) ? -2 : -1;
      sub = new ContourSubElement_Tri_3( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // check quality of subelements
  // Here we assume that the linear tri is always decomposed into reasonable quality sub-tris
  int success = true;

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    for ( unsigned i=0; i<my_subelements.size(); ++i )
    {
      if ( my_subelements[i]->physical_quality() < 0.0 )
      {
        krinolog << "low quality subelement: " << i << "\n"
                << *my_subelements[i] << "\n"
                << "parent:" << "\n"
                << *this << "\n";
      }
    }
  }

  return success;
}

int
ContourSubElement_Tri_3::process_edge( const int i0,
				const int i1,
				const int i2,
				std::array<int,6> & is_on_surf,
				std::array<stk::math::Vector3d,6> & lcoords,
				const std::array<double,3> & ldist )
{
  int edge_node_id = i2;
  is_on_surf[i2] = LevelSet::sign_change( ldist[i0], ldist[i1] );
  if ( is_on_surf[i2] )
    {
      // tolerance chosen very small since degeneracies should already be eliminated
      const double tol = std::numeric_limits<double>::epsilon();

      const double d0 = std::fabs(ldist[i0]);
      const double d1 = std::fabs(ldist[i1]);

      // make calculation completely symmetric
      if (d0 < d1)
	{
	  const double alpha = my_owner->dist_is_linear() ? d0/(d0+d1) : find_quadratic_crossing(ldist[i0],ldist[i1],my_owner->distance(0.5*(lcoords[i0]+lcoords[i1])));
	  if (alpha < tol)
	    {
	      edge_node_id = i0;
	      lcoords[i2] = lcoords[edge_node_id];
	    }
	  else
	    {
	      lcoords[i2] = (1.-alpha) * lcoords[i0] + alpha * lcoords[i1];
	    }
	}
      else
	{
	  const double alpha = my_owner->dist_is_linear() ? d1/(d1+d0) : find_quadratic_crossing(ldist[i1],ldist[i0],my_owner->distance(0.5*(lcoords[i1]+lcoords[i0])));
	  if (alpha < tol)
	    {
	      edge_node_id = i1;
	      lcoords[i2] = lcoords[edge_node_id];
	    }
	  else
	    {
	      lcoords[i2] = (1.-alpha) * lcoords[i1] + alpha * lcoords[i0];
	    }
	}
    }
  else
    {
      // eliminate side node by sliding to one end or the other
      const double d0 = std::fabs(ldist[i0]);
      const double d1 = std::fabs(ldist[i1]);
      const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon()) * (d0 + d1);

      if ( d0 > d1 + epsilon )
	{
	  edge_node_id = i0;
	}
      else if ( d1 > d0 + epsilon )
	{
	  edge_node_id = i1;
	}
      else
	{
	  // tie breaker
	  const stk::math::Vector3d phys0 = my_owner->coordinates(lcoords[i0]);
	  const stk::math::Vector3d phys1 = my_owner->coordinates(lcoords[i1]);

	  if ( is_more(phys1,phys0) )
	    {
	      edge_node_id = i0;
	    }
	  else
	    {
	      edge_node_id = i1;
	    }
	}
      lcoords[i2] = lcoords[edge_node_id];
    }
  return edge_node_id;
}

bool
ContourSubElement_Tri_3::is_degenerate( const std::array<int,6> & edge_node_ids,
				 const int i0, const int i1, const int i2 )
{

  // DRN:  This is really ugly, hand-optimized code for looking for degenerate tris
  //       Basically, it checks if any of the edges are degenerate. Then it has to look for
  //       the entire tri being degerate because it consists of 3 colinear points.
  //       This is handled by checking against the 3 specific bad cases.

  if ( edge_node_ids[i0] == edge_node_ids[i1] ||
       edge_node_ids[i0] == edge_node_ids[i2] ||
       edge_node_ids[i1] == edge_node_ids[i2] )
    {
      // this tri is degenerate with two coincident nodes
      return true;
    }

  if ( edge_node_ids[i0]==i0 &&
       edge_node_ids[i1]==i1 &&
       edge_node_ids[i2]==i2 )
    {
      // this tri is not degenerate since is has no degenerate nodes
      return false;
    }

  // look for a colinear triangle
  std::vector<int> is_used(6); // initializes to zero (false);
  is_used[edge_node_ids[i0]] = true;
  is_used[edge_node_ids[i1]] = true;
  is_used[edge_node_ids[i2]] = true;

  if ((is_used[0] && ((is_used[1] && is_used[3]) || (is_used[2] && is_used[5]))) ||
      (is_used[1] && is_used[2] && is_used[4]))
    {
      // this tri is colinear
      return true;
    }

  return false;
}

int
ContourSubElement_Tri_3::side_facets( FacetedSurfaceBase & facets,
			       int side ) const
{
  STK_ThrowAssert( get_side_ids()[side] == -2 );

  // just one linear facet per side
  const int num_facets = 1;

  const unsigned * const lnn = get_side_node_ordinals(topology(), side);

  if ( LevelSet::sign_change(0.0, (double) my_sign) )
    facets.emplace_back_2d( my_owner->coordinates(myCoords[lnn[0]]), my_owner->coordinates(myCoords[lnn[1]]) );
  else
    facets.emplace_back_2d( my_owner->coordinates(myCoords[lnn[1]]), my_owner->coordinates(myCoords[lnn[0]]) );

  return( num_facets );
}

double
ContourSubElement_Tri_3::side_area( int side ) const
{
  STK_ThrowAssert( get_side_ids()[side] == -2 );

  const unsigned * const lnn = get_side_node_ordinals(topology(), side);

  return (my_owner->coordinates(myCoords[lnn[0]]) - my_owner->coordinates(myCoords[lnn[1]])).length();
}

const int ContourSubElement_Adaptive_Tri_3::MAX_REFINMENT_LEVELS = 6;

ContourSubElement_Adaptive_Tri_3::ContourSubElement_Adaptive_Tri_3(
  const std::array<stk::math::Vector3d,3> & coords,
  const std::array<int,3> &  sideIds,
  const ContourElement * in_owner,
  const int in_subelement_depth )
    : ContourSubElementWithTopology<stk::topology::TRI_3_2D>( coords,
                  sideIds,
                  in_owner,
                  in_subelement_depth,
                  0 )
{
  my_edge_age.fill(0);
  non_conformal_decomposition();
}

ContourSubElement_Adaptive_Tri_3::ContourSubElement_Adaptive_Tri_3(
  const std::array<stk::math::Vector3d,3> & coords,
  const std::array<int,3> & side_ids,
  const std::array<int,3> &  edge_age,
  const ContourElement * in_owner,
  const int in_subelement_depth )
    : ContourSubElementWithTopology<stk::topology::TRI_3_2D>( coords,
                  side_ids,
                  in_owner,
                  in_subelement_depth,
                  0 )
{
  my_edge_age = edge_age;
  non_conformal_decomposition();
}

int
ContourSubElement_Adaptive_Tri_3::non_conformal_decomposition()
{
  int success = true; // optimism

  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_conservative_nonlinear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign != 0)
  {
    return success;
  }

  int longest_bad_edge = -1;

  std::array<stk::math::Vector3d,6> lcoords = {myCoords[0], myCoords[1], myCoords[2], 0.5*(myCoords[0] + myCoords[1]), 0.5*(myCoords[1] + myCoords[2]), 0.5*(myCoords[2] + myCoords[0])};

  std::array<stk::math::Vector3d,6> lphyscoords;
  for (int n = 0; n < 6; ++n)
    lphyscoords[n] = my_owner->coordinates( lcoords[n] );

  std::array<double,6> ldist = {myDist[0], myDist[1], myDist[2], 0., 0., 0.};
  for (int n = 3; n < 6; ++n)
    ldist[n] = my_owner->distance( lcoords[n] );

  const stk::topology Top = stk::topology::TRIANGLE_6_2D;
  int num_edges = Top.num_edges();

  std::vector<int> bad_edges;
  bad_edges.reserve(num_edges);

  std::vector<double> edge_lengths(num_edges);

  for ( int edge = 0; edge < num_edges; edge++ )
    {
      const unsigned * const lnn = get_edge_node_ordinals(Top, edge);

      STK_ThrowAssert(Top.edge_topology(edge).num_nodes() == 3);

      const double edge_straight_length = (lphyscoords[lnn[0]] - lphyscoords[lnn[1]]).length();
      STK_ThrowRequire(edge_straight_length > 0.0);
      edge_lengths[edge] = edge_straight_length;

      const double edge_curve_error = (lphyscoords[lnn[2]] - 0.5*(lphyscoords[lnn[0]] + lphyscoords[lnn[1]])).length();

      const double edge_dist_error = std::fabs(ldist[lnn[2]] - 0.5*(ldist[lnn[0]]+ldist[lnn[1]]));

      const double scale = std::min(std::sqrt(std::numeric_limits<double>::max()),std::fabs(ldist[lnn[0]]) + std::fabs(ldist[lnn[1]]) + my_owner->length_scale());

      const double edge_error = (edge_curve_error + edge_dist_error)*edge_straight_length/(scale*scale);

      if (edge_error > my_owner->edge_nonlinear_tolerance() && my_edge_age[edge] < MAX_REFINMENT_LEVELS)
	{
	  bad_edges.push_back(edge);
	}
    }

  double max_length = 0.0;
  for (auto edge : bad_edges)
  {
    const double edge_length = edge_lengths[edge];
    STK_ThrowRequire(edge_length > 0.0);

    // we need an absolute mechanism for selecting the edge to bisect so that all elements that share
    // common edges will make the same decisions
    if (utility::is_more(edge_length,max_length))
    {
      longest_bad_edge = edge;
      max_length = edge_length;
    }
    else if (!utility::is_less(edge_length,max_length)) // tie breaker
    {
      const stk::math::Vector3d & edge_midside_coords = lphyscoords[get_edge_node_ordinals(Top, edge)[2]];
      // note that it is safe to assume that longest_bad_edge is already assigned if edge_length == max_length
      const stk::math::Vector3d longest_edge_midside_coords = lphyscoords[get_edge_node_ordinals(Top, longest_bad_edge)[2]];

      STK_ThrowAssert((utility::is_not_equal(edge_midside_coords[0],longest_edge_midside_coords[0]) ||
                   utility::is_not_equal(edge_midside_coords[1],longest_edge_midside_coords[1])));

      if (utility::is_more(edge_midside_coords[0],longest_edge_midside_coords[0]) ||
          (!utility::is_less(edge_midside_coords[0],longest_edge_midside_coords[0]) &&
           (utility::is_more(edge_midside_coords[1],longest_edge_midside_coords[1]))))
        {
          longest_bad_edge = edge;
          max_length = edge_length;
        }
      }
    }

  if ( longest_bad_edge == -1 )
    {
      // no bad edges

      // use a single nonconformal linear tet subelement
      my_subelements.clear();
      my_subelements.reserve(1);

      ContourSubElement *sub = new ContourSubElement_Tri_3( myCoords, mySideIds, my_owner, my_subelement_depth+1 );
      my_subelements.push_back( sub );
    }
  else
    {
      //
      // create 2, adaptive, 3-noded triangles by cutting the longest_bad_edge
      //

      my_subelements.clear();
      my_subelements.reserve(2);
      ContourSubElement *sub  = NULL;
      std::array<stk::math::Vector3d,3> sub_coords;
      std::array<int,3> sub_ids;
      std::array<int,3> sub_edge_age;

      static const unsigned permute_0[] = { 0,1,2 };
      static const unsigned permute_1[] = { 1,2,0 };
      static const unsigned permute_2[] = { 2,0,1 };
      static const unsigned * permute_table[] = { permute_0, permute_1, permute_2 };

      const unsigned * lnn = permute_table[longest_bad_edge];
      const unsigned * lsn = lnn; // side permutation mirrors node permutation

      const stk::math::Vector3d edge_node = 0.5 * (myCoords[lnn[0]] + myCoords[lnn[1]]);

      // tri #1
      sub_coords[0] = myCoords[lnn[0]];
      sub_coords[1] = edge_node;
      sub_coords[2] = myCoords[lnn[2]];
      sub_ids[0] = mySideIds[lsn[0]];
      sub_ids[1] = -1; /* not on any parent side */
      sub_ids[2] = mySideIds[lsn[2]];
      sub_edge_age[0] = my_edge_age[lsn[0]]+1;
      sub_edge_age[1] = my_edge_age[lsn[0]]+1;
      sub_edge_age[2] = my_edge_age[lsn[2]];
      sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, sub_edge_age, my_owner, my_subelement_depth+1 );
      my_subelements.push_back( sub );

      // tri #2
      sub_coords[0] = edge_node;
      sub_coords[1] = myCoords[lnn[1]];
      sub_coords[2] = myCoords[lnn[2]];
      sub_ids[0] = mySideIds[lsn[0]];
      sub_ids[1] = mySideIds[lsn[1]];
      sub_edge_age[0] = my_edge_age[lsn[0]]+1;
      sub_edge_age[1] = my_edge_age[lsn[1]];
      sub_edge_age[2] = my_edge_age[lsn[0]]+1;
      sub_ids[2] = -1; /* not on any parent side */
      sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, sub_ids, sub_edge_age, my_owner, my_subelement_depth+1 );
      my_subelements.push_back( sub );
    }

  return success;
}

ContourSubElement_Tri_6::ContourSubElement_Tri_6(
  const std::array<stk::math::Vector3d,6> & coords,
  const std::array<int,3> &  sideIds,
  const ContourElement * in_owner,
  const int in_subelement_depth,
  const int subelement_sign )
    : ContourSubElementWithTopology<stk::topology::TRI_6_2D> ( coords,
		  sideIds,
		  in_owner,
		  in_subelement_depth,
		  subelement_sign )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_conservative_nonlinear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    // use a single non-conformal, adaptive 3-noded tri
    const std::array<stk::math::Vector3d,3> sub_coords = { myCoords[0], myCoords[1], myCoords[2] };
    ContourSubElement *sub = new ContourSubElement_Adaptive_Tri_3( sub_coords, mySideIds, my_owner, my_subelement_depth+1 );
    my_subelements.push_back( sub );
  }
}

ContourSubElement_Hex_8::ContourSubElement_Hex_8(
  const std::array<stk::math::Vector3d,8> & coords,
  const std::array<int,6> &  sideIds,
  const ContourElement * in_owner )
    : ContourSubElementWithTopology<stk::topology::HEX_8> ( coords,
		  sideIds,
		  in_owner,
		  0, /* in_subelement_depth*/
		  0  /* subelement_sign=0 for now, correct this below if this element is entirely on one side */ )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_linear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    // create 24, 4-noded, adaptive tetrahedra
    my_subelements.reserve(24);

    // Conceptually, hex is broken into 6 prisms, with the
    // bases of the prisms corresponding to a face of the hex.
    int success = true; // optimism
    for ( int face = 0; face < 6 && success; ++face )
    {
      success &= subpyramid_non_conformal_decomposition( face );
    }
  }
}

int
ContourSubElement_Hex_8::subpyramid_non_conformal_decomposition( const int face )
{
  int success = true; // optimism
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,4> sub_coords;
  std::array<int,4> sub_ids;

  static const unsigned face_0[] = { 0,1,5,4 };
  static const unsigned face_1[] = { 1,2,6,5 };
  static const unsigned face_2[] = { 2,3,7,6 };
  static const unsigned face_3[] = { 0,4,7,3 };
  static const unsigned face_4[] = { 0,3,2,1 };
  static const unsigned face_5[] = { 4,5,6,7 };

  static const unsigned * face_table[] = { face_0 , face_1 , face_2 , face_3 , face_4 , face_5 };

  const unsigned * lnn = face_table[face];

  //
  // create 4, 4-noded adaptive tetrahedra
  //
  // The advantage of 4 tets per face over 2 tets per face is that all corners
  // will be bisected.  This eliminates some of the pathologies that occur when
  // 3 nodes have the same value while the 4th node on the face has a different
  // sign. Note that this problem can be mitigated, however, if the non-conformal
  // refinement of the sub-tets will do longest edge bisection rather than the
  // self-similar 8 subtet refinement.

  stk::math::Vector3d vol_center = 0.125*(myCoords[0]+myCoords[1]+myCoords[2]+myCoords[3]+
			      myCoords[4]+myCoords[5]+myCoords[6]+myCoords[7]);
  stk::math::Vector3d face_center = 0.25*(myCoords[lnn[0]]+myCoords[lnn[1]]+myCoords[lnn[2]]+myCoords[lnn[3]]);

  // tet #1
  sub_coords[0] = myCoords[lnn[0]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[1]];
  sub_coords[3] = vol_center;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #2
  sub_coords[0] = myCoords[lnn[1]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[2]];
  sub_coords[3] = vol_center;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #3
  sub_coords[0] = myCoords[lnn[2]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[3]];
  sub_coords[3] = vol_center;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #4
  sub_coords[0] = myCoords[lnn[3]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[0]];
  sub_coords[3] = vol_center;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  return success;
}

ContourSubElement_Hex_27::ContourSubElement_Hex_27(
  const std::array<stk::math::Vector3d,27> & coords,
  const std::array<int,6> &  sideIds,
  const ContourElement * in_owner )
    : ContourSubElementWithTopology<stk::topology::HEX_27> ( coords,
		  sideIds,
		  in_owner,
		  0, /* in_subelement_depth*/
		  0  /* subelement_sign=0 for now, correct this below if this element is entirely on one side */ )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_conservative_nonlinear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    // create 24, 4-noded, adaptive tetrahedra
    my_subelements.reserve(24);

    // Conceptually, hex is broken into 6 prisms, with the
    // bases of the prisms corresponding to a face of the hex.
    int success = true; // optimism
    for ( int face = 0; face < 6 && success; ++face )
    {
      success &= subpyramid_non_conformal_decomposition( face );
    }
  }
}

int
ContourSubElement_Hex_27::subpyramid_non_conformal_decomposition( const int face )
{
  int success = true; // optimism
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,4> sub_coords;
  std::array<int,4> sub_ids;

  static const unsigned face_0[] = { 0,1,5,4, 25 };
  static const unsigned face_1[] = { 1,2,6,5, 24 };
  static const unsigned face_2[] = { 2,3,7,6, 26 };
  static const unsigned face_3[] = { 0,4,7,3, 23 };
  static const unsigned face_4[] = { 0,3,2,1, 21 };
  static const unsigned face_5[] = { 4,5,6,7, 22 };

  static const unsigned * face_table[] = { face_0 , face_1 , face_2 , face_3 , face_4 , face_5 };

  const unsigned * lnn = face_table[face];

  //
  // create 4, 4-noded adaptive tetrahedra
  //
  // The advantage of 4 tets per face over 2 tets per face is that all corners
  // will be bisected.  This eliminates some of the pathologies that occur when
  // 3 nodes have the same value while the 4th node on the face has a different
  // sign. Note that this problem can be mitigated, however, if the non-conformal
  // refinement of the sub-tets will do longest edge bisection rather than the
  // self-similar 8 subtet refinement.

  // tet #1
  sub_coords[0] = myCoords[lnn[0]];
  sub_coords[1] = myCoords[lnn[4]];
  sub_coords[2] = myCoords[lnn[1]];
  sub_coords[3] = myCoords[20];
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #2
  sub_coords[0] = myCoords[lnn[1]];
  sub_coords[1] = myCoords[lnn[4]];
  sub_coords[2] = myCoords[lnn[2]];
  sub_coords[3] = myCoords[20];
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #3
  sub_coords[0] = myCoords[lnn[2]];
  sub_coords[1] = myCoords[lnn[4]];
  sub_coords[2] = myCoords[lnn[3]];
  sub_coords[3] = myCoords[20];
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #4
  sub_coords[0] = myCoords[lnn[3]];
  sub_coords[1] = myCoords[lnn[4]];
  sub_coords[2] = myCoords[lnn[0]];
  sub_coords[3] = myCoords[20];
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  return success;
}

ContourSubElement_Wedge_6::ContourSubElement_Wedge_6(
  const std::array<stk::math::Vector3d,6> & coords,
  const std::array<int,5> &  sideIds,
  const ContourElement * in_owner )
    : ContourSubElementWithTopology<stk::topology::WEDGE_6>( coords,
                  sideIds,
                  in_owner,
                  0, /* in_subelement_depth*/
                  0  /* subelement_sign=0 for now, correct this below if this element is entirely on one side */ )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_linear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    // create 12, 4-noded, adaptive tetrahedra
    my_subelements.reserve(12);

    // Conceptually, hex is broken into 6 prisms, with the
    // bases of the prisms corresponding to a face of the hex.
    int success = true; // optimism
    for ( int face = 0; face < 3 && success; ++face )
    {
      success &= subpyramid_non_conformal_decomposition( face );
    }
  }
}

int
ContourSubElement_Wedge_6::subpyramid_non_conformal_decomposition( const int face )
{
  int success = true; // optimism
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,4> sub_coords;
  std::array<int,4> sub_ids;

  static const unsigned face_0[] = {0, 1, 4, 3};
  static const unsigned face_1[] = {1, 2, 5, 4};
  static const unsigned face_2[] = {0, 3, 5, 2};
  static const unsigned * face_table[] = { face_0 , face_1 , face_2 };

  const unsigned * lnn = face_table[face];

  //
  // create 4, 4-noded adaptive tetrahedra
  //
  // The advantage of 4 tets per face over 2 tets per face is that all corners
  // will be bisected.  This eliminates some of the pathologies that occur when
  // 3 nodes have the same value while the 4th node on the face has a different
  // sign. Note that this problem can be mitigated, however, if the non-conformal
  // refinement of the sub-tets will do longest edge bisection rather than the
  // self-similar 8 subtet refinement.

  // Not guaranteed to be within a highly deformed wedge
  const stk::math::Vector3d centroid = (myCoords[0]+myCoords[1]+myCoords[2]+myCoords[3]+myCoords[4]+myCoords[5])/6.;
  const stk::math::Vector3d face_center = 0.25*(myCoords[lnn[0]]+myCoords[lnn[1]]+myCoords[lnn[2]]+myCoords[lnn[3]]);

  // tet #1
  sub_coords[0] = myCoords[lnn[0]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[1]];
  sub_coords[3] = centroid;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #2
  sub_coords[0] = myCoords[lnn[1]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[2]];
  sub_coords[3] = centroid;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #3
  sub_coords[0] = myCoords[lnn[2]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[3]];
  sub_coords[3] = centroid;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  // tet #4
  sub_coords[0] = myCoords[lnn[3]];
  sub_coords[1] = face_center;
  sub_coords[2] = myCoords[lnn[0]];
  sub_coords[3] = centroid;
  sub_ids[0] = mySideIds[face];
  sub_ids[1] = -1; /* not on any parent side */
  sub_ids[2] = -1; /* not on any parent side */
  sub_ids[3] = -1; /* not on any parent side */
  sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1 );
  my_subelements.push_back( sub );

  return success;
}

ContourSubElement_Tet_4::ContourSubElement_Tet_4(
  const std::array<stk::math::Vector3d,4> & coords,
  const std::array<int,4> &  sideIds,
  const ContourElement * in_owner,
  const int in_subelement_depth,
  const int subelement_sign )
    : ContourSubElementWithTopology<stk::topology::TET_4>( coords,
		  sideIds,
		  in_owner,
		  in_subelement_depth,
		  subelement_sign )
{
  // if this is a conformal element, return quickly
  if ( subelement_sign != 0 )
    {
      return;
    }

  // snap to mesh
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  for (int n = 0; n < 4; ++n)
    {
      if (std::fabs(myDist[n]) < snapTol)
	{
	  myDist[n] = 0.0;
	}
    }

  my_sign = ContourElement::compute_linear_distance_sign(snapTol, myDist.size(), myDist.data());

  if ( my_sign == 0 )
    {
      // attempt conformal decomposition
      int success = conformal_decomposition();
      STK_ThrowErrorMsgIf(!success, " Conformal decomposition failed.\n");
    }
}

int
ContourSubElement_Tet_4::conformal_decomposition()
{

  // attempt to create 8 conforming tetrahedral subelements
  // This attempt may unsuccessful if the resulting subelements
  // are of poor quality
  int success = true; // optimism

  // create 8, 4-noded tets
  my_subelements.clear();
  my_subelements.reserve(8);
  ContourSubElement *sub  = NULL;
  std::array<stk::math::Vector3d,4> sub_coords;
  std::array<int,4> sub_ids;

  // For any edge with a crossing, we will move the
  // mid side node for that egdge to the crossing
  // we will keep the modified locations of the nodes
  // in a local vector of nodes (lcoords).
  // We will also create local vectors for the distance and side_ids
  // so that we can reorient the tet as discussed below.
  std::array<stk::math::Vector3d,10> lcoords = {myCoords[0], myCoords[1], myCoords[2], myCoords[3],
      stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO, stk::math::Vector3d::ZERO };
  std::array<int,4> lsides = mySideIds;
  std::array<double,4> ldist = { myDist[0], myDist[1], myDist[2], myDist[3] };
  std::array<int,10> is_on_surf = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  int sub_sign;

  // Find orientation of tet
  // Specifically, orient such that we don't have nodes 0 and 2 on one side and
  // nodes 1 and 3 on the other
  if ( LevelSet::sign_change(myDist[0],myDist[1]) &&
       LevelSet::sign_change(myDist[1],myDist[2]) &&
       LevelSet::sign_change(myDist[2],myDist[3]) )
    {
      lcoords[0] = myCoords[0];
      lcoords[1] = myCoords[3];
      lcoords[2] = myCoords[1];
      lcoords[3] = myCoords[2];
      ldist[0] = myDist[0];
      ldist[1] = myDist[3];
      ldist[2] = myDist[1];
      ldist[3] = myDist[2];
      lsides[0] = mySideIds[2];
      lsides[1] = mySideIds[1];
      lsides[2] = mySideIds[3];
      lsides[3] = mySideIds[0];
    }

  std::array<int,10> edge_node_ids;
  edge_node_ids[0] = 0;
  edge_node_ids[1] = 1;
  edge_node_ids[2] = 2;
  edge_node_ids[3] = 3;
  edge_node_ids[4] = process_edge( 0, 1, 4, is_on_surf, lcoords, ldist );
  edge_node_ids[5] = process_edge( 1, 2, 5, is_on_surf, lcoords, ldist );
  edge_node_ids[6] = process_edge( 2, 0, 6, is_on_surf, lcoords, ldist );
  edge_node_ids[7] = process_edge( 0, 3, 7, is_on_surf, lcoords, ldist );
  edge_node_ids[8] = process_edge( 1, 3, 8, is_on_surf, lcoords, ldist );
  edge_node_ids[9] = process_edge( 2, 3, 9, is_on_surf, lcoords, ldist );

  const int zero_sign = LevelSet::sign(0.0);
  std::vector<int> sub_degenerate(8); // initializes to zero (false)

  sub_degenerate[0] = is_degenerate(edge_node_ids,0,4,6,7);
  sub_degenerate[1] = is_degenerate(edge_node_ids,4,1,5,8);
  sub_degenerate[2] = is_degenerate(edge_node_ids,6,5,2,9);
  sub_degenerate[3] = is_degenerate(edge_node_ids,7,8,9,3);
  sub_degenerate[4] = is_degenerate(edge_node_ids,8,7,6,4);
  sub_degenerate[5] = is_degenerate(edge_node_ids,6,9,8,5);
  sub_degenerate[6] = is_degenerate(edge_node_ids,9,8,7,6);
  sub_degenerate[7] = is_degenerate(edge_node_ids,5,6,4,8);

  // tet #1
  if (!sub_degenerate[0])
    {
      sub_coords[0] = lcoords[0];
      sub_coords[1] = lcoords[4];
      sub_coords[2] = lcoords[6];
      sub_coords[3] = lcoords[7];
      sub_sign = LevelSet::sign(ldist[0]);
      sub_ids[0] = lsides[0];
      sub_ids[1] = ((is_on_surf[4] && is_on_surf[6] && is_on_surf[7]) && (zero_sign != sub_sign || sub_degenerate[4])) ? -2 : -1;
      sub_ids[2] = lsides[2];
      sub_ids[3] = lsides[3];
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #2
  if (!sub_degenerate[1])
    {
      sub_coords[0] = lcoords[4];
      sub_coords[1] = lcoords[1];
      sub_coords[2] = lcoords[5];
      sub_coords[3] = lcoords[8];
      sub_sign = LevelSet::sign(ldist[1]);
      sub_ids[0] = lsides[0];
      sub_ids[1] = lsides[1];
      sub_ids[2] = ((is_on_surf[4] && is_on_surf[5] && is_on_surf[8]) && (zero_sign != sub_sign || sub_degenerate[7])) ? -2 : -1;
      sub_ids[3] = lsides[3];
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #3
  if (!sub_degenerate[2])
    {
      sub_coords[0] = lcoords[6];
      sub_coords[1] = lcoords[5];
      sub_coords[2] = lcoords[2];
      sub_coords[3] = lcoords[9];
      sub_sign = LevelSet::sign(ldist[2]);
      sub_ids[0] = ((is_on_surf[5] && is_on_surf[6] && is_on_surf[9]) && (zero_sign != sub_sign || sub_degenerate[5])) ? -2 : -1;
      sub_ids[1] = lsides[1];
      sub_ids[2] = lsides[2];
      sub_ids[3] = lsides[3];
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #4
  if (!sub_degenerate[3])
    {
      sub_coords[0] = lcoords[7];
      sub_coords[1] = lcoords[8];
      sub_coords[2] = lcoords[9];
      sub_coords[3] = lcoords[3];
      sub_sign = LevelSet::sign(ldist[3]);
      sub_ids[0] = lsides[0];
      sub_ids[1] = lsides[1];
      sub_ids[2] = lsides[2];
      sub_ids[3] = ((is_on_surf[7] && is_on_surf[8] && is_on_surf[9]) && (zero_sign != sub_sign || sub_degenerate[6])) ? -2 : -1;
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #5
  if (!sub_degenerate[4])
    {
      sub_coords[0] = lcoords[8];
      sub_coords[1] = lcoords[7];
      sub_coords[2] = lcoords[6];
      sub_coords[3] = lcoords[4];
      sub_sign = LevelSet::sign( (is_on_surf[8] ? 0.0 : ldist[1]+ldist[3]) +
				 (is_on_surf[7] ? 0.0 : ldist[0]+ldist[3]) +
				 (is_on_surf[6] ? 0.0 : ldist[0]+ldist[2]) +
				 (is_on_surf[4] ? 0.0 : ldist[0]+ldist[1]) );
      sub_ids[0] = lsides[0]; // 8-7-4
      sub_ids[1] = ((is_on_surf[7] && is_on_surf[6] && is_on_surf[4]) && (zero_sign != sub_sign || sub_degenerate[0])) ? -2 : -1;
      sub_ids[2] = ((is_on_surf[8] && is_on_surf[6] && is_on_surf[4]) && (zero_sign != sub_sign || sub_degenerate[7])) ? -2 : -1;
      sub_ids[3] = ((is_on_surf[8] && is_on_surf[7] && is_on_surf[6]) && (zero_sign != sub_sign || sub_degenerate[6])) ? -2 : -1;
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #6
  if (!sub_degenerate[5])
    {
      sub_coords[0] = lcoords[6];
      sub_coords[1] = lcoords[9];
      sub_coords[2] = lcoords[8];
      sub_coords[3] = lcoords[5];
      sub_sign = LevelSet::sign( (is_on_surf[6] ? 0.0 : ldist[0]+ldist[2]) +
				 (is_on_surf[9] ? 0.0 : ldist[2]+ldist[3]) +
				 (is_on_surf[8] ? 0.0 : ldist[1]+ldist[3]) +
				 (is_on_surf[5] ? 0.0 : ldist[1]+ldist[2]) );
      sub_ids[0] = ((is_on_surf[6] && is_on_surf[9] && is_on_surf[5]) && (zero_sign != sub_sign || sub_degenerate[2])) ? -2 : -1;
      sub_ids[1] = lsides[1]; // 8-9-5
      sub_ids[2] = ((is_on_surf[6] && is_on_surf[8] && is_on_surf[5]) && (zero_sign != sub_sign || sub_degenerate[7])) ? -2 : -1;
      sub_ids[3] = ((is_on_surf[6] && is_on_surf[9] && is_on_surf[8]) && (zero_sign != sub_sign || sub_degenerate[6])) ? -2 : -1;
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #7
  if (!sub_degenerate[6])
    {
      sub_coords[0] = lcoords[9];
      sub_coords[1] = lcoords[8];
      sub_coords[2] = lcoords[7];
      sub_coords[3] = lcoords[6];
      sub_sign = LevelSet::sign( (is_on_surf[9] ? 0.0 : ldist[2]+ldist[3]) +
				 (is_on_surf[8] ? 0.0 : ldist[1]+ldist[3]) +
				 (is_on_surf[7] ? 0.0 : ldist[0]+ldist[3]) +
				 (is_on_surf[6] ? 0.0 : ldist[0]+ldist[2]) );
      sub_ids[0] = ((is_on_surf[9] && is_on_surf[8] && is_on_surf[6]) && (zero_sign != sub_sign || sub_degenerate[5])) ? -2 : -1;
      sub_ids[1] = ((is_on_surf[8] && is_on_surf[7] && is_on_surf[6]) && (zero_sign != sub_sign || sub_degenerate[4])) ? -2 : -1;
      sub_ids[2] = lsides[2]; // 9-7-6
      sub_ids[3] = ((is_on_surf[9] && is_on_surf[8] && is_on_surf[7]) && (zero_sign != sub_sign || sub_degenerate[3])) ? -2 : -1;
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // tet #8
  if (!sub_degenerate[7])
    {
      sub_coords[0] = lcoords[5];
      sub_coords[1] = lcoords[6];
      sub_coords[2] = lcoords[4];
      sub_coords[3] = lcoords[8];
      sub_sign = LevelSet::sign( (is_on_surf[5] ? 0.0 : ldist[1]+ldist[2]) +
				 (is_on_surf[6] ? 0.0 : ldist[0]+ldist[2]) +
				 (is_on_surf[4] ? 0.0 : ldist[0]+ldist[1]) +
				 (is_on_surf[8] ? 0.0 : ldist[1]+ldist[3]) );
      sub_ids[0] = ((is_on_surf[5] && is_on_surf[6] && is_on_surf[8]) && (zero_sign != sub_sign || sub_degenerate[5])) ? -2 : -1;
      sub_ids[1] = ((is_on_surf[6] && is_on_surf[4] && is_on_surf[8]) && (zero_sign != sub_sign || sub_degenerate[4])) ? -2 : -1;
      sub_ids[2] = ((is_on_surf[5] && is_on_surf[4] && is_on_surf[8]) && (zero_sign != sub_sign || sub_degenerate[1])) ? -2 : -1;
      sub_ids[3] = lsides[3]; // 4-5-6
      sub = new ContourSubElement_Tet_4( sub_coords, sub_ids, my_owner, my_subelement_depth+1, sub_sign );
      my_subelements.push_back( sub );
    }

  // check quality of subelements
  // Here we assume that the linear tet is always decomposed into reasonable quality sub-tets
  success = true;

  if (krinolog.shouldPrint(LOG_DEBUG))
  {
    for ( unsigned i=0; i<my_subelements.size(); ++i )
    {
      if ( my_subelements[i]->physical_quality() < 0.0 )
      {
        krinolog << "low quality subelement: " << i << "\n"
                << *my_subelements[i] << "\n"
                << "parent:" << "\n"
                << *this << "\n";
      }
    }
  }

  return success;
}

int
ContourSubElement_Tet_4::process_edge( const int i0,
				const int i1,
				const int i2,
				std::array<int,10> & is_on_surf,
				std::array<stk::math::Vector3d,10> & lcoords,
				const std::array<double,4> & ldist )
{
  int edge_node_id = i2;
  is_on_surf[i2] = LevelSet::sign_change( ldist[i0], ldist[i1] );
  if ( is_on_surf[i2] )
    {
      // tolerance chosen very small since degeneracies should already be eliminated
      const double tol = std::numeric_limits<double>::epsilon();

      const double d0 = std::fabs(ldist[i0]);
      const double d1 = std::fabs(ldist[i1]);

      // make calculation completely symmetric
      if (d0 < d1)
	{
	  const double linear_alpha = d0/(d0+d1);

	  if (linear_alpha < tol)
	    {
	      edge_node_id = i0;
	      lcoords[i2] = lcoords[edge_node_id];
	    }
	  else
	    {
	      const double alpha = my_owner->dist_is_linear() ? linear_alpha :
							      find_quadratic_crossing(ldist[i0],ldist[i1],my_owner->distance(0.5*(lcoords[i0]+lcoords[i1])));
	      lcoords[i2] = (1.-alpha) * lcoords[i0] + alpha * lcoords[i1];
	    }
	}
      else
	{
	  const double linear_alpha = d1/(d1+d0);

	  if (linear_alpha < tol)
	    {
	      edge_node_id = i1;
	      lcoords[i2] = lcoords[edge_node_id];
	    }
	  else
	    {
	      const double alpha = my_owner->dist_is_linear() ? linear_alpha :
							      find_quadratic_crossing(ldist[i1],ldist[i0],my_owner->distance(0.5*(lcoords[i1]+lcoords[i0])));
	      lcoords[i2] = (1.-alpha) * lcoords[i1] + alpha * lcoords[i0];
	    }
	}
    }
  else
    {
      // eliminate side node by sliding to one end or the other
      const double d0 = std::fabs(ldist[i0]);
      const double d1 = std::fabs(ldist[i1]);
      const double epsilon = std::sqrt(std::numeric_limits<double>::epsilon()) * (d0 + d1);

      if ( d0 > d1 + epsilon )
	{
	  edge_node_id = i0;
	}
      else if ( d1 > d0 + epsilon )
	{
	  edge_node_id = i1;
	}
      else
	{
	  // tie breaker
	  const stk::math::Vector3d phys0 = my_owner->coordinates(lcoords[i0]);
	  const stk::math::Vector3d phys1 = my_owner->coordinates(lcoords[i1]);

	  if ( is_more(phys1,phys0) )
	    {
	      edge_node_id = i0;
	    }
	  else
	    {
	      edge_node_id = i1;
	    }
	}
      lcoords[i2] = lcoords[edge_node_id];
    }

  return edge_node_id;
}

bool
ContourSubElement_Tet_4::is_degenerate( const std::array<int,10> & edge_node_ids,
				 const int i0, const int i1, const int i2, const int i3 )
{

  // DRN:  This is really ugly, hand-optimized code for looking for degenerate tets
  //       Basically, it checks if any of the edges are degenerate. Then it has to look for degenerate
  //       faces that consist of 3 colinear points.  These are handled by checking against
  //       the 6 specific bad cases.

  if ( edge_node_ids[i0] == edge_node_ids[i1] ||
       edge_node_ids[i0] == edge_node_ids[i2] ||
       edge_node_ids[i0] == edge_node_ids[i3] ||
       edge_node_ids[i1] == edge_node_ids[i2] ||
       edge_node_ids[i1] == edge_node_ids[i3] ||
       edge_node_ids[i2] == edge_node_ids[i3] )
    {
      // this tet is degenerate with two coincident nodes
      return true;
    }

  if ( edge_node_ids[i0]==i0 &&
       edge_node_ids[i1]==i1 &&
       edge_node_ids[i2]==i2 &&
       edge_node_ids[i3]==i3 )
    {
      // this tet is not degenerate since is has no degenerate nodes
      return false;
    }

  // look for a colinear face
  std::vector<int> is_used(10); // initializes to zero (false);
  is_used[edge_node_ids[i0]] = true;
  is_used[edge_node_ids[i1]] = true;
  is_used[edge_node_ids[i2]] = true;
  is_used[edge_node_ids[i3]] = true;

  if ((is_used[0] && ((is_used[1] && is_used[4]) || (is_used[2] && is_used[6]) || (is_used[3] && is_used[7]))) ||
      (is_used[1] && ((is_used[2] && is_used[5]) || (is_used[3] && is_used[8]))) ||
      (is_used[2] && is_used[3] && is_used[9]))
    {
      // this tet has a colinear face
      return true;
    }

  // look for all nodes on the same face
  if ((!is_used[0] && !is_used[4] && !is_used[6] && !is_used[7]) || // all on face 1
      (!is_used[4] && !is_used[1] && !is_used[5] && !is_used[8]) || // all on face 2
      (!is_used[6] && !is_used[5] && !is_used[2] && !is_used[9]) || // all on face 0
      (!is_used[7] && !is_used[8] && !is_used[9] && !is_used[3]))   // all on face 3
    {
      // this tet has all nodes on one face
      return true;
    }

  return false;
}

int
ContourSubElement_Tet_4::side_facets( FacetedSurfaceBase & facets,
			       int side ) const
{
  STK_ThrowAssert( mySideIds[side] == -2 );

  // just one linear facet per linear triangle
  const int num_facets = 1;

  const unsigned * const lnn = get_side_node_ordinals(topology(), side);

  if ( LevelSet::sign_change(0.0, (double) my_sign) )
    facets.emplace_back_3d( my_owner->coordinates(myCoords[lnn[0]]), my_owner->coordinates(myCoords[lnn[1]]), my_owner->coordinates(myCoords[lnn[2]]) );
  else
    facets.emplace_back_3d( my_owner->coordinates(myCoords[lnn[0]]), my_owner->coordinates(myCoords[lnn[2]]), my_owner->coordinates(myCoords[lnn[1]]) );

  return( num_facets );
}

double
ContourSubElement_Tet_4::side_area( int side ) const
{
  STK_ThrowAssert( mySideIds[side] == -2 );

  const unsigned * const lnn = get_side_node_ordinals(topology(), side);

  const std::array<stk::math::Vector3d,3> facetOwnerCoords = { my_owner->coordinates(myCoords[lnn[0]]), my_owner->coordinates(myCoords[lnn[1]]), my_owner->coordinates(myCoords[lnn[2]]) };
  return compute_tri_volume(facetOwnerCoords);
}

ContourSubElement_Tet_10::ContourSubElement_Tet_10(
  const std::array<stk::math::Vector3d,10> & coords,
  const std::array<int,4> &  sideIds,
  const ContourElement * in_owner,
  const int in_subelement_depth,
  const int subelement_sign )
    : ContourSubElementWithTopology<stk::topology::TET_10>( coords,
		  sideIds,
		  in_owner,
		  in_subelement_depth,
		  subelement_sign )
{
  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_conservative_nonlinear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign == 0)
  {
    const std::array<stk::math::Vector3d,4> sub_coords = { myCoords[0], myCoords[1], myCoords[2], myCoords[3] };
    ContourSubElement *sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, mySideIds, my_owner, my_subelement_depth+1 );
    my_subelements.push_back( sub );
  }
}

const int ContourSubElement_Adaptive_Tet_4::MAX_REFINMENT_LEVELS = 6;

ContourSubElement_Adaptive_Tet_4::ContourSubElement_Adaptive_Tet_4(
  const std::array<stk::math::Vector3d,4> & coords,
  const std::array<int,4> &  sideIds,
  const ContourElement * in_owner,
  const int in_subelement_depth )
    : ContourSubElementWithTopology<stk::topology::TET_4>( coords,
                  sideIds,
                  in_owner,
                  in_subelement_depth,
                  0 )
{
  my_edge_age.fill(0);
  non_conformal_decomposition();
}

ContourSubElement_Adaptive_Tet_4::ContourSubElement_Adaptive_Tet_4(
  const std::array<stk::math::Vector3d,4> & coords,
  const std::array<int,4> &  sideIds,
  const std::array<int,6> &  edge_age,
  const ContourElement * in_owner,
  const int in_subelement_depth )
    : ContourSubElementWithTopology<stk::topology::TET_4>( coords,
		  sideIds,
		  in_owner,
		  in_subelement_depth,
		  0 )
{
  my_edge_age = edge_age;
  non_conformal_decomposition();
}

int
ContourSubElement_Adaptive_Tet_4::non_conformal_decomposition()
{
  int success = true; // optimism

  const double snapTol = my_owner->edge_linear_tolerance() * my_owner->length_scale();
  my_sign = ContourElement::compute_conservative_nonlinear_distance_sign(snapTol, myDist.size(), myDist.data());

  if (my_sign != 0)
  {
    return success;
  }

  int longest_bad_edge = -1;

  const std::array<stk::math::Vector3d,10> lcoords = {myCoords[0], myCoords[1], myCoords[2], myCoords[3],
      0.5*(myCoords[0] + myCoords[1]), 0.5*(myCoords[1] + myCoords[2]), 0.5*(myCoords[2] + myCoords[0]),
      0.5*(myCoords[0] + myCoords[3]), 0.5*(myCoords[1] + myCoords[3]), 0.5*(myCoords[2] + myCoords[3])};

  std::array<stk::math::Vector3d,10> lphyscoords;
  for (int n = 0; n < 10; ++n)
    lphyscoords[n] = my_owner->coordinates( lcoords[n] );

  std::array<double,10> ldist = {myDist[0], myDist[1], myDist[2], myDist[3], 0., 0., 0., 0., 0., 0.};
  for (int n = 4; n < 10; ++n)
    ldist[n] = my_owner->distance( lcoords[n] );

  const stk::topology Top = stk::topology::TETRAHEDRON_10;
  int num_edges = Top.num_edges();

  std::vector<int> bad_edges;
  bad_edges.reserve(num_edges);

  std::vector<double> edge_lengths(num_edges);

  for ( int edge = 0; edge < num_edges; edge++ )
    {
      const unsigned * const lnn = get_edge_node_ordinals(Top, edge);

      STK_ThrowAssert(Top.edge_topology(edge).num_nodes() == 3);

      const double edge_straight_length = (lphyscoords[lnn[0]] - lphyscoords[lnn[1]]).length();
      STK_ThrowRequire(edge_straight_length > 0.0);
      edge_lengths[edge] = edge_straight_length;

      const double edge_curve_error = (lphyscoords[lnn[2]] - 0.5*(lphyscoords[lnn[0]] + lphyscoords[lnn[1]])).length();

      const double edge_dist_error = std::fabs(ldist[lnn[2]] - 0.5*(ldist[lnn[0]]+ldist[lnn[1]]));

      const double scale = std::min(std::sqrt(std::numeric_limits<double>::max()),std::fabs(ldist[lnn[0]]) + std::fabs(ldist[lnn[1]]) + my_owner->length_scale());

      const double edge_error = (edge_curve_error + edge_dist_error)*edge_straight_length/(scale*scale);

      if (edge_error > my_owner->edge_nonlinear_tolerance() && my_edge_age[edge] < MAX_REFINMENT_LEVELS)
        {
          bad_edges.push_back(edge);
        }
    }

  double max_length = 0.0;
  for (auto edge : bad_edges)
  {
    const double edge_length = edge_lengths[edge];
    STK_ThrowRequire(edge_length > 0.0);

    // we need an absolute mechanism for selecting the edge to bisect so that all elements that share
    // common edges will make the same decisions
    if (utility::is_more(edge_length,max_length))
    {
      longest_bad_edge = edge;
      max_length = edge_length;
    }
    else if (!utility::is_less(edge_length,max_length)) // tie breaker
    {
      const stk::math::Vector3d & edge_midside_coords = lphyscoords[get_edge_node_ordinals(Top, edge)[2]];
      // note that it is safe to assume that longest_bad_edge is already assigned if edge_length == max_length
      const stk::math::Vector3d longest_edge_midside_coords = lphyscoords[get_edge_node_ordinals(Top, longest_bad_edge)[2]];

      STK_ThrowAssert((utility::is_not_equal(edge_midside_coords[0],longest_edge_midside_coords[0]) ||
                   utility::is_not_equal(edge_midside_coords[1],longest_edge_midside_coords[1]) ||
                   utility::is_not_equal(edge_midside_coords[2],longest_edge_midside_coords[2])));

      if (utility::is_more(edge_midside_coords[0],longest_edge_midside_coords[0]) ||
          (!utility::is_less(edge_midside_coords[0],longest_edge_midside_coords[0]) &&
           (utility::is_more(edge_midside_coords[1],longest_edge_midside_coords[1]) ||
            (!utility::is_less(edge_midside_coords[1],longest_edge_midside_coords[1]) &&
             (utility::is_more(edge_midside_coords[2],longest_edge_midside_coords[2]))))))
      {
        longest_bad_edge = edge;
        max_length = edge_length;
      }
    }
  }

  if ( longest_bad_edge == -1 )
    {
      // no bad edges

      // use a single nonconformal linear tet subelement
      my_subelements.clear();
      my_subelements.reserve(1);

      ContourSubElement *sub = new ContourSubElement_Tet_4( myCoords, mySideIds, my_owner, my_subelement_depth+1 );
      my_subelements.push_back( sub );
    }
  else
    {
      //
      // create 2, adaptive, 4-noded tetrahedra by cutting the longest_bad_edge
      //

      my_subelements.clear();
      my_subelements.reserve(2);
      ContourSubElement *sub  = NULL;
      std::array<stk::math::Vector3d,4> sub_coords;
      std::array<int, 4> sub_ids;
      std::array<int, 6> sub_edge_age;

      static const unsigned node_permute_0[] = { 0,1,2,3 };
      static const unsigned node_permute_1[] = { 1,2,0,3 };
      static const unsigned node_permute_2[] = { 2,0,1,3 };
      static const unsigned node_permute_3[] = { 3,0,2,1 };
      static const unsigned node_permute_4[] = { 1,3,2,0 };
      static const unsigned node_permute_5[] = { 3,2,1,0 };
      static const unsigned side_permute_0[] = { 0,1,2,3 };
      static const unsigned side_permute_1[] = { 1,2,0,3 };
      static const unsigned side_permute_2[] = { 2,0,1,3 };
      static const unsigned side_permute_3[] = { 0,3,1,2 };
      static const unsigned side_permute_4[] = { 0,2,3,1 };
      static const unsigned side_permute_5[] = { 2,3,0,1 };
      static const unsigned edge_permute_0[] = { 0,1,2,3,4,5 };
      static const unsigned edge_permute_1[] = { 1,2,0,4,5,3 };
      static const unsigned edge_permute_2[] = { 2,0,1,5,3,4 };
      static const unsigned edge_permute_3[] = { 3,2,5,4,0,1 };
      static const unsigned edge_permute_4[] = { 4,5,1,0,3,2 };
      static const unsigned edge_permute_5[] = { 5,1,4,3,2,0 };
      static const unsigned * node_permute_table[] = { node_permute_0, node_permute_1, node_permute_2, node_permute_3, node_permute_4, node_permute_5 };
      static const unsigned * side_permute_table[] = { side_permute_0, side_permute_1, side_permute_2, side_permute_3, side_permute_4, side_permute_5 };
      static const unsigned * edge_permute_table[] = { edge_permute_0, edge_permute_1, edge_permute_2, edge_permute_3, edge_permute_4, edge_permute_5 };

      const unsigned * lnn = node_permute_table[longest_bad_edge];
      const unsigned * lsn = side_permute_table[longest_bad_edge];
      const unsigned * len = edge_permute_table[longest_bad_edge];

      const stk::math::Vector3d edge_node = 0.5 * (myCoords[lnn[0]] + myCoords[lnn[1]]);

      // tet #1
      sub_coords[0] = myCoords[lnn[0]];
      sub_coords[1] = edge_node;
      sub_coords[2] = myCoords[lnn[2]];
      sub_coords[3] = myCoords[lnn[3]];
      sub_ids[0] = mySideIds[lsn[0]];
      sub_ids[1] = -1; /* not on any parent side */
      sub_ids[2] = mySideIds[lsn[2]];
      sub_ids[3] = mySideIds[lsn[3]];
      sub_edge_age[0] = my_edge_age[len[0]]+1;
      sub_edge_age[1] = my_edge_age[len[0]]+1;
      sub_edge_age[2] = my_edge_age[len[2]];
      sub_edge_age[3] = my_edge_age[len[3]];
      sub_edge_age[4] = my_edge_age[len[0]]+1;
      sub_edge_age[5] = my_edge_age[len[5]];
      sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, sub_edge_age, my_owner, my_subelement_depth+1 );
      my_subelements.push_back( sub );

      // tet #2
      sub_coords[0] = edge_node;
      sub_coords[1] = myCoords[lnn[1]];
      sub_coords[2] = myCoords[lnn[2]];
      sub_coords[3] = myCoords[lnn[3]];
      sub_ids[0] = mySideIds[lsn[0]];
      sub_ids[1] = mySideIds[lsn[1]];
      sub_ids[2] = -1; /* not on any parent side */
      sub_ids[3] = mySideIds[lsn[3]];
      sub_edge_age[0] = my_edge_age[len[0]]+1;
      sub_edge_age[1] = my_edge_age[len[1]];
      sub_edge_age[2] = my_edge_age[len[0]]+1;
      sub_edge_age[3] = my_edge_age[len[0]]+1;
      sub_edge_age[4] = my_edge_age[len[4]];
      sub_edge_age[5] = my_edge_age[len[5]];
      sub = new ContourSubElement_Adaptive_Tet_4( sub_coords, sub_ids, sub_edge_age, my_owner, my_subelement_depth+1 );
      my_subelements.push_back( sub );
    }

  return success;
}

} // namespace krino
