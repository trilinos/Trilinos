// $Id$
#ifndef radial_inline_mesh_descH
#define radial_inline_mesh_descH

namespace PAMGEN_NEVADA {

class Vector;
class Radial_Inline_Mesh_Desc : public Inline_Mesh_Desc
{
public:
  Radial_Inline_Mesh_Desc(long long dim){dimension = dim;};
  Radial_Inline_Mesh_Desc(){};
  virtual ~Radial_Inline_Mesh_Desc(){};
  virtual long long Set_Up();
  virtual void Calc_Intervals();
  virtual void calculateSize(long long & total_el_count, 
			     long long & total_node_count, 
			     long long & total_edge_count);
  virtual long long Calc_Coord_Vectors();
  Vector calc_coords_periodic(double total_theta,
			      long long i, 
			      long long j, 
			      long long k);
  virtual void Populate_Coords(double * coords,   
		       std::vector<long long> & global_node_vector, 
		       std::map <long long, long long> & global_node_map,
		       long long num_nodes);
};
}// end namespace
#endif
