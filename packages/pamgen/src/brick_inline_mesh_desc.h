// $Id$
#ifndef brick_inline_mesh_descH
#define brick_inline_mesh_descH

namespace PAMGEN_NEVADA {


class Brick_Inline_Mesh_Desc : public Inline_Mesh_Desc
{
public:
  Brick_Inline_Mesh_Desc(int dim){
    dimension = dim;
  };
  virtual ~Brick_Inline_Mesh_Desc(){};
  virtual int Set_Up();
  virtual void Calc_Intervals();
  virtual void calculateSize(long long & total_el_count, 
			     long long & total_node_count, 
			     long long & total_edge_count);
  virtual int Calc_Coord_Vectors();
  virtual void Populate_Coords(double * coords,   
		       std::vector<int> & global_node_vector, 
		       std::map <int, int> & global_node_map,
		       int num_nodes);
};
}// end namespace PAMGEN_NEVADA
#endif
