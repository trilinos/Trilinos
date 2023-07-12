#include "mesh_io.hpp"

#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <vector>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

void print_tin(const std::string& fname, std::shared_ptr<Mesh> mesh)
{
  std::ofstream fv(fname + "_verts.txt");
  fv << std::setprecision(16);

  for (auto& vert : mesh->get_vertices())
    fv << vert->get_point_orig(0).get_x() << " " << vert->get_point_orig(0).get_y() << " "
       << vert->get_point_orig(0).get_z() << "\n";

  fv.close();

  std::ofstream ft(fname + "_triangles.txt");
  MeshEntityPtr verts[MAX_DOWN];
  for (auto& triI : mesh->get_elements())
  {
    get_downward(triI, 0, verts);
    ft << verts[0]->get_id() << " " << verts[1]->get_id() << " " << verts[2]->get_id() << "\n";
  }

  ft.close();
}

void print_tin(const std::string& fname, std::vector<MeshEntityPtr>& els)
{
  // get unique verts
  std::map<MeshEntityPtr, int, MeshEntityCompare> verts;
  int id = 0;
  MeshEntityPtr vertsEl[MAX_DOWN];
  for (auto& el : els)
  {
    int nverts = get_downward(el, 0, vertsEl);
    for (int i = 0; i < nverts; ++i)
      if (verts.count(vertsEl[i]) == 0)
        verts[vertsEl[i]] = id++;
  }

  // get verts in order
  std::vector<MeshEntityPtr> vertsFlat(verts.size());
  for (auto& p : verts)
    vertsFlat[p.second] = p.first;

  std::ofstream fv(fname + "_verts.txt");
  fv << std::setprecision(16);

  for (auto& vert : vertsFlat)
    fv << vert->get_point_orig(0).get_x() << " " << vert->get_point_orig(0).get_y() << " "
       << vert->get_point_orig(0).get_z() << "\n";

  fv.close();

  std::ofstream ft(fname + "_triangles.txt");
  MeshEntityPtr vertsI[MAX_DOWN];
  for (auto& triI : els)
  {
    int ndown = get_downward(triI, 0, vertsI);
    for (int i = 0; i < ndown; ++i)
      ft << verts[vertsI[i]] << " ";
    ft << "\n";
  }

  ft.close();
}

void print_vert_edges(const std::string& fname, std::shared_ptr<Mesh> meshIn)
{
  std::ofstream fv(fname + "_verts.txt");
  fv << std::setprecision(16);

  for (auto& vert : meshIn->get_vertices())
  {
    utils::Point pt = vert->get_point_orig(0);
    fv << pt.x << " " << pt.y << " " << pt.z << "\n";
  }

  fv.close();

  std::ofstream fe(fname + "_edges.txt");
  for (auto& edgeI : meshIn->get_edges())
    if (edgeI)
      fe << edgeI->get_down(0)->get_id() << " " << edgeI->get_down(1)->get_id() << "\n";

  fe.close();
}

void print_vert_edges(const std::string& fname, const std::vector<MeshEntityPtr>& els)
{
  // get unique verts
  std::map<MeshEntityPtr, int> verts;
  int id = 0;
  MeshEntityPtr vertsEl[MAX_DOWN];
  for (auto& el : els)
  {
    int nverts = get_downward(el, 0, vertsEl);
    for (int i = 0; i < nverts; ++i)
      if (verts.count(vertsEl[i]) == 0)
        verts[vertsEl[i]] = id++;
  }

  // get verts in order
  std::vector<MeshEntityPtr> vertsFlat(verts.size());
  for (auto& p : verts)
    vertsFlat[p.second] = p.first;

  // get unique edges
  std::set<MeshEntityPtr, MeshEntityCompare> edges;
  id = 0;
  for (auto& el : els)
    for (int i = 0; i < el->count_down(); ++i)
      if (edges.count(el->get_down(i)) == 0)
        edges.insert(el->get_down(i));

  // print verts
  std::ofstream fv(fname + "_verts.txt");
  fv << std::setprecision(16);

  for (auto& vert : vertsFlat)
  {
    utils::Point pt = vert->get_point_orig(0);
    fv << pt.x << " " << pt.y << " " << pt.z << "\n";
  }

  fv.close();

  // print edges
  std::ofstream fe(fname + "_edges.txt");
  for (auto& edgeI : edges)
    fe << verts[edgeI->get_down(0)] << " " << verts[edgeI->get_down(1)] << "\n";

  fe.close();
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
