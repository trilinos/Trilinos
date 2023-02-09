#include "CDT.h"
#include <fstream>
#include <iostream>

void print_entities(CDT::Triangulation<double>& tri)
{
  std::cout << "number of vertices = " << tri.vertices.size() << std::endl;
  std::cout << "number of triangles = " << tri.triangles.size() << std::endl;
  std::cout << "number of constrained edges = " << tri.fixedEdges.size() << std::endl;
}

void print_tin(const std::string& fname, CDT::Triangulation<double>& tri)
{
  std::ofstream fv(fname + "_verts.txt");

  for (auto& vert : tri.vertices)
    fv << vert.pos.x << " " << vert.pos.y << std::endl;

  fv.close();

  std::ofstream ft(fname + "_triangles.txt");
  for (auto& triI : tri.triangles)
    ft << triI.vertices[0] << " " << triI.vertices[1] << " " << triI.vertices[2] << std::endl;

  ft.close();
}

int main(int argc, char* argv[])
{
  using PT = CDT::V2d<double>;
  std::vector<PT> verts{PT::make(0, 0),       PT::make(1, 0),      PT::make(1, 1),       PT::make(0, 1),
                        PT::make(0.25, 0.25), PT::make(1.0, 0.25), PT::make(0.25, 0.75), PT::make(0.75, 0.75)};

  // the algorithm creates a super-triangle that encloses the given
  // points.  There vertices of this super-triangle are always inserted
  // at the beginning of tri.vertices, and are removed by tri.eraseSuperTriangle
  // Also, when creating edges, the vertex IDs are the indices in the
  // vector verts above, not in tri::verticies (ie. they do not include the
  // offset of the 3 super-triangle vertices
  CDT::Triangulation<double> tri(CDT::FindingClosestPoint::ClosestRandom);
  tri.insertVertices(verts);

  print_entities(tri);
  print_tin("mesh", tri);

  std::cout << "inserting constrained edge" << std::endl;
  int offset = 0;
  std::vector<CDT::Edge> edges{CDT::Edge(4 + offset, 7 + offset), CDT::Edge(3 + offset, 7 + offset)};
  tri.insertEdges(edges);

  // tri.eraseSuperTriangle();
  print_entities(tri);
  print_tin("mesh2", tri);

  return 0;
}
