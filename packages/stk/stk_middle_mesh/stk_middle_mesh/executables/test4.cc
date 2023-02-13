#include "CDT.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

using Point = CDT::V2d<double>;

std::vector<Point> read_verts(const std::string& fname)
{
  std::vector<Point> verts;

  std::string fnameFull = fname + "_verts.txt";
  std::ifstream fv(fnameFull);
  if (!fv)
    throw std::invalid_argument("could not open file " + fnameFull);

  double x, y;
  std::string line;
  while (std::getline(fv, line))
  {
    std::istringstream buff(line);
    buff >> x;
    buff >> y;
    verts.push_back(Point::make(x, y));
  }

  return verts;
}

std::vector<CDT::Edge> read_edges(const std::string& fname)
{
  std::vector<CDT::Edge> edges;

  std::string fnameFull = fname + "_edges.txt";
  std::ifstream fv(fnameFull);
  if (!fv)
    throw std::invalid_argument("could not open file " + fnameFull);

  int v1, v2;
  std::string line;
  while (std::getline(fv, line))
  {
    std::istringstream buff(line);
    buff >> v1;
    buff >> v2;
    if (v1 < v2)
      edges.emplace_back(v1, v2);
    else
      edges.emplace_back(v2, v1);
  }

  return edges;
}

double dot(const Point& v1, const Point& v2)
{
  return v1.x * v2.x + v1.y * v2.y;
}

void compute_angles(CDT::VerticesArr3& verts, std::vector<Point>& vertCoords, double angles[3])
{
  int nverts = 3;
  for (int i = 0; i < 3; ++i)
  {
    Point v1 = vertCoords[verts[(i - 1 + nverts) % nverts]];
    Point v2 = vertCoords[verts[i]];
    Point v3 = vertCoords[verts[(i + 1) % nverts]];

    Point r1   = Point::make(v1.x - v2.x, v1.y - v2.y);
    Point r2   = Point::make(v3.x - v2.x, v3.y - v2.y);
    double arg = dot(r1, r2) / (std::sqrt(dot(r1, r1)) * std::sqrt(dot(r2, r2)));
    angles[i]  = std::acos(arg);
  }
}

void print_angles(CDT::Triangulation<double>& tri, std::vector<Point>& vertCoords)
{
  double angles[3];
  double pi = std::atan(1) * 4;
  for (unsigned int i = 0; i < tri.triangles.size(); ++i)
  {
    auto& vertsI = tri.triangles[i].vertices;
    compute_angles(vertsI, vertCoords, angles);
    std::cout << "triangle " << i << " has angles " << angles[0] * 180.0 / pi << ", " << angles[1] * 180.0 / pi << ", "
              << angles[2] * 180.0 / pi << std::endl;
  }
}

void print_tin(const std::string& fname, CDT::Triangulation<double>& tri)
{
  std::ofstream fv(fname + "_verts.txt");
  fv << std::setprecision(16);

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
  auto verts = read_verts("mesh_constraints");
  auto edges = read_edges("mesh_constraints");

  // print input
  std::cout << "input:" << std::endl;
  for (unsigned int i = 0; i < verts.size(); ++i)
    std::cout << "vertex " << i << " = (" << verts[i].x << ", " << verts[i].y << ")" << std::endl;

  for (unsigned int i = 0; i < edges.size(); ++i)
    std::cout << "edge " << i << " has vertices " << edges[i].v1() << ", " << edges[i].v2() << std::endl;

  // triangulate
  CDT::Triangulation<double> tri(CDT::FindingClosestPoint::ClosestRandom);

  tri.insertVertices(verts);
  tri.insertEdges(edges);
  tri.eraseSuperTriangle();

  // output
  std::cout << "\noutput:" << std::endl;
  print_angles(tri, verts);
  print_tin("mesh_output", tri);

  return 0;
}
