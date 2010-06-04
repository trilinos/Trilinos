/* Copyright 2001 Sandia Corporation, Albuquerque, NM. */
#ifndef SIERRA_Ioss_TriShell3_h
#define SIERRA_Ioss_TriShell3_h

#include <Ioss_CodeTypes.h>
#include <Ioss_ElementTopology.h>
#include <string>

// STL Includes
#include <vector>

namespace Ioss {
  class TriShell3 : public Ioss::ElementTopology {

  public:
    static void factory();
    ~TriShell3();

    int spatial_dimension()           const;
    int parametric_dimension()       const;
    bool is_element()                 const {return true;}
    int order()               const;

    int number_corner_nodes() const;
    int number_nodes()        const;
    int number_edges()        const;
    int number_faces()        const;

    int number_nodes_edge(int edge= 0)   const;
    int number_nodes_face(int face= 0)   const;
    int number_edges_face(int face= 0)   const;

    Ioss::IntVector edge_connectivity(int edge_number) const;
    Ioss::IntVector face_connectivity(int face_number) const;
    Ioss::IntVector element_connectivity()             const;

    Ioss::IntVector face_edge_connectivity(int face_number) const;

    Ioss::ElementTopology* face_type(int face_number = 0) const;
    Ioss::ElementTopology* edge_type(int edge_number = 0) const;

  protected:
    TriShell3();

  private:
    static TriShell3 instance_;

    TriShell3(const TriShell3&); // Do not implement
  };
}
#endif // SIERRA_Ioss_TriShell3_h
