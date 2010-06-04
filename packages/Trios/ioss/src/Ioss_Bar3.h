/*--------------------------------------------------------------------*/
/*    Copyright 2004 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_Bar3_h
#define SIERRA_Ioss_Bar3_h

#include <Ioss_CodeTypes.h>
#include <Ioss_ElementTopology.h>
#include <string>

// STL Includes
#include <vector>

namespace Ioss {
  class Bar3 : public Ioss::ElementTopology {

  public:
    static void factory();
    ~Bar3();

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

    Ioss::ElementTopology* face_type(int face_number = 0) const;
    Ioss::ElementTopology* edge_type(int edge_number = 0) const;

  protected:
    Bar3();

  private:
    static Bar3 instance_;

    Bar3(const Bar3&); // Do not implement
  };
}

#endif // SIERRA_Ioss_Bar3_h
