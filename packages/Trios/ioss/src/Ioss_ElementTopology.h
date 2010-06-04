/*--------------------------------------------------------------------*/
/*    Copyright 2000 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

// -*- Mode: c++ -*-
#ifndef SIERRA_Ioss_Element_Topology_h
#define SIERRA_Ioss_Element_Topology_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <vector>
#include <map>

namespace Ioss {
  class ElementTopology;
  typedef std::vector<int> IntVector;
  typedef std::vector<std::string> NameList;
  typedef std::map<std::string, ElementTopology*,
    std::less<std::string> > ElementTopologyMap;
  typedef ElementTopologyMap::value_type ETM_VP;

  // ========================================================================
  class ElementTopology {
  public:
    void alias(const std::string& base, const std::string& syn);
    bool is_alias(const std::string &alias) const;

    virtual ~ElementTopology();

    const std::string &name() const {return name_;}

    //: Return the Sierra master element name corresponding to this
    //: element topology.  Somewhat klugy coupling between IO subsystem
    //: and Sierra, but least klugy I could think of...
    std::string master_element_name() const {return masterElementName_;}

    //: Return whether the topology describes an "element". If it
    //: isn't an element, then it is a component of an element.  For
    //example, a quadrilater Shell is an element, but a QuadFace is
    //not.
    //
    // Default implementation returns true if spatial_dimension() ==
    // parametric_dimension(), otherwise returns false;
    // "Structural" elements (shells, rods, trusses, particles) need
    // to override.
    virtual bool is_element() const; 
    virtual int spatial_dimension() const = 0;
    virtual int parametric_dimension() const = 0;
    virtual int order()     const = 0;

    virtual bool edges_similar() const; // true if all edges have same topology
    virtual bool faces_similar() const; // true if all faces have same topology

    virtual int number_corner_nodes() const = 0;
    virtual int number_nodes() const = 0;
    virtual int number_edges() const = 0;
    virtual int number_faces() const = 0;
    int number_boundaries() const;

    virtual int number_nodes_edge(int edge = 0) const = 0;
    virtual int number_nodes_face(int face = 0) const = 0;
    virtual int number_edges_face(int face = 0) const = 0;

    IntVector boundary_connectivity(int edge_number) const;
    virtual IntVector edge_connectivity(int edge_number)     const = 0;
    virtual IntVector face_connectivity(int face_number)     const = 0;
    virtual IntVector element_connectivity()                 const = 0;

    // These have default implementations in ElementTopology.
    // The defaults simply fill in the vector with 0..num.
    // For 'face_edge_connectivity', this is sufficient for 2d
    // elements, 3d need to override.
    // For 'element_edge_connectivity', this works for all elements.
    virtual IntVector face_edge_connectivity(int face_number) const;
    IntVector element_edge_connectivity()             const;

    ElementTopology* boundary_type(int face_number = 0) const;
    virtual ElementTopology* face_type(int face_number = 0)     const = 0;
    virtual ElementTopology* edge_type(int edge_number = 0)     const = 0;

    static ElementTopology* factory(const std::string& type, bool ok_to_fail = false);
    static ElementTopology* factory(unsigned int unique_id);
    static unsigned int get_unique_id(const std::string& type);
    static int describe(NameList *names);

  protected:
    ElementTopology(const std::string& type, const std::string& master_elem_name);

  private:
    ElementTopology(const ElementTopology&); // Do not implement...
    ElementTopology& operator=(const ElementTopology&); // Do not implement...

    const std::string name_;
    const std::string masterElementName_;
    
    static ElementTopologyMap* registry();
  };
}
#endif
