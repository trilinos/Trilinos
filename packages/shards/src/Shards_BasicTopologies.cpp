// @HEADER
// *****************************************************************************
//                Shards : Shared Discretization Tools
//
// Copyright 2008-2011 NTESS and the Shards contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>

#include <Shards_BasicTopologies.hpp>

extern "C" {
typedef struct CellTopologyData_Subcell      Subcell ;
typedef struct CellTopologyData_Permutation  Permutation ;
}

namespace shards {

enum { MAXIMUM_INDICES = 256 };

const unsigned * index_identity_array()
{
  static unsigned self[ MAXIMUM_INDICES ];
  static int init = 1 ;

  if ( init ) {
    for ( unsigned i = 0 ; i < MAXIMUM_INDICES ; ++i ) { self[i] = i ; }
    init = 0 ;
  }

  return self ;
}

const Subcell * subcell_nodes_array()
{
  static Subcell self[ MAXIMUM_INDICES ];
  static int init = 1 ;

  if ( init ) {

    const CellTopologyData * const top = getCellTopologyData<Node>();
    const unsigned     * const ident = index_identity_array();

    for ( int i = 0 ; i < MAXIMUM_INDICES ; ++i ) {
      self[i].topology = top ;
      self[i].node     = ident + i ;
    }

    init = 0 ;
  }

  return self ;
}

namespace {

template< class IList >
const unsigned * index_list( const IList & )
{
  static const unsigned self[] = {
    (unsigned) IndexListAt< IList ,  0 >::value ,
    (unsigned) IndexListAt< IList ,  1 >::value ,
    (unsigned) IndexListAt< IList ,  2 >::value ,
    (unsigned) IndexListAt< IList ,  3 >::value ,
    (unsigned) IndexListAt< IList ,  4 >::value ,
    (unsigned) IndexListAt< IList ,  5 >::value ,
    (unsigned) IndexListAt< IList ,  6 >::value ,
    (unsigned) IndexListAt< IList ,  7 >::value ,
    (unsigned) IndexListAt< IList ,  8 >::value ,
    (unsigned) IndexListAt< IList ,  9 >::value ,
    (unsigned) IndexListAt< IList , 10 >::value ,
    (unsigned) IndexListAt< IList , 11 >::value ,
    (unsigned) IndexListAt< IList , 12 >::value ,
    (unsigned) IndexListAt< IList , 13 >::value ,
    (unsigned) IndexListAt< IList , 14 >::value ,
    (unsigned) IndexListAt< IList , 15 >::value ,
    (unsigned) IndexListAt< IList , 16 >::value ,
    (unsigned) IndexListAt< IList , 17 >::value ,
    (unsigned) IndexListAt< IList , 18 >::value ,
    (unsigned) IndexListAt< IList , 19 >::value ,
    (unsigned) IndexListAt< IList , 20 >::value ,
    (unsigned) IndexListAt< IList , 21 >::value ,
    (unsigned) IndexListAt< IList , 22 >::value ,
    (unsigned) IndexListAt< IList , 23 >::value ,
    (unsigned) IndexListAt< IList , 24 >::value ,
    (unsigned) IndexListAt< IList , 25 >::value ,
    (unsigned) IndexListAt< IList , 26 >::value ,
    (unsigned) IndexListAt< IList , 27 >::value ,
    (unsigned) IndexListAt< IList , 28 >::value ,
    (unsigned) IndexListAt< IList , 29 >::value ,
    (unsigned) IndexListAt< IList , 30 >::value ,
    (unsigned) IndexListAt< IList , 31 >::value
  };

  return self ;
}

//----------------------------------------------------------------------

template< class TList , class IList , unsigned N >
struct SubcellValue ;

template< class TList , class IList >
struct SubcellValue<TList,IList,0>
{ static void assign( Subcell * ) {} };

template< class TList , class IList , unsigned N >
struct SubcellValue {
  static void assign( Subcell * s )
    {
      enum { I = N - 1 };
      SubcellValue<TList,IList,I>::assign( s );
      s[I].topology = getCellTopologyData< typename TypeListAt<TList,I>::type >();
      s[I].node     = index_list( typename TypeListAt<IList,I>::type() );
    }
};

template< class TList , class IList >
struct SubcellArray {

  enum { N = TypeListLength<TList>::value };

  Subcell array[ N ];

  SubcellArray() { SubcellValue<TList,IList,N>::assign( array ); }
};

//----------------------------------------------------------------------

template< class IList , class PList, unsigned N >
struct PermutationValue ;

template< class IList , class PList>
struct PermutationValue<IList,PList,0>
{ static void assign( Permutation * , Permutation * ) {} };

template< class IList , class PList, unsigned N >
struct PermutationValue {
  static void assign( Permutation * forward , Permutation * inverse )
    {
      enum { I = N - 1 };
      typedef typename TypeListAt<IList,I>::type ForwardType ;
      enum { polarity = IndexListAt<PList,I>::value };
      typedef typename IndexListInverse< ForwardType >::type InverseType ;
      PermutationValue<IList,PList,I>::assign( forward , inverse );
      forward[I].node = index_list( ForwardType() );
      forward[I].polarity = polarity;
      inverse[I].node = index_list( InverseType() );
      inverse[I].polarity = polarity;
    }
};

template< class IList , class PList,
          unsigned NPerm = TypeListLength< IList >::value >
struct PermutationArray ;

template< class IList , class PList>
struct PermutationArray< IList , PList, 0 > {
  Permutation * forward ;
  Permutation * inverse ;
  PermutationArray() : forward(NULL), inverse(NULL) {}
};

template< class IList , class PList, unsigned NPerm >
struct PermutationArray {

  enum { N = NPerm };

  Permutation forward[ N ];
  Permutation inverse[ N ];

  PermutationArray()
  { PermutationValue< IList , PList, N >::assign( forward , inverse ); }
};
  
//----------------------------------------------------------------------

template< class Traits > struct Descriptor ;

template< unsigned Number_Vertex ,
          unsigned Number_Node ,
          class    EdgeList ,
          class    EdgeMaps ,
          class    FaceList ,
          class    FaceMaps >
struct Descriptor<
  CellTopologyTraits< 3 , Number_Vertex , Number_Node ,
                      EdgeList , EdgeMaps , FaceList , FaceMaps ,
                      TypeListEnd > >
{
  typedef CellTopologyTraits< 3 , Number_Vertex , Number_Node ,
                              EdgeList , EdgeMaps ,
                              FaceList , FaceMaps ,
                              TypeListEnd > Traits ;

  typedef SubcellArray< EdgeList , EdgeMaps > EdgeArray ;
  typedef SubcellArray< FaceList , FaceMaps > FaceArray ;

  EdgeArray edges ;
  FaceArray faces ;

  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name )
    : edges(), faces()
    {
      typedef typename Traits::template subcell<0> subcell_0 ;
      typedef typename Traits::template subcell<1> subcell_1 ;
      typedef typename Traits::template subcell<2> subcell_2 ;
      typedef typename Traits::template subcell<3> subcell_3 ;

      self.topology = & top ;
      self.node     = index_identity_array();

      top.base             = base ? base : & top ;
      top.name             = name ;
      top.key              = Traits::key ;
      top.dimension        = 3 ;
      top.vertex_count     = Number_Vertex ;
      top.node_count       = Number_Node ;
      top.edge_count       = EdgeArray::N ;
      top.side_count       = Traits::side_count ;
      top.subcell_homogeneity[0] = subcell_0::homogeneity ;
      top.subcell_homogeneity[1] = subcell_1::homogeneity ;
      top.subcell_homogeneity[2] = subcell_2::homogeneity ;
      top.subcell_homogeneity[3] = subcell_3::homogeneity ;
      top.subcell_count[0] = Number_Node ;
      top.subcell_count[1] = EdgeArray::N ;
      top.subcell_count[2] = FaceArray::N ;
      top.subcell_count[3] = 1 ;
      top.subcell[0]       = subcell_nodes_array();
      top.subcell[1]       = edges.array ;
      top.subcell[2]       = faces.array ;
      top.subcell[3]       = & self ;
      top.side             = faces.array ;
      top.edge             = edges.array ;
    };
};

template< unsigned Number_Vertex ,
          unsigned Number_Node ,
          class    EdgeList ,
          class    EdgeMaps ,
          class    PermutationMaps ,
          class    PermutationPolarity >
struct Descriptor<
  CellTopologyTraits< 2 , Number_Vertex , Number_Node ,
                      EdgeList , EdgeMaps , TypeListEnd , TypeListEnd ,
                      PermutationMaps, PermutationPolarity > >
{
  typedef CellTopologyTraits< 2 , Number_Vertex , Number_Node ,
                              EdgeList , EdgeMaps ,
                              TypeListEnd , TypeListEnd ,
                              PermutationMaps, PermutationPolarity > Traits ;
  
  typedef SubcellArray< EdgeList , EdgeMaps > EdgeArray ;
  typedef PermutationArray< PermutationMaps, PermutationPolarity> PermArray ;

  EdgeArray edges ;
  PermArray perm ;

  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name ) : edges()
    {
      typedef typename Traits::template subcell<0> subcell_0 ;
      typedef typename Traits::template subcell<1> subcell_1 ;
      typedef typename Traits::template subcell<2> subcell_2 ;
      typedef typename Traits::template subcell<3> subcell_3 ;

      self.topology = & top ;
      self.node     = index_identity_array();

      top.base              = base ? base : & top ;
      top.name              = name ;
      top.key               = Traits::key ;
      top.dimension         = 2 ;
      top.vertex_count      = Number_Vertex ;
      top.node_count        = Number_Node ;
      top.edge_count        = EdgeArray::N ;
      top.side_count        = Traits::side_count ;
      top.permutation_count = Traits::permutation_count ;
      top.subcell_homogeneity[0] = subcell_0::homogeneity ;
      top.subcell_homogeneity[1] = subcell_1::homogeneity ;
      top.subcell_homogeneity[2] = subcell_2::homogeneity ;
      top.subcell_homogeneity[3] = subcell_3::homogeneity ;
      top.subcell_count[0]    = Number_Node ;
      top.subcell_count[1]    = EdgeArray::N ;
      top.subcell_count[2]    = 1 ;
      top.subcell_count[3]    = 0 ;
      top.subcell[0]          = subcell_nodes_array();
      top.subcell[1]          = edges.array ;
      top.subcell[2]          = & self ;
      top.subcell[3]          = NULL ;
      top.side                = edges.array ;
      top.edge                = edges.array ;
      top.permutation         = perm.forward ;
      top.permutation_inverse = perm.inverse ;
    };
};

template< unsigned Number_Node , unsigned Number_Vertex ,
          class PermutationMaps , class PermutationPolarity >
struct Descriptor<
  CellTopologyTraits< 1 , Number_Vertex , Number_Node ,
                      TypeListEnd , TypeListEnd ,
                      TypeListEnd , TypeListEnd ,
                      PermutationMaps, PermutationPolarity > >
{
  typedef CellTopologyTraits< 1 , Number_Vertex , Number_Node ,
                              TypeListEnd , TypeListEnd ,
                              TypeListEnd , TypeListEnd ,
                              PermutationMaps, PermutationPolarity > Traits ;
  
  typedef PermutationArray< PermutationMaps, PermutationPolarity > PermArray ;

  PermArray perm ;

  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name )
    {
      typedef typename Traits::template subcell<0> subcell_0 ;
      typedef typename Traits::template subcell<1> subcell_1 ;
      typedef typename Traits::template subcell<2> subcell_2 ;
      typedef typename Traits::template subcell<3> subcell_3 ;

      self.topology = & top ;
      self.node     = index_identity_array();

      top.base              = base ? base : & top ;
      top.name              = name ;
      top.key               = Traits::key ;
      top.dimension         = 1 ;
      top.vertex_count      = Number_Vertex ;
      top.node_count        = Number_Node ;
      top.edge_count        = 0 ;
      top.side_count        = 0 ;
      top.permutation_count = Traits::permutation_count ;
      top.subcell_homogeneity[0] = subcell_0::homogeneity ;
      top.subcell_homogeneity[1] = subcell_1::homogeneity ;
      top.subcell_homogeneity[2] = subcell_2::homogeneity ;
      top.subcell_homogeneity[3] = subcell_3::homogeneity ;
      top.subcell_count[0]    = Number_Node ;
      top.subcell_count[1]    = 1 ;
      top.subcell_count[2]    = 0 ;
      top.subcell_count[3]    = 0 ;
      top.subcell[0]          = subcell_nodes_array();
      top.subcell[1]          = & self ;
      top.subcell[2]          = NULL ;
      top.subcell[3]          = NULL ;
      top.side                = NULL ;
      top.edge                = NULL ;
      top.permutation         = perm.forward ;
      top.permutation_inverse = perm.inverse ;
    };
};

template< unsigned Number_Node , unsigned Number_Vertex>
struct Descriptor<
  CellTopologyTraits< 0 , Number_Node , Number_Vertex ,
                      TypeListEnd , TypeListEnd ,
                      TypeListEnd , TypeListEnd ,
                      TypeListEnd > >
{
  // Two cases: Nodes=Vertices=0 for topo=Node, and Nodes=Vertices=1 for topo=Particle
  static_assert (Number_Node==0 || Number_Node==1,
      "Invalid number of nodes for 0-dimensional topology.");
  static_assert (Number_Node==Number_Vertex,
      "Incompatible number of nodes/vertices for 0-dimensional topology.");

  typedef CellTopologyTraits< 0 , Number_Node , Number_Vertex ,
                              TypeListEnd , TypeListEnd ,
                              TypeListEnd , TypeListEnd ,
                              TypeListEnd > Traits ;
  
  Subcell self ;

  CellTopologyData top ;

  Descriptor( const CellTopologyData * base , const char * name )
    {
      self.topology = & top ;
      self.node     = index_identity_array();

      top.base              = base ? base : & top ;
      top.name              = name ;
      top.key               = Traits::key ;
      top.dimension         = 0 ;
      top.vertex_count      = Number_Vertex ;
      top.node_count        = Number_Node ;
      top.edge_count        = 0 ;
      top.side_count        = 0 ;
      top.permutation_count = 0 ;
      top.subcell_homogeneity[0] = true ;
      top.subcell_homogeneity[1] = false ;
      top.subcell_homogeneity[2] = false ;
      top.subcell_homogeneity[3] = false ;
      top.subcell_count[0]       = 1 ;
      top.subcell_count[1]       = 0 ;
      top.subcell_count[2]       = 0 ;
      top.subcell_count[3]       = 0 ;
      top.subcell[0]             = & self ;
      top.subcell[1]             = NULL ;
      top.subcell[2]             = NULL ;
      top.subcell[3]             = NULL ;
      top.side                   = NULL ;
      top.edge                   = NULL ;
      top.permutation            = NULL ;
      top.permutation_inverse    = NULL ;
    };
};

}

//----------------------------------------------------------------------

std::ostream & operator << ( std::ostream & s , const CellTopologyData & t )
{
  s << t.name ;
  s << " { D = " << t.dimension ;
  s << " , NV = " << t.vertex_count ;
  s << " , K = 0x" << std::hex << t.key << std::dec ;
  s << std::endl ;

  for ( unsigned d = 0 ; d < 4 ; ++d ) {
    for ( unsigned i = 0 ; i < t.subcell_count[d] ; ++i ) {

      const Subcell & sub = t.subcell[d][i] ;

      s << "  subcell[" << d << "][" << i << "] = { " ;

      s << sub.topology->name ;
      s << " ," ;
      for ( unsigned j = 0 ; j < sub.topology->node_count ; ++j ) {
        s << " " << sub.node[j] ;
      }
      s << " }" << std::endl ;
    }
  }

  s << "}" << std::endl ;
  return s ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Node>()
{
  static const char name[] = "Node" ;
  static const Descriptor< Node::Traits > self( NULL , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData< Particle >()
{
  static const char name[] = "Particle" ;
  static const Descriptor< Particle::Traits > self( NULL , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Line<2> >()
{
  static const char name[] = "Line_2" ;
  static const Descriptor< Line<2>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< Line<3> >()
{
  static const char name[] = "Line_3" ;
  static const Descriptor< Line<3>::Traits >
    self( getCellTopologyData<Line<2> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< Beam<2> >()
{
  static const char name[] = "Beam_2" ;
  static const Descriptor< Beam<2>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< Beam<3> >()
{
  static const char name[] = "Beam_3" ;
  static const Descriptor< Beam<3>::Traits > 
    self( getCellTopologyData<Beam<2> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< ShellLine<2> >()
{
  static const char name[] = "ShellLine_2" ;
  static const Descriptor< ShellLine<2>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData< ShellLine<3> >()
{
  static const char name[] = "ShellLine_3" ;
  static const Descriptor< ShellLine<3>::Traits >
    self( getCellTopologyData< ShellLine<2> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Triangle<3> >()
{
  static const char name[] = "Triangle_3" ;
  static const Descriptor< Triangle<3>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Triangle<4> >()
{
  static const char name[] = "Triangle_4" ;
  static const Descriptor< Triangle<4>::Traits >
    self( getCellTopologyData< Triangle<3> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Triangle<6> >()
{
  static const char name[] = "Triangle_6" ;
  static const Descriptor< Triangle<6>::Traits >
    self( getCellTopologyData< Triangle<3> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<ShellTriangle<3> >()
{
  static const char name[] = "ShellTriangle_3" ;
  static const Descriptor< ShellTriangle<3>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<ShellTriangle<6> >()
{
  static const char name[] = "ShellTriangle_6" ;
  static const Descriptor< ShellTriangle<6>::Traits >
    self( getCellTopologyData< ShellTriangle<3> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Quadrilateral<4> >()
{
  static const char name[] = "Quadrilateral_4" ;
  static const Descriptor< Quadrilateral<4>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Quadrilateral<8> >()
{
  static const char name[] = "Quadrilateral_8" ;
  static const Descriptor< Quadrilateral<8>::Traits >
    self( getCellTopologyData<Quadrilateral<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Quadrilateral<9> >()
{
  static const char name[] = "Quadrilateral_9" ;
  static const Descriptor< Quadrilateral<9>::Traits >
    self( getCellTopologyData<Quadrilateral<4> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<ShellQuadrilateral<4> >()
{
  static const char name[] = "ShellQuadrilateral_4" ;
  static const Descriptor< ShellQuadrilateral<4>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<ShellQuadrilateral<8> >()
{
  static const char name[] = "ShellQuadrilateral_8" ;
  static const Descriptor< ShellQuadrilateral<8>::Traits >
    self( getCellTopologyData<ShellQuadrilateral<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<ShellQuadrilateral<9> >()
{
  static const char name[] = "ShellQuadrilateral_9" ;
  static const Descriptor< ShellQuadrilateral<9>::Traits >
    self( getCellTopologyData<ShellQuadrilateral<4> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Hexahedron<8> >()
{
  static const char name[] = "Hexahedron_8" ;
  static const Descriptor< Hexahedron<8>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Hexahedron<20> >()
{
  static const char name[] = "Hexahedron_20" ;
  static const Descriptor< Hexahedron<20>::Traits >
    self( getCellTopologyData<Hexahedron<8> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Hexahedron<27> >()
{
  static const char name[] = "Hexahedron_27" ;
  static const Descriptor< Hexahedron<27>::Traits >
    self( getCellTopologyData<Hexahedron<8> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<4> >()
{
  static const char name[] = "Tetrahedron_4" ;
  static const Descriptor< Tetrahedron<4>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<10> >()
{
  static const char name[] = "Tetrahedron_10" ;
  static const Descriptor< Tetrahedron<10>::Traits >
    self( getCellTopologyData<Tetrahedron<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<11> >()
{
  static const char name[] = "Tetrahedron_11" ;
  static const Descriptor< Tetrahedron<11>::Traits >
    self( getCellTopologyData<Tetrahedron<4> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Tetrahedron<8> >()
{
  static const char name[] = "Tetrahedron_8" ;
  static const Descriptor< Tetrahedron<8>::Traits >
    self( getCellTopologyData<Tetrahedron<4> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Pyramid<5> >()
{
  static const char name[] = "Pyramid_5" ;
  static const Descriptor< Pyramid<5>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Pyramid<13> >()
{
  static const char name[] = "Pyramid_13" ;
  static const Descriptor< Pyramid<13>::Traits >
    self( getCellTopologyData<Pyramid<5> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Pyramid<14> >()
{
  static const char name[] = "Pyramid_14" ;
  static const Descriptor< Pyramid<14>::Traits >
    self( getCellTopologyData<Pyramid<5> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Wedge<6> >()
{
  static const char name[] = "Wedge_6" ;
  static const Descriptor< Wedge<6>::Traits > self( NULL , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Wedge<15> >()
{
  static const char name[] = "Wedge_15" ;
  static const Descriptor< Wedge<15>::Traits >
    self( getCellTopologyData<Wedge<6> >() , name );
  return & self.top ;
}

template<>
const CellTopologyData * getCellTopologyData<Wedge<18> >()
{
  static const char name[] = "Wedge_18" ;
  static const Descriptor< Wedge<18>::Traits >
    self( getCellTopologyData<Wedge<6> >() , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Pentagon<5> >()
{
  static const char name[] = "Pentagon_5" ;
  static const Descriptor< Pentagon<5>::Traits > self( NULL , name );
  return & self.top ;
}

//----------------------------------------------------------------------

template<>
const CellTopologyData * getCellTopologyData<Hexagon<6> >()
{
  static const char name[] = "Hexagon_6" ;
  static const Descriptor< Hexagon<6>::Traits > self( NULL , name );
  return & self.top ;
}


}//namespace shards

