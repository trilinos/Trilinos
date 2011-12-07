/* $Header: /src4/opennurbs/example_brep/example_brep.cpp 3     8/29/06 3:34p Dalelear $ */
/* $NoKeywords: $ */
//
// Copyright (c) 1993-1998 Robert McNeel & Associates. All rights reserved.
// Rhinoceros is a registered trademark of Robert McNeel & Assoicates.
//
// THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY.
// ALL IMPLIED WARRANTIES OF FITNESS FOR ANY PARTICULAR PURPOSE AND OF
// MERCHANTABILITY ARE HEREBY DISCLAIMED.
//				
////////////////////////////////////////////////////////////////

// uncomment the "ON_DLL_IMPORTS" define to use opennurbs as a Windows DLL
//#define ON_DLL_IMPORTS
#include "../opennurbs.h"
#include "../examples_linking_pragmas.h"

// This example demonstrates two things:
//
// 1) How to dig through a ON_Brep face  and get at the surface
//    and trimming information.  See TraverseBrepFace() for
//    details.
//
// 2) How to write an OpenNURBS B-rep.  See MakeTwistedCube() for details.

void TraverseBrepFace( 
       const ON_Brep& brep,
       int fi,                  // brep face index
       ON_TextLog& error_log
       )
{
  if ( fi < 0 || fi >= brep.m_F.Count() ) 
  {
    error_log.Print("Invalid face index\n");
    return;
  }

  const ON_BrepFace& face = brep.m_F[fi];

  // pSrf = underlying untrimmed surface
  const ON_Surface* pSrf = NULL;
  if ( face.m_si < 0 || face.m_si >= brep.m_S.Count() )
    error_log.Print("ERROR: invalid brep.m_F[%d].m_si\n", fi );
  else {
    pSrf = brep.m_S[face.m_si];
    if ( !pSrf )
      error_log.Print("ERROR: invalid brep.m_S[%d] is NULL\n", face.m_si );
  }

  // The face is trimmed with one or more trimming loops.
  //
  // All the 2d trimming curves are oriented so that the
  // active region of the trimmed surface lies to the left
  // of the 2d trimming curve.  
  //
  // If face.m_bRev is true, the orientations of the face in
  // the b-rep is opposited the natural parameteric orientation
  // of the surface.

  // loop_count = number of trimming loops on this face (>=1)
  const int loop_count = face.m_li.Count(); 

  int fli; // face's loop index
  for ( fli = 0; fli < loop_count; fli++ ) {
    const int li = face.m_li[fli]; // li = brep loop index
    const ON_BrepLoop& loop = brep.m_L[li];

    // loop_edge_count = number of trimming edges in this loop
    const int loop_trim_count = loop.m_ti.Count();

    int lti; // loop's trim index
    for ( lti = 0; lti < loop_trim_count; lti++ ) {
      const int ti = loop.m_ti[lti]; // ti = brep trim index
      const ON_BrepTrim& trim = brep.m_T[ti];

      //////////////////////////////////////////////////////
      // 2d trimming information
      //
      // Each trim has a 2d parameter space curve.
      const ON_Curve* p2dCurve = NULL;
      const int c2i = trim.m_c2i; // c2i = brep 2d curve index
      if ( c2i < 0 || c2i >= brep.m_C2.Count() ) {
        error_log.Print("ERROR: invalid brep.m_T[%d].m_c2i\n", ti );
      }
      else {
        p2dCurve = brep.m_C2[c2i];
        if ( !p2dCurve )
          error_log.Print("ERROR: invalid brep.m_C2[%d] is NULL\n", c2i );
      }


      //////////////////////////////////////////////////////
      // topology and 3d geometry information
      //

      // Trim starts at v0 and ends at v1.  When the trim
      // is a loop or on a singular surface side, v0i and v1i
      // will be equal.
      //const int v0i = trim.m_vi[0]; // v0i = brep vertex index
      //const int v1i = trim.m_vi[1]; // v1i = brep vertex index
      //const ON_BrepVertex& v0 = brep.m_V[v0i];
      //const ON_BrepVertex& v1 = brep.m_V[v1i];
      // The vX.m_ei[] array contains the brep.m_E[] indices of
      // the edges that begin or end at vX.
      
      const int ei = trim.m_ei;
      if ( ei == -1 ) {
        // This trim lies on a portion of a singular surface side.
        // The vertex indices are still valid and will be equal.
      }
      else {
        // If trim.m_bRev3d is false, the orientations of the 3d edge
        // and the 3d curve obtained by composing the surface and 2d
        // curve agree.
        //
        // If trim.m_bRev3d is true, the orientations of the 3d edge
        // and the 3d curve obtained by composing the surface and 2d
        // curve are opposite.
        const ON_BrepEdge& edge = brep.m_E[ei];
        const int c3i = edge.m_c3i;
        const ON_Curve* p3dCurve = NULL;

        if ( c3i < 0 || c3i >= brep.m_C3.Count() ) {
          error_log.Print("ERROR: invalid brep.m_E[%d].m_c3i\n", ei );
        }
        else {
          p3dCurve = brep.m_C3[c3i];
          if ( !p3dCurve )
            error_log.Print("ERROR: invalid brep.m_C3[%d] is NULL\n", c3i );
        }

        // The edge.m_ti[] array contains the brep.m_T[] indices
        // for the other trims that are joined to this edge.
      }
    }
  }
}

// symbolic vertex index constants to make code more readable
static const int 
  A = 0,
  B = 1,
  C = 2,
  D = 3,
  E = 4,
  F = 5,
  G = 6,
  H = 7;

// symbolic edge index constants to make code more readable
static const int
  AB =  0,
  BC =  1,
  CD =  2,
  AD =  3,
  EF =  4,
  FG =  5,
  GH =  6,
  EH =  7,
  AE =  8,
  BF =  9,
  CG = 10,
  DH = 11;

// symbolic face index constants to make code more readable
static const int
  ABCD =  0,
  BCGF =  1,
  CDHG =  2,
  ADHE =  3,
  ABFE =  4,
  EFGH =  5;

static ON_Curve* TwistedCubeTrimmingCurve(
              const ON_Surface& s,
              int side // 0 = SW to SE
                       // 1 = SE to NE
                       // 2 = NE to NW
                       // 3 = NW to SW
              )
{
  // A trimming curve is a 2d curve whose image lies in the surface's domain.
  // The "active" portion of the surface is to the left of the trimming curve.
  // An outer trimming loop consists of a simple closed curve running 
  // counter-clockwise around the region it trims.

  ON_2dPoint from, to;
  double u0, u1, v0, v1;

  s.GetDomain( 0, &u0, &u1 );
  s.GetDomain( 1, &v0, &v1 );

  switch ( side ) {
  case 0:  // SW to SE
    from.x = u0; from.y = v0;
    to.x   = u1; to.y   = v0;
    break;
  case 1: // SE to NE
    from.x = u1; from.y = v0;
    to.x   = u1; to.y   = v1;
    break;
  case 2: // NE to NW
    from.x = u1; from.y = v1;
    to.x   = u0; to.y   = v1;
    break;
  case 3: // NW to SW
    from.x = u0; from.y = v1;
    to.x   = u0; to.y   = v0;
    break;
  default:
    return 0;
  }

  ON_Curve* c2d = new ON_LineCurve( from, to );
  c2d->SetDomain(0.0,1.0);

  return c2d;
}


static ON_Curve* TwistedCubeEdgeCurve( const ON_3dPoint& from, const ON_3dPoint& to )
{
  // creates a 3d line segment to be used as a 3d curve in a ON_Brep
  ON_Curve* c3d = new ON_LineCurve( from, to );
  c3d->SetDomain( 0.0, 1.0 );
  return c3d;
}

static ON_Surface* TwistedCubeSideSurface( 
                             const ON_3dPoint& SW, const ON_3dPoint& SE,
                             const ON_3dPoint& NE, const ON_3dPoint& NW
                             )
{
  ON_NurbsSurface* pNurbsSurface = new ON_NurbsSurface(
                                        3,     // dimension
                                        false, // not rational
                                        2,     // "u" order
                                        2,     // "v" order
                                        2,     // number of control vertices in "u" dir
                                        2      // number of control vertices in "v" dir
                                        );
  // corner CVs in counter clockwise order starting in the south west
  pNurbsSurface->SetCV( 0,0, SW );
  pNurbsSurface->SetCV( 1,0, SE );
  pNurbsSurface->SetCV( 1,1, NE );
  pNurbsSurface->SetCV( 0,1, NW );
  // "u" knots
  pNurbsSurface->SetKnot( 0,0, 0.0 );
  pNurbsSurface->SetKnot( 0,1, 1.0 );
  // "v" knots
  pNurbsSurface->SetKnot( 1,0, 0.0 );
  pNurbsSurface->SetKnot( 1,1, 1.0 );

  return pNurbsSurface;
}

static void MakeTwistedCubeEdge( ON_Brep& brep,
                                 int vi0, // index of start vertex
                                 int vi1, // index of end vertex
                                 int c3i  // index of 3d curve
                                 )
{
  ON_BrepVertex& v0 = brep.m_V[vi0];
  ON_BrepVertex& v1 = brep.m_V[vi1];
  ON_BrepEdge& edge = brep.NewEdge(v0,v1,c3i);
  edge.m_tolerance = 0.0;  // this simple example is exact - for models with
                           // non-exact data, set tolerance as explained in
                           // definition of ON_BrepEdge.
}

static void MakeTwistedCubeEdges( ON_Brep& brep )
{

  // In this simple example, the edge indices exactly match the 3d
  // curve indices.  In general,the correspondence between edge and
  // curve indices can be arbitrary.  It is permitted for multiple
  // edges to use different portions of the same 3d curve.  The 
  // orientation of the edge always agrees with the natural 
  // parametric orientation of the curve.
  
  // edge that runs from A to B
  MakeTwistedCubeEdge( brep, A, B, AB );
  
  // edge that runs from B to C
  MakeTwistedCubeEdge( brep, B, C, BC );

  // edge that runs from C to D
  MakeTwistedCubeEdge( brep, C, D, CD );

  // edge that runs from A to D
  MakeTwistedCubeEdge( brep, A, D, AD );

  // edge that runs from E to F
  MakeTwistedCubeEdge( brep, E, F, EF );

  // edge that runs from F to G
  MakeTwistedCubeEdge( brep, F, G, FG );

  // edge that runs from G to H
  MakeTwistedCubeEdge( brep, G, H, GH );

  // edge that runs from E to H
  MakeTwistedCubeEdge( brep, E, H, EH );

  // edge that runs from A to E
  MakeTwistedCubeEdge( brep, A, E, AE );

  // edge that runs from B to F
  MakeTwistedCubeEdge( brep, B, F, BF );

  // edge that runs from C to G
  MakeTwistedCubeEdge( brep, C, G, CG );

  // edge that runs from D to H
  MakeTwistedCubeEdge( brep, D, H, DH );
}

static int MakeTwistedCubeTrimmingLoop(  ON_Brep& brep, // returns index of loop
     ON_BrepFace& face,  // face loop is on
     //int vSWi, int vSEi, int vNEi, int vNWi, // Indices of corner vertices listed in SW,SE,NW,NE order
     int eSi,     // index of edge on south side of surface
     int eS_dir,  // orientation of edge with respect to surface trim
     int eEi,     // index of edge on south side of surface
     int eE_dir,  // orientation of edge with respect to surface trim
     int eNi,     // index of edge on south side of surface
     int eN_dir,  // orientation of edge with respect to surface trim
     int eWi,     // index of edge on south side of surface
     int eW_dir   // orientation of edge with respect to surface trim
                                )
{
  const ON_Surface& srf = *brep.m_S[face.m_si];

  ON_BrepLoop& loop = brep.NewLoop( ON_BrepLoop::outer, face );

  // Create trimming curves running counter clockwise around the surface's domain.
  // Start at the south side
  ON_Curve* c2;
  int c2i, ei=0, bRev3d=0;
  ON_2dPoint q;
  ON_Surface::ISO iso = ON_Surface::not_iso;

  for ( int side = 0; side < 4; side++ ) {
    // side: 0=south, 1=east, 2=north, 3=west
    
    c2 = TwistedCubeTrimmingCurve( srf, side );
    c2i = brep.m_C2.Count();
    brep.m_C2.Append(c2);

    switch ( side ) {
    case 0: // south
      ei = eSi;
      bRev3d = (eS_dir == -1);
      iso = ON_Surface::S_iso;
      break;
    case 1: // east
      ei = eEi;
      bRev3d = (eE_dir == -1);
      iso = ON_Surface::E_iso;
      break;
    case 2: // north
      ei = eNi;
      bRev3d = (eN_dir == -1);
      iso = ON_Surface::N_iso;
      break;
    case 3: // west
      ei = eWi;
      bRev3d = (eW_dir == -1);
      iso = ON_Surface::W_iso;
      break;
    }

    ON_BrepTrim& trim = brep.NewTrim( brep.m_E[ei], bRev3d, loop, c2i );
    q = c2->PointAtStart();
    //trim.m_P[0] = srf.PointAt( q.x, q.y );
    q = c2->PointAtEnd();
    //trim.m_P[1] = srf.PointAt( q.x, q.y );
    trim.m_iso = iso;
    trim.m_type = ON_BrepTrim::mated; // This b-rep is closed, so all trims
                                         // have mates.
    trim.m_tolerance[0] = 0.0; // This simple example is exact - for models with
    trim.m_tolerance[1] = 0.0; // non-exact data, set tolerance as explained in
                               // definition of ON_BrepTrim.
  }

  return loop.m_loop_index;
}

static void MakeTwistedCubeFace( ON_Brep& brep,
     int si,      // index of 3d surface
     int s_dir,   // orientation of surface with respect to brep
     //int vSWi, int vSEi, int vNEi, int vNWi, // Indices of corner vertices listed in SW,SE,NW,NE order
     int eSi,     // index of edge on south side of surface
     int eS_dir,  // orientation of edge with respect to surface trim
     int eEi,     // index of edge on south side of surface
     int eE_dir,  // orientation of edge with respect to surface trim
     int eNi,     // index of edge on south side of surface
     int eN_dir,  // orientation of edge with respect to surface trim
     int eWi,     // index of edge on south side of surface
     int eW_dir   // orientation of edge with respect to surface trim
                                )
{
  ON_BrepFace& face = brep.NewFace(si);

  MakeTwistedCubeTrimmingLoop( brep, face,
                //vSWi, vSEi, vNEi, vNWi, 
                eSi, eS_dir, 
                eEi, eE_dir, 
                eNi, eN_dir, 
                eWi, eW_dir 
                );

  face.m_bRev = (s_dir == -1);
}

static void MakeTwistedCubeFaces( ON_Brep& brep )
{

  MakeTwistedCubeFace( brep,
    ABCD,       // Index of surface ABCD
    +1,         // orientation of surface with respect to brep
    //A, B, C, D, // Indices of vertices listed in SW,SE,NW,NE order
    AB,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (AB)
    BC,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (BC)
    CD,+1,      // South side edge and its orientation with respect to
                // to the trimming curve   (CD)
    AD,-1       // South side edge and its orientation with respect to
                // to the trimming curve   (AD)
    );

  MakeTwistedCubeFace( brep,
    BCGF,       // Index of surface BCGF
    -1,         // orientation of surface with respect to brep
    //B, C, G, F, // Indices of vertices listed in SW,SE,NW,NE order
    BC,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (BC)
    CG,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (CG)
    FG,-1,      // South side edge and its orientation with respect to
                // to the trimming curve   (FG)
    BF,-1       // South side edge and its orientation with respect to
                // to the trimming curve   (BF)
    );

  MakeTwistedCubeFace( brep,
    CDHG,       // Index of surface CDHG
    -1,         // orientation of surface with respect to brep
    //C, D, H, G, // Indices of vertices listed in SW,SE,NW,NE order
    CD,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (CD)
    DH,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (DH)
    GH,-1,      // South side edge and its orientation with respect to
                // to the trimming curve   (GH)
    CG,-1       // South side edge and its orientation with respect to
                // to the trimming curve   (CG)
    );

  MakeTwistedCubeFace( brep,
    ADHE,       // Index of surface ADHE
    +1,         // orientation of surface with respect to brep
    //A, D, H, E, // Indices of vertices listed in SW,SE,NW,NE order
    AD,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (AD)
    DH,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (DH)
    EH,-1,      // South side edge and its orientation with respect to
                // to the trimming curve   (EH)
    AE,-1       // South side edge and its orientation with respect to
                // to the trimming curve   (AE)
    );

  MakeTwistedCubeFace( brep,
    ABFE,       // Index of surface ABFE
    -1,         // orientation of surface with respect to brep
    //A, B, F, E, // Indices of vertices listed in SW,SE,NW,NE order
    AB,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (AB)
    BF,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (BF)
    EF,-1,      // South side edge and its orientation with respect to
                // to the trimming curve   (EF)
    AE,-1       // South side edge and its orientation with respect to
                // to the trimming curve   (AE)
    );

  MakeTwistedCubeFace( brep,
    EFGH,       // Index of surface EFGH
    -1,         // orientation of surface with respect to brep
    //E, F, G, H, // Indices of vertices listed in SW,SE,NW,NE order
    EF,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (EF)
    FG,+1,      // South side edge and its orientation with respect to
                // to the trimming curve.  (FG)
    GH,+1,      // South side edge and its orientation with respect to
                // to the trimming curve   (GH)
    EH,-1       // South side edge and its orientation with respect to
                // to the trimming curve   (EH)
    );
}


static ON_Brep* MakeTwistedCube( ON_TextLog& error_log )
{
  // This example demonstrates how to construct a ON_Brep
  // with the topology shown below.
  //
  //
  //             H-------e6-------G
  //            /                /|
  //           / |              / |
  //          /  e7            /  e5
  //         /   |            /   |
  //        /                e10  |
  //       /     |          /     |
  //      e11    E- - e4- -/- - - F
  //     /                /      /
  //    /      /         /      /
  //   D---------e2-----C      e9
  //   |     /          |     /
  //   |    e8          |    /
  //   e3  /            e1  /
  //   |                |  /
  //   | /              | /
  //   |                |/
  //   A-------e0-------B
  //
  //

  ON_3dPoint point[8] = {
    ON_3dPoint(  0.0,  0.0,  0.0 ),  // point A = geometry for vertex 0
    ON_3dPoint( 10.0,  0.0,  0.0 ),  // point B = geometry for vertex 1
    ON_3dPoint( 10.0,  8.0, -1.0 ),  // point C = geometry for vertex 2
    ON_3dPoint(  0.0,  6.0,  0.0 ),  // point D = geometry for vertex 3
    ON_3dPoint(  1.0,  2.0, 11.0 ),  // point E = geometry for vertex 4
    ON_3dPoint( 10.0,  0.0, 12.0 ),  // point F = geometry for vertex 5
    ON_3dPoint( 10.0,  7.0, 13.0 ),  // point G = geometry for vertex 6
    ON_3dPoint(  0.0,  6.0, 12.0 )   // point H = geometry for vertex 7
  };

  ON_Brep* brep = new ON_Brep();

  // create eight vertices located at the eight points
  int vi;
  for ( vi = 0; vi < 8; vi++ ) {
    ON_BrepVertex& v = brep->NewVertex(point[vi]);
    v.m_tolerance = 0.0; // this simple example is exact - for models with
                         // non-exact data, set tolerance as explained in
                         // definition of ON_BrepVertex.
  }

  // Create 3d curve geometry - the orientations are arbitrarily chosen
  // so that the end vertices are in alphabetical order.
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[A], point[B] ) ); // line AB
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[B], point[C] ) ); // line BC
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[C], point[D] ) ); // line CD
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[A], point[D] ) ); // line AD
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[E], point[F] ) ); // line EF
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[F], point[G] ) ); // line FG
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[G], point[H] ) ); // line GH
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[E], point[H] ) ); // line EH
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[A], point[E] ) ); // line AE
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[B], point[F] ) ); // line BF
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[C], point[G] ) ); // line CG
  brep->m_C3.Append( TwistedCubeEdgeCurve( point[D], point[H] ) ); // line DH

  // Create the 12 edges that connect the corners of the cube.
  MakeTwistedCubeEdges( *brep );

  // Create 3d surface geometry - the orientations are arbitrarily chosen so
  // that some normals point into the cube and others point out of the cube.
  brep->m_S.Append( TwistedCubeSideSurface( point[A], point[B], point[C], point[D] ) ); // ABCD
  brep->m_S.Append( TwistedCubeSideSurface( point[B], point[C], point[G], point[F] ) ); // BCGF
  brep->m_S.Append( TwistedCubeSideSurface( point[C], point[D], point[H], point[G] ) ); // CDHG
  brep->m_S.Append( TwistedCubeSideSurface( point[A], point[D], point[H], point[E] ) ); // ADHE
  brep->m_S.Append( TwistedCubeSideSurface( point[A], point[B], point[F], point[E] ) ); // ABFE
  brep->m_S.Append( TwistedCubeSideSurface( point[E], point[F], point[G], point[H] ) ); // EFGH


  // Create the CRhinoBrepFaces
  MakeTwistedCubeFaces( *brep );

  if ( !brep->IsValid() ) 
  {
    error_log.Print("Twisted cube b-rep is not valid.\n");
    delete brep;
    brep = NULL;
  }

  //ON_BOOL32 bIsManifold;
  //ON_BOOL32 bHasBoundary;
  //ON_BOOL32 b = brep->IsManifold( &bIsManifold,&bHasBoundary );

  return brep;
}

//int main( int argc, const char *argv[] )
int main()
{
  ON::Begin();

  ON_TextLog error_log;

  // Before working through this example, you should understand
  // the example_write.cpp example.

  ON_Brep* brep = MakeTwistedCube(error_log);
  if ( !brep )
    return 1;

  ONX_Model model;

  ONX_Model_Object& mo = model.m_object_table.AppendNew();
  mo.m_object = brep;
  mo.m_bDeleteObject = true; // ~ONX_Model will delete brep
  brep = 0;
  mo.m_attributes.m_name = "Twisted b-rep";

  // OPTIONAL - change values from defaults
  model.m_properties.m_Notes.m_notes = "File created by OpenNURBS example_brep.cpp";
  model.m_properties.m_Notes.m_bVisible = true;

  model.m_properties.m_Application.m_application_name 
    = "OpenNURBS example_brep.cpp";
  model.m_properties.m_Application.m_application_URL 
    = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details 
    = "OpenNURBS example showing how to create and write a simple b-rep";


  int version = 4; // File can be read by Rhino 4 and Rhino 5
  //int version = 5; // File can be read by Rhino 5
  model.Polish();
  const char* filename = "my_brep.3dm";
  bool rc = model.Write( filename, 
               version,
               __FILE__ " example_brep.cpp " __DATE__,
               &error_log
               );

  if (rc)
    printf("Wrote %s.\n",filename);
  else
    printf("Errors writing %s.\n",filename);

  ON::End();

  return 0;
}

