/* $Header: /src4/opennurbs/example_write/example_write.cpp 6     12/05/06 9:39a Dalelear $ */
/* $NoKeywords: $ */
//
// Copyright (c) 1993-2000 Robert McNeel & Associates. All rights reserved.
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

#include "../example_userdata/example_ud.h"

static bool write_simple_file_example(
            FILE* fp,
            int version,
            ON_TextLog& error_log,
            const char* sNotes,
            const ON_3dmSettings* settings,
            int material_count, const ON_Material* material,  // optional rendering material
            int layer_count,    const ON_Layer* layer,        // optional layer definitions
            int light_count,
            const ON_3dmObjectAttributes* light_attributes, // optional light attributes
            ON_Light*                     light,            // lights
            int object_count,
            const ON_3dmObjectAttributes* object_attributes, // optional object attributes
            ON_Object**                   object      // objects
            )
{
  ONX_Model model;
  int i;
  ON_3dmObjectAttributes attribs;


  // some notes
  model.m_properties.m_Notes.m_notes = sNotes;
  model.m_properties.m_Notes.m_bVisible = (model.m_properties.m_Notes.m_notes.Length() > 0);

  // set revision history information
  model.m_properties.m_RevisionHistory.NewRevision();
  
  // set application information
  model.m_properties.m_Application.m_application_name = "OpenNURBS write_simple_file_example() function";
  model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";


  if ( 0 != settings )
    model.m_settings = *settings;

  if ( 0 != material && material_count > 0 )
  {
    model.m_material_table.Reserve(material_count);
    for ( i = 0; i < material_count; i++ )
      model.m_material_table.Append(material[i]);
  }

  // layer table
  {
    // Each object in the object table (written below)
    // should to be on a defined layer.  There should be
    // at least one layer with layer index 0 in every file.

    // layer table indices begin at 0
    ON_Layer default_layer;
    default_layer.SetLayerIndex(0);
    default_layer.SetLayerName("Default");
    if ( 0 == layer || layer_count <= 0 || layer[0].LayerIndex() != 0 ) {
      layer = &default_layer;
      layer_count = 1;
    }

    model.m_layer_table.Reserve(layer_count);
    for ( i = 0; i < layer_count; i++ ) 
    {
      // check that layer index is correct
      if ( layer[i].LayerIndex() != i ) 
      {
        error_log.Print("error: layer[%d].LayerIndex() == %d\n",i,layer[i].LayerIndex());
        layer_count = i;
        break;
      }
      // check that layer's material index is correct
      if (    layer[i].RenderMaterialIndex() < -1 
           || layer[i].RenderMaterialIndex() >= material_count ) 
      {
        error_log.Print("error: layer[%d].RenderMaterialIndex() == %d\n",i,layer[i].RenderMaterialIndex());
        layer_count = i;
        break;
      }
      model.m_layer_table.Append(layer[i]);
    }
  }

  if ( 0 != light && light_count > 0 )
  {
    for ( i = 0; i < light_count; i++ ) 
    {
      ONX_Model_RenderLight& mrl = model.m_light_table.AppendNew();
      mrl.m_light = light[i];
      if ( light_attributes )
        mrl.m_attributes = light_attributes[i];
    }
  }

  if ( 0 != object && object_count > 0 )
  {
    for ( i = 0; i < object_count; i++ ) 
    {
      // get object attributes and make sure layer and material indices are legit
      if ( object[i] ) 
      {
        ONX_Model_Object& mo = model.m_object_table.AppendNew();
        mo.m_object = object[i];
        mo.m_bDeleteObject = false;
        if ( object_attributes )
          mo.m_attributes = object_attributes[i];
      }
    }
  }

  // archive to write to
  ON_BinaryFile archive( ON::write3dm, fp );

  // Set uuid's, indices, etc.
  model.Polish();
  // writes model to archive
  bool ok = model.Write( archive,
                         version, 
                         __FILE__ " write_simple_file_example() " __DATE__, 
                         &error_log );

  return ok;
}

static bool write_points_example( FILE* fp, int version, ON_TextLog& error_log  )
{
  // example demonstrates how to write a singe points and point clouds
  ONX_Model model;


  // file properties (notes, preview image, revision history, ...)
  {
    // set revision history information
    model.m_properties.m_RevisionHistory.NewRevision();
    
    // set application information
    model.m_properties.m_Application.m_application_name = "OpenNURBS write_points_example() function";
    model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
    model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";

    {
      // OPTIONAL - add some notes
      model.m_properties.m_Notes.m_notes = "This file was made with the OpenNURBS write_points_example() function.";
      model.m_properties.m_Notes.m_bVisible = true;
    }
  }

  // file settings (units, tolerances, views, ...)
  {
    // OPTIONAL - change values from defaults
    model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::meters;
    model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.01;
    model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
    model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%
  }

  // layer table
  {
    // define some layers
    ON_Layer layer[3];

    layer[0].SetLayerName("Default");
    layer[0].SetVisible(true);
    layer[0].SetLocked(false);
    layer[0].SetLayerIndex(0);
    layer[0].SetColor( ON_Color(0,0,0) );

    layer[1].SetLayerName("red points");
    layer[1].SetVisible(true);
    layer[1].SetLocked(false);
    layer[1].SetLayerIndex(1);
    layer[1].SetColor( ON_Color(255,0,0) );

    layer[2].SetLayerName("one blue point");
    layer[2].SetVisible(true);
    layer[2].SetLocked(false);
    layer[2].SetLayerIndex(2);
    layer[2].SetColor( ON_Color(0,0,255) );

    model.m_layer_table.Append(layer[0]);
    model.m_layer_table.Append(layer[1]);
    model.m_layer_table.Append(layer[2]);
  }

  // group table
  {
    // we'll put 2 red and one blue point in a group
    ON_Group group;
    group.SetGroupName("group of points");
    group.SetGroupIndex(0);
    model.m_group_table.Append(group);
  }

  // object table

  // single point at (1,4,5) on default layer
  ON_Point point1(ON_3dPoint( 1.0, 4.0, 5.0 ));
  point1.AttachUserData( new CExampleWriteUserData("write_points_example()-point1") );
  {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = &point1;
    mo.m_bDeleteObject = false; // point1 is on the stack
    mo.m_attributes.m_layer_index = 0;
    mo.m_attributes.m_name = "first point";
  }

  // point "cloud" with 3 points on red point cloud layer
  ON_PointCloud* pointcloud = new ON_PointCloud();
  pointcloud->AppendPoint(ON_3dPoint( 1.0, 6.0, 5.0 ));
  pointcloud->AppendPoint(ON_3dPoint( 1.5, 4.5, 6.0 ));
  pointcloud->AppendPoint(ON_3dPoint( 2.0, 5.0, 7.0 ));

  pointcloud->AttachUserData( new CExampleWriteUserData("write_points_example()-pointcloud") );
  {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = pointcloud;
    mo.m_bDeleteObject = true; // ~ONX_Model will delete pointcloud.
    mo.m_attributes.m_layer_index = 1;
    mo.m_attributes.AddToGroup(0); // put these points in the group
    mo.m_attributes.m_name = "3 points";
  }

  // single point at (3,2,4) on red point layer
  ON_Point point2(ON_3dPoint( 3.0, 2.0, 4.0  ));
  point2.AttachUserData( new CExampleWriteUserData("write_points_example()-point2") );
  {
    ONX_Model_Object& mo = model.m_object_table.AppendNew();
    mo.m_object = &point2;
    mo.m_bDeleteObject = false;
    mo.m_attributes.m_layer_index = 2;
    mo.m_attributes.AddToGroup(0); // put this points in the group
    mo.m_attributes.m_name = "last point";
  }

  ON_BinaryFile archive( ON::write3dm, fp ); // fp = pointer from fopoen(...,"wb")

  // start section comment
  const char* sStartSectionComment = __FILE__ "write_points_example()" __DATE__;

  // Set uuid's, indices, etc.
  model.Polish();
  // writes model to archive
  bool ok = model.Write( archive, version, sStartSectionComment, &error_log );

  return ok;
}



static bool write_curves_example( FILE* fp, int version, ON_TextLog& error_log )
{
  // example demonstrates how to write a NURBS curve, line, and circle
  ONX_Model model;

  // file properties (notes, preview image, revision history, ...)

  // set revision history information
  model.m_properties.m_RevisionHistory.NewRevision();

  // set application information
  model.m_properties.m_Application.m_application_name = "OpenNURBS write_curves_example() function";
  model.m_properties.m_Application.m_application_URL = "http://www.opennurbs.org";
  model.m_properties.m_Application.m_application_details = "Example program in OpenNURBS toolkit.";

  // some notes
  model.m_properties.m_Notes.m_notes = "This file was made with the OpenNURBS write_curves_example() function.";
  model.m_properties.m_Notes.m_bVisible = true;


  // file settings (units, tolerances, views, ...)
  model.m_settings.m_ModelUnitsAndTolerances.m_unit_system = ON::inches;
  model.m_settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.001;
  model.m_settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
  model.m_settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%


  // layer table
  {
    // OPTIONAL - define some layers
    ON_Layer layer[3];

    layer[0].SetLayerName("Default");
    layer[0].SetVisible(true);
    layer[0].SetLocked(false);
    layer[0].SetColor( ON_Color(0,0,0) );

    layer[1].SetLayerName("green NURBS wiggle");
    layer[1].SetVisible(true);
    layer[1].SetLocked(false);
    layer[1].SetLayerIndex(1);
    layer[1].SetColor( ON_Color(0,255,0) );

    layer[2].SetLayerName("blue circles");
    layer[2].SetVisible(true);
    layer[2].SetLocked(false);
    layer[2].SetLayerIndex(2);
    layer[2].SetColor( ON_Color(0,0,255) );

    model.m_layer_table.Append(layer[0]);
    model.m_layer_table.Append(layer[1]);
    model.m_layer_table.Append(layer[2]);
  }


  // object table
  {

    {
      // write a line on the default layer
      ONX_Model_Object& mo = model.m_object_table.AppendNew();
      mo.m_object = new ON_LineCurve( ON_Line( ON_3dPoint(1.0,2.0,-1.5), ON_3dPoint(5.0,3.0,2.0) ) );
      mo.m_bDeleteObject = true;
      mo.m_attributes.m_layer_index = 0;
      mo.m_attributes.m_name = "straight line curve";

    }

    {
      // write a wiggly cubic curve on the "green NURBS wiggle" layer
      ON_NurbsCurve* wiggle = new ON_NurbsCurve(
        3, // dimension
        false, // true if rational
        4,     // order = degree+1
        6      // number of control vertices
        );
      int i;
      for ( i = 0; i < wiggle->CVCount(); i++ ) {
        ON_3dPoint pt( 2*i, -i, (i-3)*(i-3) ); // pt = some 3d point
        wiggle->SetCV( i, pt );
      }

      // ON_NurbsCurve's have order+cv_count-2 knots.
      wiggle->SetKnot(0, 0.0);
      wiggle->SetKnot(1, 0.0);
      wiggle->SetKnot(2, 0.0);
      wiggle->SetKnot(3, 1.5);
      wiggle->SetKnot(4, 2.3);
      wiggle->SetKnot(5, 4.0);
      wiggle->SetKnot(6, 4.0);
      wiggle->SetKnot(7, 4.0);

      
      if ( wiggle->IsValid() ) 
      {
        ONX_Model_Object& mo = model.m_object_table.AppendNew();
        mo.m_object = wiggle;
        mo.m_bDeleteObject = true;
        mo.m_attributes.m_layer_index = 1;
        mo.m_attributes.m_name = "wiggly cubic curve";
      }
      else
        delete wiggle;
    }

    {
      // write two circles on the "blue circles"
      ONX_Model_Object& circle1 = model.m_object_table.AppendNew();
      circle1.m_object = new ON_ArcCurve( ON_Circle( ON_3dPoint(1.0,2.0,-1.5), 3.0 ) );
      circle1.m_bDeleteObject = true;
      circle1.m_attributes.m_layer_index = 2;
      circle1.m_attributes.m_name = "radius 3 circle";

      ONX_Model_Object& circle2 = model.m_object_table.AppendNew();
      circle2.m_object = new ON_ArcCurve( ON_Circle( ON_3dPoint(1.0,2.0,-1.5), 5.0 ) );
      circle2.m_bDeleteObject = true;
      circle2.m_attributes.m_layer_index = 2;
      circle2.m_attributes.m_name = "radius 5 circle";
    }

  }

  // start section comments
  const char* sStartSectionComment = __FILE__ "write_curves_example()" __DATE__;

  ON_BinaryFile archive( ON::write3dm, fp ); // fp = pointer from fopoen(...,"wb")

  // Set uuid's, indices, etc.
  model.Polish();
  // writes model to archive
  bool ok = model.Write(archive, version, sStartSectionComment, &error_log );

  return ok;
}


static bool write_surfaces_example( FILE* fp, int version )
{
  // example demonstrates how to write a NURBS surface

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  // The code between the comment bands has nothing to do with I/O.
  // It is simply an easy way to get a NURBS surface to write.
  const int bIsRational = false;
  const int dim = 3;
  const int u_degree = 2;
  const int v_degree = 3;
  const int u_cv_count = 3;
  const int v_cv_count = 5;

  // The knot vectors do NOT have the 2 superfluous knots
  // at the start and end of the knot vector.  If you are
  // coming from a system that has the 2 superfluous knots,
  // just ignore them when writing a 3dm file.
  double u_knot[ u_cv_count + u_degree - 1 ];
  double v_knot[ v_cv_count + v_degree - 1 ];

  // make up a quadratic knot vector with no interior knots
  u_knot[0] = u_knot[1] = 0.0;
  u_knot[2] = u_knot[3] = 1.0;

  // make up a cubic knot vector with one simple interior knot
  v_knot[0] = v_knot[1] = v_knot[2] = 0.0;
  v_knot[3] = 1.5;
  v_knot[4] = v_knot[5] = v_knot[6] = 2.0;

  // Rational control points can be in either homogeneous
  // or euclidean form. Non-rational control points do not
  // need to specify a weight.  
  ON_3dPoint CV[u_cv_count][v_cv_count];

  int i, j;
  for ( i = 0; i < u_cv_count; i++ ) {
    for ( j = 0; j < v_cv_count; j++ ) {
      CV[i][j].x = i;
      CV[i][j].y = j;
      CV[i][j].z = i-j;
    }
  }

  // write a line on the default layer
  ON_NurbsSurface nurbs_surface( dim, bIsRational, 
                        u_degree+1, v_degree+1,
                        u_cv_count, v_cv_count );

  for ( i = 0; i < nurbs_surface.KnotCount(0); i++ )
    nurbs_surface.SetKnot( 0, i, u_knot[i] );

  for ( j = 0; j < nurbs_surface.KnotCount(1); j++ )
    nurbs_surface.SetKnot( 1, j, v_knot[j] );

  for ( i = 0; i < nurbs_surface.CVCount(0); i++ ) {
    for ( j = 0; j < nurbs_surface.CVCount(1); j++ ) {
      nurbs_surface.SetCV( i, j, CV[i][j] );
    }
  }

  bool ok = false;
  if ( nurbs_surface.IsValid() ) 
  {
    ON_BinaryFile archive( ON::write3dm, fp );
    ok = ON_WriteOneObjectArchive( archive, version, nurbs_surface );
  }

  return ok;
}


static bool write_mesh_example( FILE* fp, int version )
{
  // example demonstrates how to create and write a mesh

  // create a mesh to write
  // The mesh is a pyramid with 4 triangular sides and a quadranglar 
  // base.  The mesh has 5 vertices and 5 faces.  
  // The side faces share normals at their common vertices.  The
  // quadrangular base has normals different from the side normal.
  // Coincident vertices that have distinct normals must be
  // duplicated in the vertex list.
  //
  // The apex will be at (1,1.5,4) with normal (0,0,1).
  // The base corners will be at (0,0,0), (0,2,0), (2,3,0), (0,3,0).


  bool bHasVertexNormals = true; // we will specify vertex normals
  bool bHasTexCoords = false;    // we will not specify texture coordinates
  const int vertex_count = 5+4;  // 4 duplicates for different base normals
  const int face_count = 5; // 4 triangle sides and a quad base
  ON_Mesh mesh( face_count, vertex_count, bHasVertexNormals, bHasTexCoords);

  // The SetVertex(), SetNormal(), SetTCoord() and SetFace() functions
  // return true if successful and false if input is illegal.  It is
  // a good idea to inspect this returned value.

  // vertex #0: apex location and normal
  mesh.SetVertex( 0, ON_3dPoint(1.0,  1.5,  5.0) );
  mesh.SetVertexNormal( 0, ON_3dVector(0.0,  0.0,  1.0) );

  // vertex #1: SW corner vertex for sides
  mesh.SetVertex( 1, ON_3dPoint(0.0,  0.0,  0.0) );
  mesh.SetVertexNormal( 1, ON_3dVector(-1.0, -1.0,  0.0) ); // set normal will unitize if needed

  // vertex #2: SE corner vertex for sides
  mesh.SetVertex( 2, ON_3dPoint(2.0,  0.0,  0.0) );
  mesh.SetVertexNormal( 2, ON_3dVector(+1.0, -1.0,  0.0) );

  // vertex #3: NE corner vertex for sides
  mesh.SetVertex( 3, ON_3dPoint(2.0,  3.0,  0.0) );
  mesh.SetVertexNormal( 3, ON_3dVector(+1.0, +1.0,  0.0) );

  // vertex #4: NW corner vertex for sides
  mesh.SetVertex( 4, ON_3dPoint(0.0,  3.0,  0.0) );
  mesh.SetVertexNormal( 4, ON_3dVector(-1.0, +1.0,  0.0) );

  // vertex #5: SW corner vertex for base
  mesh.SetVertex( 5, ON_3dPoint(0.0,  0.0,  0.0) ); // == location of v1
  mesh.SetVertexNormal( 5, ON_3dVector(0.0,  0.0, -1.0) );

  // vertex #6: SE corner vertex for base
  mesh.SetVertex( 6, ON_3dPoint(2.0,  0.0,  0.0) ); // == location of v2
  mesh.SetVertexNormal( 6, ON_3dVector(0.0,  0.0, -1.0) );

  // vertex #7: SW corner vertex for base
  mesh.SetVertex( 7, ON_3dPoint(2.0,  3.0,  0.0) ); // == location of v3
  mesh.SetVertexNormal( 7, ON_3dVector(0.0,  0.0, -1.0) );

  // vertex #8: SW corner vertex for base
  mesh.SetVertex( 8, ON_3dPoint(0.0,  3.0,  0.0) ); // == location of v4
  mesh.SetVertexNormal( 8, ON_3dVector(0.0,  0.0, -1.0) );

  // faces have vertices ordered counter-clockwise

  // South side triangle
  mesh.SetTriangle( 0,   1, 2, 0 );

  // East side triangle
  mesh.SetTriangle( 1,   2, 3, 0 );

  // North side triangle
  mesh.SetTriangle( 2,   3, 4, 0 );

  // West side triangle
  mesh.SetTriangle( 3,   4, 1, 0 );

  // last face is quadrangular base
  mesh.SetQuad( 4,   5, 8, 7, 6 );

  //////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////

  bool ok = false;
  if ( mesh.IsValid() ) 
  {
    // Most applications expect vertex normals.
    // If they are not present, ComputeVertexNormals sets
    // them by averaging face normals.
    if ( !mesh.HasVertexNormals() )
      mesh.ComputeVertexNormals();
    ON_BinaryFile archive( ON::write3dm, fp );
    ok = ON_WriteOneObjectArchive( archive, version, mesh );
  }

  return ok;
}

static bool write_mesh_with_material_example( FILE* fp, int version, ON_TextLog& error_log )
{
  // example demonstrates how to create and write a mesh that uses
  // a rendering material.  You may want to study write_mesh_example() 
  // before examining this function.
  //
  // The key to attaching a texture is to create a mesh
  // with texture coordinates, create a material using the
  // texture as a bitmap, and to attach the material to the
  // mesh before writing the mesh.


  bool bHasVertexNormals = false; // we will NOT specify vertex normals
  bool bHasTexCoords = true;      // we will specify texture coordinates
  const int vertex_count = 40;
  const int face_count = 28;
  ON_Mesh mesh( face_count, vertex_count, bHasVertexNormals, bHasTexCoords);

  // The SetVertex(), SetNormal(), SetTextureCoord() and SetQuad() functions
  // return true if successful and false if input is illegal.  It is
  // a good idea to inspect this returned value.

  // cook up a 5 x 8 grid of vertices
  int vertex_index = 0;
  int i, j;
  for ( i = 0; i <= 4; i++ ) {
    for ( j = 0; j <= 7; j++ ) {
      ON_3fPoint v( (float)i, (float)j, (float)(sin( 2.0*3.14159*j/7.0 ) + cos( 3.14159*i/4.0 )) );
      mesh.SetVertex( vertex_index, v ); // 3d location

      // normalized texture coordinate
      double tcoord_u = i/4.0;
      double tcoord_v = j/7.0;
      mesh.SetTextureCoord( vertex_index, tcoord_u, tcoord_v ); // 2d texture coordinate

      vertex_index++;
    }
  }

  // faces have vertices ordered counter-clockwise

  // cook up a 4 x 7 grid of quadrangular faces
  int face_index = 0;
  for ( i = 0; i < 4; i++ ) {
    for ( j = 0; j < 7; j++ ) {
      int vi[4]; // indices of corner vertices
      vi[0] = i*8 + j;  // vertex at "lower left" corner
      vi[1] = vi[0]+8;  // vertex at "lower right" corner
      vi[2] = vi[1]+1;  // vertex at "upper left" corner
      vi[3] = vi[0]+1;  // vertex at "upper right" corner
      mesh.SetQuad( face_index, vi[0], vi[1], vi[2], vi[3] );
      face_index++;
    }
  }


  // rendering material with a texture map.
  ON_Material material;
  material.SetMaterialIndex(0);
  material.SetAmbient(  ON_Color(  40,  40,  40 ) );
  material.SetDiffuse(  ON_Color( 220, 220, 220 ) );
  material.SetEmission( ON_Color(   0,   0,   0 ) );
  material.SetSpecular( ON_Color( 180, 180, 180 ) );

  material.SetShine( 0.35*ON_Material::MaxShine() ); // 0 = flat
                                                        // MaxShine() = shiney

  material.SetTransparency( 0.2 );  // 0 = opaque, 1 = transparent

  // Texture and bump bitmaps can be Windows bitmap (.BMP), Targa (.TGA),
  // JPEG (.JPG), PCX or PNG files.  Version 1 of Rhino will not support
  // filenames using unicode or multibyte character sets.  As soon as
  // Rhino supports these character sets, the const char* filename will
  // changed to const _TCHAR*.

  // For Rhino to find the texture bitmap, the .3dm file and the
  // .bmp file need to be in the same directory.
  ON_Texture texture;
  texture.m_filename = L"example_texture.bmp";
  material.AddTexture( texture );

  // The render material name is a string used to identify rendering
  // materials in RIB, POV, OBJ, ..., files.  In Rhino, the render
  // material name is set with the SetObjectMaterial command and can
  // be viewed in the Info tab of the dialog displayed by the
  // Properties command.
  material.SetMaterialName( L"my render material" );

  bool ok = false;
  if ( mesh.IsValid() ) 
  {
    // Most applications expect vertex normals.
    // If they are not present, ComputeVertexNormals sets
    // them by averaging face normals.
    if ( !mesh.HasVertexNormals() )
      mesh.ComputeVertexNormals();

    ON_Object* object[1];
    object[0] = &mesh;
    ON_3dmObjectAttributes attributes;
    attributes.m_name = "my mesh with material";
    attributes.m_material_index = 0;
    attributes.SetMaterialSource( ON::material_from_object );
    ok = write_simple_file_example( fp, version, error_log,
         "OpenNURBS write_mesh_with_material_example()", // notes
         NULL,               // default settings
         1, &material,       // render material table
         0, NULL,            // layer table
         0,                  // light count
         NULL,               // light attributes
         NULL,               // lights
         1,                  // object count
         &attributes,        // array of object_count attributes
         object              // array of object_count objects
         );
  }

  return ok;
}

static bool write_spot_light_example( FILE* fp, int version, ON_TextLog& error_log )
{
  // create a blue spotlight shining on a white plane

  // white material for surface
  ON_Material material;
  ON_Color white(255,255,255);
  ON_Color black(0,0,0);
  material.SetMaterialIndex(0);
  material.SetAmbient(  white );
  material.SetDiffuse(  white );
  material.SetEmission( black );
  material.SetSpecular( white );
  material.SetShine( 0.35*ON_Material::MaxShine() ); // 0 = flat
                                                        // MaxShine() = shiney
  material.SetTransparency( 0.0 );  // 0 = opaque, 1 = transparent
  material.SetMaterialName( L"white material" );


  // 2 layers
  ON_Layer layer[2];
  layer[0].SetLayerName(L"surfaces");
  layer[0].SetVisible(true);
  layer[0].SetLocked(false);
  layer[0].SetLayerIndex(0);
  layer[0].SetColor( ON_Color(0, 0, 0) );

  layer[1].SetLayerName(L"lights");
  layer[1].SetVisible(true);
  layer[1].SetLocked(false);
  layer[1].SetLayerIndex(1);
  layer[1].SetColor( ON_Color(0, 0, 0) );


  // spotlight
  ON_Light light;
  light.SetLightIndex(0);
  light.SetLocation( ON_3dPoint(2.0, 3.0, 10.0) );
  light.SetDirection( ON_3dVector(-1.0, -1.0, -10.0) );
  light.SetDiffuse( ON_Color( 0, 0, 255 ) );
  light.SetAmbient( ON_Color( 0, 0, 0 ) );
  light.SetSpecular( ON_Color( 0, 0, 255 ) );
  light.SetSpotExponent( 60 );    // 0 = hard, 128 = soft
  light.SetSpotAngleDegrees( 30.0 );
  light.SetStyle(ON::world_spot_light);

  light.SetLightName( "Blue spot light" );

  ON_3dmObjectAttributes light_attributes;
  light_attributes.Default();
  light_attributes.m_layer_index = 1;     // spotlights layer we defined above
  light_attributes.m_name = "Blue spot light";


  // quick and dirty plane
  ON_PlaneSurface plane (ON_xy_plane);
  plane.SetDomain( 0, -10.0, +10.0 );
  plane.SetDomain( 1, -10.0, +10.0 );

  ON_3dmObjectAttributes object_attributes;
  object_attributes.Default();
  object_attributes.m_layer_index = 0; // surfaces layer we defined above
  object_attributes.m_material_index = 0; // white material we defined above
  object_attributes.SetMaterialSource(ON::material_from_object);
  object_attributes.m_name = "20x20 plane";

  ON_Object* object[1];
  object[0] = &plane;
  bool ok = write_simple_file_example( fp, version, error_log,
       "OpenNURBS write_spot_light_example()", // notes
       NULL,               // default settings
       0, NULL,            // render material table
       0, NULL,            // layer table
       1,                  // light_count
       &light_attributes,  // array of light_count attributes
       &light,             // array of light_count lights
       1,                  // object_count
       &object_attributes, // array of object_count attributes
       object              // array of object_count objects
       );

  return ok;
}


static bool write_viewport_example( FILE* fp, int version, ON_TextLog& error_log,
      const ON_Sphere& sphere // sphere containing region to be viewed
       )
{
  // Writes a 7 viewport layout - 3 along the right side,
  // 3 along the left side, and 1 big on in the middle
  // that displays the space inside the sphere.

  // Viewports have a "target" point inside of the view frustum.
  // This target is the center of view rotations.

  // Viewports have a "construction plane".  This plane is
  // (optionally) displayed as a grid.
  //

  ON_3dmSettings settings;
  {
    // OPTIONAL - change values from defaults
    settings.m_ModelUnitsAndTolerances.m_unit_system = ON::millimeters;
    settings.m_ModelUnitsAndTolerances.m_absolute_tolerance = 0.01;
    settings.m_ModelUnitsAndTolerances.m_angle_tolerance = ON_PI/180.0; // radians
    settings.m_ModelUnitsAndTolerances.m_relative_tolerance = 0.01; // 1%
  }

  // reserve room for 7 views
  settings.m_views.Reserve(7);

  // some values needed to fill in view information
  const double pos_x[4]  = {0.0,0.25,0.75,1.0};       // x: 0 = left, 1 = right
  const double pos_y[4]  = {0.0,1.0/3.0,2.0/3.0,1.0}; // y: 0 = top, 1 = bottom
  ON_3dVector camDir;
  double fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far;
  double target_distance;

  fr_left   = -sphere.radius;
  fr_right  =  sphere.radius;
  fr_bottom = -sphere.radius;
  fr_top    =  sphere.radius;
  fr_near         = 2.0*sphere.radius; // Rhino's default
  target_distance = 3.0*sphere.radius;
  fr_far          = 4.0*sphere.radius;

  // view number 1
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    view.m_vp.SetProjection( ON::parallel_view );
    camDir = -ON_zaxis;
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_yaxis );
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "+X+Y parallel";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[0];
    view.m_position.m_wnd_right = pos_x[1];
    view.m_position.m_wnd_top = pos_y[2];
    view.m_position.m_wnd_bottom = pos_y[3];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane = ON_xy_plane;
  }

  // view number 2
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    view.m_vp.SetProjection( ON::parallel_view );
    camDir = ON_yaxis;
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_zaxis );
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "+X+Z parallel";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[0];
    view.m_position.m_wnd_right = pos_x[1];
    view.m_position.m_wnd_top = pos_y[1];
    view.m_position.m_wnd_bottom = pos_y[2];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane = ON_zx_plane;
  }

  // view number 3
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    view.m_vp.SetProjection( ON::parallel_view );
    camDir = -ON_xaxis;
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_zaxis );
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "+Y+Z parallel";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[0];
    view.m_position.m_wnd_right = pos_x[1];
    view.m_position.m_wnd_top = pos_y[0];
    view.m_position.m_wnd_bottom = pos_y[1];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane = ON_yz_plane;
  }

  // view number 4
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    view.m_vp.SetProjection( ON::parallel_view );
    camDir = ON_zaxis;
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_yaxis );
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "-X+Y parallel";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[2];
    view.m_position.m_wnd_right = pos_x[3];
    view.m_position.m_wnd_top = pos_y[2];
    view.m_position.m_wnd_bottom = pos_y[3];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane.CreateFromFrame( ON_origin, -ON_xaxis, ON_yaxis );
  }

  // view number 5
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    view.m_vp.SetProjection( ON::parallel_view );
    camDir = -ON_yaxis;
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_zaxis );
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "-X+Z parallel";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[2];
    view.m_position.m_wnd_right = pos_x[3];
    view.m_position.m_wnd_top = pos_y[1];
    view.m_position.m_wnd_bottom = pos_y[2];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane.CreateFromFrame( ON_origin, -ON_xaxis, ON_zaxis );
  }

  // view number 6
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    view.m_vp.SetProjection( ON::parallel_view );
    camDir = ON_xaxis;
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_zaxis );
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "-Y+Z parallel";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[2];
    view.m_position.m_wnd_right = pos_x[3];
    view.m_position.m_wnd_top = pos_y[0];
    view.m_position.m_wnd_bottom = pos_y[1];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane.CreateFromFrame( ON_origin, -ON_yaxis, ON_zaxis );
  }

  // view number 7
  {
    ON_3dmView& view = settings.m_views.AppendNew();

    // set primary view transformation information first
    target_distance = 10.0*sphere.radius;
    const double tan_half_angle = sphere.radius / target_distance;

    view.m_vp.SetProjection( ON::perspective_view );
    camDir = ON_3dVector(-40.0,75.0,-50.0);
    view.m_vp.SetCameraLocation( sphere.Center() - target_distance*camDir );
    view.m_vp.SetCameraDirection( camDir );
    view.m_vp.SetCameraUp( ON_zaxis );
    fr_near   = (target_distance - sphere.radius)/10.0;
    fr_far    = target_distance + 1.5*sphere.radius;
    double d  = fr_near*tan_half_angle;
    fr_left   = -d;
    fr_right  =  d;
    fr_bottom = -d;
    fr_top    =  d;
    view.m_vp.SetFrustum( fr_left, fr_right, fr_bottom, fr_top, fr_near, fr_far );

    view.SetTargetPoint(sphere.Center());

    // secondary view "fluff"
    view.m_name = "skew perspective";

    // position of viewport inside main window
    view.m_position.m_wnd_left = pos_x[1];
    view.m_position.m_wnd_right = pos_x[2];
    view.m_position.m_wnd_top = pos_y[0];
    view.m_position.m_wnd_bottom = pos_y[3];

    // construction plane
    view.m_cplane.Default(); // default grid settings
    view.m_cplane.m_plane = ON_xy_plane;
  }


  bool ok = write_simple_file_example( fp, version, error_log,
         "OpenNURBS write_viewport_example()", // notes
         &settings,          // default settings
         0, NULL,            // render material table
         0, NULL,            // layer table
         0,                  // light count
         NULL,               // light attributes
         NULL,               // lights
         0,                  // object count
         NULL,               // array of object_count attributes
         NULL                // array of object_count objects
         );

  return ok;
}

static void make_trimming_curves( ON_Brep& brep, 
                                  const ON_2dPoint& A2, // start point in parameter space
                                  const ON_2dPoint& B2, // end point in parameter space
                                  const ON_3dPoint& A3, // start point in parameter space
                                  const ON_3dPoint& B3  // end point in parameter space
                                  )
{
  ON_LineCurve* p2dCurve = new ON_LineCurve( A2, B2 );
  ON_LineCurve* p3dCurve = new ON_LineCurve( A3, B3 );

  // it is not necessary for the domains of the 2d and 3d curves
  // to match, but it makes it easier to understand the brep
  ON_Interval domain = p3dCurve->Domain();
  p2dCurve->SetDomain( domain.Min(), domain.Max() );

  brep.m_C2.Append(p2dCurve);

  brep.m_C3.Append(p3dCurve);
}


static bool write_trimmed_surface_example( FILE* fp, int version, ON_TextLog& error_log )
{
  // trimmed surfaces are written as a CRhinoBrep that has
  // a single surface and a single CRhinoBrepFace.
  //
  // Trimming loops are simple closed curves and are oriented
  // so that the active portion of the trimmed surface's
  // domain lies to the left of the trimming curves.

  ON_Brep brep;
  ON_2dPoint q;

  // Create a 10x10 plane surface at z=3 with domain [0,1]x[0,1]
  ON_PlaneSurface* pSurface = new ON_PlaneSurface( ON_Plane( ON_3dPoint( 0, 0,3), 
                                                             ON_3dPoint(10,10,3), 
                                                             ON_3dPoint(10, 0,3) ) );
  pSurface->SetDomain(0,0.0,10.0);
  pSurface->SetDomain(1,0.0,10.0);

  // ~ON_Brep() will delete this surface
  const int si = brep.m_S.Count(); // index of surface
  brep.m_S.Append(pSurface);

  // create simple trimming triangle
  ON_2dPoint A2(1.0, 2.0); // parameter space locations of 2d trim corners
  ON_2dPoint B2(9.0, 1.5);
  ON_2dPoint C2(7.0, 8.0);

  ON_3dPoint A3 = pSurface->PointAt(A2.x,A2.y);
  ON_3dPoint B3 = pSurface->PointAt(B2.x,B2.y);
  ON_3dPoint C3 = pSurface->PointAt(C2.x,C2.y);

  make_trimming_curves( brep, A2, B2, A3, B3 ); // creates 2d and 3d curve
  make_trimming_curves( brep, B2, C2, B3, C3 );
  make_trimming_curves( brep, C2, A2, C3, A3 );

  // there are vertices at the 3 corners
  brep.NewVertex( pSurface->PointAt( A2.x, A2.y ) );
  brep.NewVertex( pSurface->PointAt( B2.x, B2.y ) );
  brep.NewVertex( pSurface->PointAt( C2.x, C2.y ) );

  // the vertices are exact since we have lines on a plane
  brep.m_V[0].m_tolerance = 0.0;
  brep.m_V[1].m_tolerance = 0.0;
  brep.m_V[2].m_tolerance = 0.0;

  // there are 3 edges along the sides of the triangle
  brep.NewEdge( brep.m_V[0], brep.m_V[1], 0 ); // start vertex, end vertex, 3d curve index
  brep.NewEdge( brep.m_V[1], brep.m_V[2], 1 ); // start vertex, end vertex, 3d curve index
  brep.NewEdge( brep.m_V[2], brep.m_V[0], 2 ); // start vertex, end vertex, 3d curve index

  // the edges are exact since we have lines on a plane
  brep.m_E[0].m_tolerance = 0.0;
  brep.m_E[1].m_tolerance = 0.0;
  brep.m_E[2].m_tolerance = 0.0;

  // there is 1 face
  ON_BrepFace& face = brep.NewFace( si );

  // outer boundary trimming loops
  ON_BrepLoop& loop = brep.NewLoop( ON_BrepLoop::outer, face );

  // geometrically, loops are made from a contiguous list of 2d parameter space
  // curves that form a simple closed curve.
  brep.NewTrim( brep.m_E[0], false, loop, 0 ); // A to B
  brep.NewTrim( brep.m_E[1], false, loop, 1 ); // B to C
  brep.NewTrim( brep.m_E[2], false, loop, 2 ); // C to A

  // the trims are exact since we have lines on a plane
  q = brep.m_C2[0]->PointAtStart();
  //brep.m_T[0].m_P[0] = pSurface->PointAt(q.x,q.y);
  q = brep.m_C2[0]->PointAtEnd();
  //brep.m_T[0].m_P[1] = pSurface->PointAt(q.x,q.y);
  brep.m_T[0].m_type = ON_BrepTrim::boundary;
  brep.m_T[0].m_tolerance[0] = 0.0;
  brep.m_T[0].m_tolerance[1] = 0.0;

  q = brep.m_C2[0]->PointAtStart();
  //brep.m_T[1].m_P[0] = pSurface->PointAt(q.x,q.y);
  q = brep.m_C2[0]->PointAtEnd();
  //brep.m_T[1].m_P[1] = pSurface->PointAt(q.x,q.y);
  brep.m_T[1].m_type = ON_BrepTrim::boundary;
  brep.m_T[1].m_tolerance[0] = 0.0;
  brep.m_T[1].m_tolerance[1] = 0.0;

  q = brep.m_C2[0]->PointAtStart();
  //brep.m_T[2].m_P[0] = pSurface->PointAt(q.x,q.y);
  q = brep.m_C2[0]->PointAtEnd();
  //brep.m_T[2].m_P[1] = pSurface->PointAt(q.x,q.y);
  brep.m_T[2].m_type = ON_BrepTrim::boundary;
  brep.m_T[2].m_tolerance[0] = 0.0;
  brep.m_T[2].m_tolerance[1] = 0.0;

    // when debugging your code, IsValid(), IsSolid(), IsManifold() are useful
  // to check.

  bool ok = false;
  if ( brep.IsValid( &error_log ) )
  {
    ON_BinaryFile archive(ON::write3dm,fp);
    ok = ON_WriteOneObjectArchive( archive, version, brep );
  }
  
  int bIsManifold, bIsOriented, bHasBoundary, bIsSolid;
  bIsManifold = brep.IsManifold( &bIsOriented, &bHasBoundary );
  bIsSolid = brep.IsSolid();

  return ok;
}

//int main ( int argc, const char* argv[] )
int main ()
{
  bool rc;
  const char* filename;

  ON::Begin();
  // If you want to learn to write b-rep models, first work through
  // this example paying close attention to write_trimmed_surface_example(),
  // then examime example_brep.cpp.

  // The OpenNURBS toolkit will write version 2 and 3 and read
  // version 1, 2 and 3 of the 3DM file format.
  //
  // version 1 is the legacy Rhino I/O tookit format and was used by Rhino 1.x.
  // version 2 is the OpenNURBS format (released 1 July 2000) and is used by Rhino 2.x
  // version 3 is the OpenNURBS format (released 1 November 2002) and is used by Rhino 3.x
  // version 4 is the OpenNURBS format (released September 2006) and is used by Rhino 4.x
  // version 5 is the OpenNURBS format (released September 2009) and is used by Rhino 5.x

  // version to write
  int version = 4; // File can be read by Rhino 4 and Rhino 5
  //int version = 5; // File can be read by Rhino 5


  // errors printed to stdout
  ON_TextLog error_log;

  // messages printed to stdout
  ON_TextLog message_log;

  // errors logged in text file
  //FILE* error_log_fp = ON::OpenFile("error_log.txt","w");
  //ON_TextLog error_log(error_log_fp);

  filename = "my_points.3dm";
  FILE* fp = ON::OpenFile( filename, "wb" );
  rc = write_points_example( fp, version, error_log );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_curves.3dm";
  fp = ON::OpenFile( filename, "wb" );
  rc = write_curves_example( fp, version, error_log );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_surfaces.3dm";
  fp = ON::OpenFile( filename, "wb" );
  rc = write_surfaces_example( fp, version );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_mesh.3dm";
  fp = ON::OpenFile( filename, "wb" );
  rc = write_mesh_example( fp, version );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_mesh_with_material.3dm";
  fp = ON::OpenFile( filename, "wb" );
  rc = write_mesh_with_material_example( fp, version, error_log );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_spot_light.3dm";
  fp = ON::OpenFile( filename, "wb" );
  rc = write_spot_light_example( fp, version, error_log );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_viewports.3dm";
  fp = ON::OpenFile( filename, "wb" );
  // views will display space inside the sphere
  ON_Sphere sphere ( ON_origin, 10.0 );
  rc = write_viewport_example( fp, version, error_log, sphere );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  filename = "my_trimmed_surface.3dm";
  fp = ON::OpenFile( filename, "wb" );
  rc = write_trimmed_surface_example( fp, version, error_log );
  ON::CloseFile( fp );
  if (rc)
    message_log.Print("Successfully wrote %s.\n",filename);
  else
    message_log.Print("Errors while writing %s.\n",filename);

  ON::End();

  return 0;
}
