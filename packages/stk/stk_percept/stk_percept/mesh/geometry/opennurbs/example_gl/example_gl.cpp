/* $Header: /src4/opennurbs/example_gl/example_gl.cpp 5     8/30/06 2:55p Dalelear $ */
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

#include "../opennurbs_gl.h"

#if defined(ON_COMPILER_MSC)

if ( _MSC_VER != 1400 )
// Tested compilers:
//   Microsoft Visual Studio 2005
//   Support for other compilers is not available.
#error The OpenGL example is not supported on this compiler.
// NOTE:
//   Visual Studio 2005 / 8.0 was the last version of Visual
//   studio to install the libraries and header files for
//   Open GL auxillary functions.
#endif

#include <GL/GLaux.h>   // Open GL auxillary functions

#define ON_EXAMPLE_GL_USE_GLAUX

#elif defined(ON_COMPILER_XCODE)

// Tested compilers:
//   Apple XCode 2.4.1
//   Support for other Apple compilers is not available.
#include <GLUT/glut.h>   // Open GL auxillary functions
#define ON_EXAMPLE_GL_USE_GLUT

#else

// Unsupported compiler:
//   Support for other compilers is not available
#error Choose between OpenGL AUX or OpenGL GLUT.

//#include <GLaux.h>   // Open GL auxillary functions
//#define ON_EXAMPLE_GL_USE_GLAUX

//#include <glut.h>   // Open GL auxillary functions
//#define ON_EXAMPLE_GL_USE_GLUT

#endif


#if defined(_WINDOWS) || defined(WIN32) || defined(_WIN32) || defined(WIN64) || defined(_WIN64)
#define MY_GL_CALLBACK CALLBACK
#define MY_USE_WINDOWS_STUFF
#else
#define MY_GL_CALLBACK
#endif

// Before working on this file, be sure to study the OpenNURBS toolkit 
// file example_read.cpp and to read chapters 1 through 11 of the 
// _Open_GL_Programming_Guide_.
//
// This file contains simple example in modeled after those found in
// the _Open_GL_Programming_Guide_.  The nuts and bolts functions that
// demonstrate how to use Open GL to display OpenNURBS geometry are in
// opennurbs_gl.cpp.

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


class CModel : public ONX_Model
{
public:
  void GetObjectMaterial( int object_index, ON_Material& material ) const;
  ON_3dmView m_view;
  ON_BoundingBox m_bbox;
};

void CModel::GetObjectMaterial( 
          int object_index,
          ON_Material& material 
          ) const
{
  material.Default();
  //const ON_Geometry* geo = 0;

  if ( object_index >= 0 && object_index <= m_object_table.Count() )
  {
    const ONX_Model_Object& mo = m_object_table[object_index];
    if ( 0 != mo.m_object )
    {
      switch( mo.m_object->ObjectType() )
      {
      case ON::surface_object:
      case ON::brep_object:
      case ON::mesh_object:
      case ON::instance_reference:
        GetRenderMaterial( mo.m_attributes, material );
        break;
      default:
        {
          // use emmissive object color for curve objects
          ON_Color c = WireframeColor( mo.m_attributes );
          ON_Color black(0,0,0);
          material.Default();
          material.SetAmbient(black);
          material.SetDiffuse(black);
          material.SetSpecular(black);
          material.SetEmission(c);
        }
        break;
      }
    }
  }
}


void GetDefaultView( const ON_BoundingBox& bbox, ON_3dmView& view )
{
  // simple parallel projection of bounding box;
  double window_height = 1.0;
  double window_width = 1.0;
  double dx, dy, dz;
  double frus_near, frus_far;
  ON_3dPoint camLoc;
  ON_3dVector camDir, camUp;
  view.m_target = 0.5*(bbox.m_min + bbox.m_max);
  dx = 1.1*(bbox.m_max[0] - bbox.m_min[0]);
  dy = 1.1*(bbox.m_max[1] - bbox.m_min[1]);
  dz = 1.1*(bbox.m_max[2] - bbox.m_min[2]);
  if ( dx <= 1.0e-6 && dy <= 1.0e-6 )
    dx = dy = 2.0;
  if ( window_height*dx < window_width*dy ) {
    dx = dy*window_width/window_height;
  }
  else {
    dy = dx*window_height/window_width;
  }
  if ( dz <= 0.1*(dx+dy) )
    dz = 0.1*(dx+dy);
  dx *= 0.5;
  dy *= 0.5;
  dz *= 0.5;


  frus_near = 1.0;
  frus_far = frus_near + 2.0*dz;
  camLoc = view.m_target + (dz + frus_near)*ON_zaxis;
  camDir = -ON_zaxis;
  camUp = ON_yaxis;

  view.m_vp.SetProjection( ON::parallel_view );
  view.m_vp.SetCameraLocation( camLoc );
  view.m_vp.SetCameraDirection( camDir );
  view.m_vp.SetCameraUp( camUp );
  view.m_vp.SetFrustum( -dx, dx, -dy, dy, frus_near, frus_far );
}

///////////////////////////////////////////////////////////////////////
//
// Globals for myDisplay() function passed to auxMainLoop()
//
//////////////////////////////////////////////////////////////////////

// GL display list "name"
static GLuint glb_display_list_number = 1;

// global pointer to active model
CModel* glb_model = 0;

///////////////////////////////////////////////////////////////////////
//
// Functions used in main()
//
//////////////////////////////////////////////////////////////////////

ON_BOOL32 myInitGL( const ON_Viewport&, GLUnurbsObj*& );

void myBuildDisplayList( 
      GLuint,                  // display_list_number,
      GLUnurbsObj*,            // pointer to GL nurbs render
      const CModel&            // geometry to render
      );

extern "C" {
void MY_GL_CALLBACK myNurbsErrorCallback( GLenum ); // for gluNurbsCallback()

void MY_GL_CALLBACK myDisplay( void );              // for auxMainLoop()


void MY_GL_CALLBACK myKeyLeftArrowEvent( void );    // for auxKeyFunc();
void MY_GL_CALLBACK myKeyRightArrowEvent( void );   // for auxKeyFunc();
void MY_GL_CALLBACK myKeyUpArrowEvent( void );      // for auxKeyFunc();
void MY_GL_CALLBACK myKeyDownArrowEvent( void );    // for auxKeyFunc();
void MY_GL_CALLBACK myKeyViewExtents( void );       // for auxKeyFunc();

#if defined(ON_EXAMPLE_GL_USE_GLAUX)
void MY_GL_CALLBACK myGLAUX_Reshape( GLsizei, GLsizei );  // for auxReshapeFunc()

void MY_GL_CALLBACK myGLAUX_MouseLeftEvent( AUX_EVENTREC* );   // for auxMouseFunc();
void MY_GL_CALLBACK myGLAUX_MouseMiddleEvent( AUX_EVENTREC* ); // for auxMouseFunc();
void MY_GL_CALLBACK myGLAUX_MouseRightEvent( AUX_EVENTREC* );  // for auxMouseFunc();

typedef void (CALLBACK* RHINO_GL_NURBS_ERROR)();
#endif

#if defined(ON_EXAMPLE_GL_USE_GLUT)
void MY_GL_CALLBACK myGLUT_Reshape( int, int );  // for glutReshapeFunc()

void MY_GL_CALLBACK myGLUT_MouseEvent( int button, int state, int x, int y );
void MY_GL_CALLBACK myGLUT_KeyboardEvent( unsigned char ch, int x, int y );
void MY_GL_CALLBACK myGLUT_SpecialKeyEvent( int ch, int x, int y );    // for auxKeyFunc();

// If you are using Apple's Xcode and you get a compile error
// on the typedef below, then try using the commented out typedef.
//
// Apple's Xcode 2.4 likes this typedef witht the (...)
typedef void (CALLBACK* RHINO_GL_NURBS_ERROR)(...);
//
// Apple's Xcode 3.2 likes this typedef witht the (...)
//typedef void (CALLBACK* RHINO_GL_NURBS_ERROR)();
#endif

}

///////////////////////////////////////////////////////////////////////
//
// used to set projections
//
void SetGLModelViewMatrix( const ON_Viewport& );
void SetGLProjectionMatrix( ON_Viewport& );
///////////////////////////////////////////////////////////////////////

int main( int argc, const char *argv[] )
{
  // reads model into global glb_model;
  ON::Begin();

  ON_TextLog error_log;

  ON_BOOL32 bOK;
  int window_width  = 500;
  int window_height = 500;
  //double port_aspect = ((double)window_width)/((double)window_height);

  // read the file into model
  if ( argc != 2 ) {
    printf("Syntax: %s filename.3dm\n",argv[0] );
    return 0;
  }
  const char* sFileName = argv[1];
  printf("\nFile:  %s\n", sFileName );

  // read the file
  CModel model;
  if ( !model.Read( sFileName, &error_log ) )
  {
    // read failed
    error_log.Print("Unable to read file %s\n",sFileName);
    return 1;
  }

  glb_model = &model;

  // set bbox = world bounding box of all the objects
  model.m_bbox = model.BoundingBox();
  if ( !model.m_bbox.IsValid() )
  {
    // nothing to look at in this model
    return 2;
  }

  // set model.m_view
  if ( model.m_settings.m_views.Count() > 0 )
  {
    // use first viewport projection in file
    double angle;
    model.m_view.m_vp = model.m_settings.m_views[0].m_vp;
    model.m_view.m_target = model.m_settings.m_views[0].m_target;
    model.m_view.m_vp.GetCameraAngle( &angle );
    model.m_view.m_vp.Extents( angle, model.m_bbox );
  }
  else 
  {
    GetDefaultView( model.m_bbox, model.m_view );
  }

  // If needed, enlarge frustum so its aspect matches the window's aspect.
  // Since the Rhino file does not store the far frustum distance in the
  // file, viewports read from a Rhil file need to have the frustum's far
  // value set by inspecting the bounding box of the geometry to be
  // displayed.

  
  ///////////////////////////////////////////////////////////////////
  //
  // GL stuff starts here
  //
  for(;;) {  
    
#if defined(ON_EXAMPLE_GL_USE_GLAUX)
    wchar_t sWindowTitleString[256];
#endif
#if defined(ON_EXAMPLE_GL_USE_GLUT)
    char sWindowTitleString[256];
#endif
    sWindowTitleString[255] = 0;
    if ( argv[0] && argv[0][0] )
    {
      int i;
      for ( i = 0; i < 254 && argv[0][i]; i++ )
        sWindowTitleString[i] = argv[0][i];
      sWindowTitleString[i] = 0;
    }

#if defined(ON_EXAMPLE_GL_USE_GLAUX)
    auxInitPosition( 0, 0, window_width, window_height );
    auxInitDisplayMode( AUX_SINGLE | AUX_RGB | AUX_DEPTH );
    auxInitWindow( sWindowTitleString );

    // register event handler functions
    auxIdleFunc( 0 );
    auxReshapeFunc( myGLAUX_Reshape );
    auxMouseFunc( AUX_LEFTBUTTON,   AUX_MOUSEDOWN, myGLAUX_MouseLeftEvent );
    auxMouseFunc( AUX_LEFTBUTTON,   AUX_MOUSEUP,   myGLAUX_MouseLeftEvent );
    auxMouseFunc( AUX_MIDDLEBUTTON, AUX_MOUSEDOWN, myGLAUX_MouseMiddleEvent );
    auxMouseFunc( AUX_MIDDLEBUTTON, AUX_MOUSEUP,   myGLAUX_MouseMiddleEvent );
    auxMouseFunc( AUX_RIGHTBUTTON,  AUX_MOUSEDOWN, myGLAUX_MouseRightEvent );
    auxMouseFunc( AUX_RIGHTBUTTON,  AUX_MOUSEUP,   myGLAUX_MouseRightEvent );
    auxKeyFunc( AUX_LEFT,  myKeyLeftArrowEvent );
    auxKeyFunc( AUX_RIGHT, myKeyRightArrowEvent );
    auxKeyFunc( AUX_UP,    myKeyUpArrowEvent );
    auxKeyFunc( AUX_DOWN,  myKeyDownArrowEvent );
    auxKeyFunc( AUX_E,  myKeyViewExtents );
    auxKeyFunc( AUX_e,  myKeyViewExtents );
    auxKeyFunc( AUX_Z,  myKeyViewExtents );
    auxKeyFunc( AUX_z,  myKeyViewExtents );
#endif

#if defined(ON_EXAMPLE_GL_USE_GLUT)
    glutInit(&argc,(char**)argv);
    glutInitWindowPosition( 0, 0);
    glutInitWindowSize( window_width, window_height );
    glutInitDisplayMode( GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH );
    glutCreateWindow( sWindowTitleString );

    // register event handler functions
    glutIdleFunc( 0 );
    glutReshapeFunc( myGLUT_Reshape );
    glutMouseFunc( myGLUT_MouseEvent );
    glutKeyboardFunc( myGLUT_KeyboardEvent );
    glutSpecialFunc( myGLUT_SpecialKeyEvent );
    glutDisplayFunc( myDisplay );
#endif

    // setup model view matrix, GL defaults, and the GL NURBS renderer
    GLUnurbsObj* pTheGLNURBSRender = NULL; // OpenGL NURBS rendering context
    bOK = myInitGL( model.m_view.m_vp, pTheGLNURBSRender );

    if ( bOK ) {
      // build display list
      myBuildDisplayList( glb_display_list_number,
                          pTheGLNURBSRender,
                          model );

      // look at it
#if defined(ON_EXAMPLE_GL_USE_GLAUX)
      auxMainLoop( myDisplay );
#endif

#if defined(ON_EXAMPLE_GL_USE_GLUT)
      glutMainLoop(  );
#endif

    }

    gluDeleteNurbsRenderer( pTheGLNURBSRender );

    break;
  }

  //
  // GL stuff ends here
  //
  ///////////////////////////////////////////////////////////////////

  ON::End();

  return 0;
}

///////////////////////////////////////////////////////////////////////
void SetGLModelViewMatrix( const ON_Viewport& viewport )
{
  ON_GL( viewport ); // updates GL model view matrix
}

void SetGLProjectionMatrix( ON_Viewport& viewport )
{
  int pl, pr, pb, pt;
  viewport.GetScreenPort( &pl, &pr, &pb, &pt, NULL, NULL );
  ON_GL( viewport, pl, pr, pb, pt ); // updates GL projection matrix
}

ON_BOOL32 myInitGL( const ON_Viewport& viewport, GLUnurbsObj*& nobj )
{
  // set the model view transform
  SetGLModelViewMatrix( viewport );

  // this stuff works with MSVC 4.2's Open GL. Changes may be needed for other
  // GLs.
  //ON_Color background_color(0,128,128);
  ON_Color background_color(0,63,127);
  //background_color = glb_model->m_settings.m_RenderSettings.m_background_color;
  glClearColor( (float)background_color.FractionRed(), 
                (float)background_color.FractionGreen(), 
                (float)background_color.FractionBlue(), 
                1.0f
                );

  glLightModeli( GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE );
  glDisable( GL_CULL_FACE );
  
  // Rhino viewports have camera "Z" pointing at the camera in a right
  // handed coordinate system.
  glClearDepth( 0.0f );
  glEnable( GL_DEPTH_TEST );
  glDepthFunc( GL_GEQUAL );

  glEnable( GL_LIGHTING );
  glEnable( GL_DITHER );
  //glEnable( GL_AUTO_NORMAL );
  //glEnable( GL_NORMALIZE );

  // default material
  ON_GL( (ON_Material*)NULL );


  // GL rendering of NURBS objects requires a GLUnurbsObj.
  nobj = gluNewNurbsRenderer();
  if ( !nobj )
    return false;
  
  gluNurbsProperty( nobj, GLU_SAMPLING_TOLERANCE,   20.0f );
  gluNurbsProperty( nobj, GLU_PARAMETRIC_TOLERANCE, 0.5f );
  gluNurbsProperty( nobj, GLU_DISPLAY_MODE,         (GLfloat)GLU_FILL );
  //gluNurbsProperty( nobj, GLU_DISPLAY_MODE,         GLU_OUTLINE_POLYGON );
  //gluNurbsProperty( nobj, GLU_DISPLAY_MODE,         GLU_OUTLINE_PATCH );
  gluNurbsProperty( nobj, GLU_SAMPLING_METHOD,      (GLfloat)GLU_PATH_LENGTH );
  //gluNurbsProperty( nobj, GLU_SAMPLING_METHOD,      GLU_PARAMETRIC_ERROR );
  //gluNurbsProperty( nobj, GLU_SAMPLING_METHOD,      GLU_DOMAIN_DISTANCE );
  gluNurbsProperty( nobj, GLU_CULLING,              (GLfloat)GL_FALSE );

  // register GL NURBS error callback
  {
    // hack to get around C vs C++ type checking trauma
    RHINO_GL_NURBS_ERROR fn;
    fn = (RHINO_GL_NURBS_ERROR)myNurbsErrorCallback;
    gluNurbsCallback( nobj, GLU_ERROR, fn );
  }

  return true;
}

///////////////////////////////////////////////////////////////////////

#if defined(ON_EXAMPLE_GL_USE_GLAUX)
void MY_GL_CALLBACK myGLAUX_Reshape( GLsizei w, GLsizei h )
{
  static GLsizei w0 = 0;
  static GLsizei h0 = 0;
  if ( w != w0 || h != h0 ) {
    h0 = h;
    w0 = w;
    ON_GL( glb_model->m_view.m_vp, 0, w-1, h-1, 0 ); // set projection transform
  }
  glViewport( 0, 0, w, h );
}
#endif

#if defined(ON_EXAMPLE_GL_USE_GLUT)
void MY_GL_CALLBACK myGLUT_Reshape( int w, int h )
{
  static int w0 = 0;
  static int h0 = 0;
  if ( w != w0 || h != h0 ) {
    h0 = h;
    w0 = w;
    ON_GL( glb_model->m_view.m_vp, 0, w-1, h-1, 0 ); // set projection transform
  }
  glViewport( 0, 0, w, h );
}
#endif


///////////////////////////////////////////////////////////////////////
static void myRotateView( ON_Viewport& viewport,
                          const ON_3dVector& axis,
                          const ON_3dPoint& center,
                          double angle )
{
  ON_Xform rot;
  ON_3dPoint camLoc;
  ON_3dVector camY, camZ;

  rot.Rotation( angle, axis, center );

  if ( !viewport.GetCameraFrame( camLoc, NULL, camY, camZ ) )
    return;

  camLoc = rot*camLoc;
  camY   = rot*camY;
  camZ   = -(rot*camZ);

  viewport.SetCameraLocation( camLoc );
  viewport.SetCameraDirection( camZ );
  viewport.SetCameraUp( camY );

  ON_GL( viewport ); // update model view
}

static void myRotateLeftRight( ON_Viewport& viewport, double angle )
{
  // ON_3dVector axis = ON_zaxis; // rotate camera about world z axis (z up feel)
  ON_3dVector axis = ON_zaxis; // rotate camera about world y axis (u up feel)
  
  ON_3dPoint center;
  if ( glb_model )
    center = glb_model->m_view.m_target;
  else
    viewport.GetFrustumCenter( center );
  myRotateView( viewport, axis, center, angle );
}

static void myRotateUpDown( ON_Viewport& viewport, double angle )
{
  // rotates camera around the screen x axis
  ON_3dVector camX;
  ON_3dPoint center;
  if ( glb_model )
    center = glb_model->m_view.m_target;
  else
    viewport.GetFrustumCenter( center );
  viewport.GetCameraFrame( NULL, camX, NULL, NULL );
  myRotateView( viewport, camX, center, angle );
}

///////////////////////////////////////////////////////////////////////


void MY_GL_CALLBACK myKeyLeftArrowEvent( void )
{
  myRotateLeftRight( glb_model->m_view.m_vp, ON_PI/12.0 );
}

void MY_GL_CALLBACK myKeyRightArrowEvent( void )
{
  myRotateLeftRight( glb_model->m_view.m_vp, -ON_PI/12.0 );
}

void MY_GL_CALLBACK myKeyUpArrowEvent( void )
{
  myRotateUpDown( glb_model->m_view.m_vp, ON_PI/12.0 );
}

void MY_GL_CALLBACK myKeyDownArrowEvent( void )
{
  myRotateUpDown( glb_model->m_view.m_vp, -ON_PI/12.0 );
}

void MY_GL_CALLBACK myKeyViewExtents( void )
{
  double half_angle = 7.5*ON_PI/180.0;
  glb_model->m_view.m_vp.Extents( half_angle, glb_model->m_bbox );
  SetGLModelViewMatrix( glb_model->m_view.m_vp );
  SetGLProjectionMatrix( glb_model->m_view.m_vp );
}

#if defined(ON_EXAMPLE_GL_USE_GLAUX)

///////////////////////////////////////////////////////////////////////
//
// Mouse event handling
//

static void myGLAUX_MouseEvent( GLint button, const AUX_EVENTREC* event )
{
  static ON_BOOL32 bMouseActive = false;
  static int mx0, my0;
  int mx, my;

  if ( !event ) {
    bMouseActive = false;
    return;
  }

  if ( event->event == AUX_MOUSEDOWN ) {
    if ( bMouseActive ) {
      bMouseActive = false;
      return;
    }
    bMouseActive = true;
    mx0 = event->data[AUX_MOUSEX];
    my0 = event->data[AUX_MOUSEY];
    return;
  }

  if ( !bMouseActive || event->event != AUX_MOUSEUP )
    return;

  mx = event->data[AUX_MOUSEX];
  my = event->data[AUX_MOUSEY];

  switch (button) {
  case AUX_LEFTBUTTON:
    // zoom
    glb_model->m_view.m_vp.ZoomToScreenRect( mx0, my0, mx, my );
    break;
  case AUX_MIDDLEBUTTON:
    break;
  case AUX_RIGHTBUTTON:
    // dolly
    {
      ON_3dVector dolly_vector;
      double d;
      ON_3dPoint camLoc;
      ON_3dVector camZ;
      glb_model->m_view.m_vp.GetCameraFrame( camLoc, NULL, NULL, camZ );
      d = (camLoc-glb_model->m_view.m_target)*camZ;
      if ( glb_model->m_view.m_vp.GetDollyCameraVector(mx0,my0,mx,my,d,dolly_vector) ) {
        glb_model->m_view.m_vp.DollyCamera( dolly_vector );
      }
    }
    break;
  }

  // update GL model view and projection matrices to match viewport changes
  SetGLModelViewMatrix( glb_model->m_view.m_vp );
  SetGLProjectionMatrix( glb_model->m_view.m_vp );

  bMouseActive = false;
}

void MY_GL_CALLBACK myGLAUX_MouseLeftEvent( AUX_EVENTREC* event )
{
  myGLAUX_MouseEvent( AUX_LEFTBUTTON, event );
}

void MY_GL_CALLBACK myGLAUX_MouseMiddleEvent( AUX_EVENTREC* event )
{
  myGLAUX_MouseEvent( AUX_MIDDLEBUTTON, event );
}

void MY_GL_CALLBACK myGLAUX_MouseRightEvent( AUX_EVENTREC* event )
{
  myGLAUX_MouseEvent( AUX_RIGHTBUTTON, event );
}
#endif

#if defined(ON_EXAMPLE_GL_USE_GLUT)

void MY_GL_CALLBACK myGLUT_KeyboardEvent( unsigned char ch, int x, int y )
{
	int m = glutGetModifiers();
	if (m != GLUT_ACTIVE_ALT)
		return;
	if (ch == 'e' || ch == 'z') {
		myKeyViewExtents();
		glutPostRedisplay();
	}
}


void MY_GL_CALLBACK myGLUT_SpecialKeyEvent( int ch, int x, int y )
{
	if (ch == GLUT_KEY_LEFT)
		myKeyLeftArrowEvent();
	if (ch == GLUT_KEY_UP)
		myKeyUpArrowEvent();
	if (ch == GLUT_KEY_RIGHT)
		myKeyRightArrowEvent();
	if (ch == GLUT_KEY_DOWN)
		myKeyDownArrowEvent();
	glutPostRedisplay();
}

void myGLUT_MouseEvent( int button, int state, int x, int y )
{
	static int mx0, my0;
	static int mButton;

	if ( state == GLUT_DOWN ) {
		switch (button) {
		case GLUT_LEFT_BUTTON:
		case GLUT_MIDDLE_BUTTON:
		case GLUT_RIGHT_BUTTON:
			mButton = button;
			mx0 = x;
			my0 = y;
			break;
		}
	}

	if ( state == GLUT_UP && button == mButton ) {
		switch (mButton) {
		case GLUT_LEFT_BUTTON:
			// zoom
			glb_model->m_view.m_vp.ZoomToScreenRect( mx0, my0, x, y );
			break;
		case GLUT_MIDDLE_BUTTON:
			break;
		case GLUT_RIGHT_BUTTON:
			// dolly
			{
				ON_3dVector dolly_vector;
				double d;
				ON_3dPoint camLoc;
				ON_3dVector camZ;
				glb_model->m_view.m_vp.GetCameraFrame( camLoc, NULL, NULL, camZ );
				d = (camLoc-glb_model->m_view.m_target)*camZ;
				if ( glb_model->m_view.m_vp.GetDollyCameraVector(mx0,my0,x,y,d,dolly_vector) ) {
					glb_model->m_view.m_vp.DollyCamera( dolly_vector );
				}
			}
			break;
		}

		// update GL model view and projection matrices to match viewport changes
		SetGLModelViewMatrix( glb_model->m_view.m_vp );
		SetGLProjectionMatrix( glb_model->m_view.m_vp );
		glutPostRedisplay();
	}
}

#endif


///////////////////////////////////////////////////////////////////////

void myDisplayObject( const ON_Object& geometry, const ON_Material& material, GLUnurbsObj* nobj )
{
  // Called from myDisplay() to show geometry.
  // Uses ON_GL() functions found in rhinoio_gl.cpp.
  const ON_Point* point=0;
  const ON_PointCloud* cloud=0;
  const ON_Brep* brep=0;
  const ON_Mesh* mesh=0;
  const ON_Curve* curve=0;
  const ON_Surface* surface=0;

  // specify rendering material
  ON_GL( material );

  brep = ON_Brep::Cast(&geometry);
  if ( brep ) 
  {
    ON_GL(*brep, nobj);
    return;
  }

  mesh = ON_Mesh::Cast(&geometry);
  if ( mesh ) 
  {
    ON_GL(*mesh);
    return;
  }

  curve = ON_Curve::Cast(&geometry);
  if ( curve ) 
  {
    ON_GL( *curve, nobj );
    return;
  }

  surface = ON_Surface::Cast(&geometry);
  if ( surface ) 
  {
    gluBeginSurface( nobj );
    ON_GL( *surface, nobj );
    gluEndSurface( nobj );
    return;
  }

  point = ON_Point::Cast(&geometry);
  if ( point ) 
  {
    ON_GL(*point);
    return;
  }

  cloud = ON_PointCloud::Cast(&geometry);
  if ( cloud ) 
  {
    ON_GL(*cloud);
    return;
  }

}

///////////////////////////////////////////////////////////////////////

void MY_GL_CALLBACK myDisplayLighting( const ON_Viewport&, // viewport, // unreferenced
                                       const CModel& model
                                     )
{
  int light_count = model.m_light_table.Count();
  if ( light_count > 0 ) {
    int maxlighti = light_count;
    if ( maxlighti > GL_MAX_LIGHTS )
      maxlighti = GL_MAX_LIGHTS;
    int lighti;
    for ( lighti = 0; lighti < maxlighti; lighti++ ) {
      ON_GL( model.m_light_table[lighti].m_light, lighti+GL_LIGHT0 );
    }
  }
  else {
    // use default headlight
    // use basic bright white head light with a bit of ambient
    ON_Light head_light;
    head_light.Default();
    ON_GL( head_light, GL_LIGHT0 );
  }
}

///////////////////////////////////////////////////////////////////////

#if defined(MY_USE_WINDOWS_STUFF)
static void myDrawAxesSprite( const ON_Viewport& viewport, HDC hdc )
{
  // Use simple Windows calls to draw world axes sprite in lower left corner.
  // Note that Windows has screen (0,0) in the upper left corner; i.e,
  // screen "y" increases downwards.
  if ( !hdc )
    return;
  const int axes_size = 30;

  int port_left, port_right, port_top, port_bottom;
  if ( !viewport.GetScreenPort( &port_left, &port_right, &port_bottom, &port_top, NULL, NULL ) )
    return;
  const int scr_width  = port_right - port_left; // no "+1" here
  const int scr_height = port_bottom - port_top; // no "+1" here

  if (4*axes_size >= scr_width )
    return;
  if (4*axes_size >= scr_height )
    return;

  int x0 = 3*axes_size/2;
  int y0 = port_bottom - 3*axes_size/2;
  int indx[3] = {0,1,2};
  double scr_coord[3][2];
  viewport.GetCoordinateSprite( axes_size, x0, y0, indx, scr_coord );

#define LXSIZE 3
#define LYSIZE 3
#define LOFF 3

  // draw 3 axes from back to front
  HPEN axis_pen[3];
  axis_pen[0] = CreatePen( PS_SOLID, 2, RGB(255,0,0) );
  axis_pen[1] = CreatePen( PS_SOLID, 2, RGB(0,255,0) );
  axis_pen[2] = CreatePen( PS_SOLID, 2, RGB(0,0,255) );
  HGDIOBJ saved_pen = SelectObject( hdc, axis_pen[0] );

  int i, k, x, y, lx, ly;
  for (i=0;i<3;i++) {
    k = indx[i];
    x = (int)scr_coord[k][0];
    y = (int)scr_coord[k][1];
    // use direction of screen vector to determine letter placement
    lx = x-x0; ly = y-y0;
    if (abs(lx) > abs(ly)) {
      // center letter to right/left of axis end
      lx = (x >= x0) ? x + LXSIZE+LOFF : x - LXSIZE-LOFF;
      ly = y;
    }
    else if (abs(ly) > abs(lx)) {
      // center letter above/below axis end
      lx = x;
      ly = (y >= y0) ? y + LYSIZE+LOFF : y - LYSIZE-LOFF;
    }
    else if (lx) {
      // diagonal axis - center letter on axis
      lx = (x >= x0) ? x + LXSIZE+LOFF : x - LXSIZE-LOFF;
      ly = (y >= y0) ? y + LYSIZE+LOFF : y - LYSIZE-LOFF;
    }
    else {
      // axis is perp to screen - center letter at axis end
      lx = x;
      ly = y;
    }
    SelectObject( hdc, axis_pen[k] );

    // draw axis
    MoveToEx( hdc, x0, y0, NULL );
    LineTo( hdc, x, y );

    // draw axis label
    switch (k) {
    case 0: // X
      MoveToEx( hdc, lx-LXSIZE, ly-LYSIZE, NULL );
      LineTo(   hdc, lx+LXSIZE, ly+LYSIZE );
      MoveToEx( hdc, lx-LXSIZE, ly+LYSIZE, NULL );
      LineTo(   hdc, lx+LXSIZE, ly-LYSIZE );
      break;
    case 1: // Y
      MoveToEx( hdc, lx-LXSIZE, ly-LYSIZE, NULL );
      LineTo(   hdc, lx,    ly    );
      LineTo(   hdc, lx+LXSIZE, ly-LYSIZE );
      MoveToEx( hdc, lx,    ly, NULL    );
      LineTo(   hdc, lx,    ly+LYSIZE );
      break;
    case 2: // Z
      MoveToEx( hdc, lx-LXSIZE, ly-LYSIZE, NULL );
      LineTo(   hdc, lx+LXSIZE, ly-LYSIZE );
      LineTo(   hdc, lx-LXSIZE, ly+LYSIZE );
      LineTo(   hdc, lx+LXSIZE, ly+LYSIZE );
      break;
    }

  }
  SelectObject( hdc, saved_pen );
  DeleteObject( axis_pen[0] );
  DeleteObject( axis_pen[1] );
  DeleteObject( axis_pen[2] );

#undef LXSIZE
#undef LYSIZE
#undef LOFF

}
#endif

///////////////////////////////////////////////////////////////////////

void myBuildDisplayList( GLuint display_list_number,
                         GLUnurbsObj* pTheGLNurbsRender,
                         const CModel& model
                         )
{
  ON_Material material;
  glNewList( display_list_number, GL_COMPILE );

  // display Rhino geometry using ON_GL() functions found in rhinoio_gl.cpp
  int i;
  const int object_count = model.m_object_table.Count();
  for ( i = 0; i < object_count; i++ ) 
  {
    const ONX_Model_Object& mo = model.m_object_table[i];
    if ( 0 != mo.m_object )
    {
      model.GetObjectMaterial( i, material );
      myDisplayObject( *mo.m_object, material, pTheGLNurbsRender );
    }
  }

  glEndList();
}

///////////////////////////////////////////////////////////////////////

void MY_GL_CALLBACK myDisplay( void )
{
  // Uses globals glb_* because the GL aux tools don't provide an
  // easy way to pass information into this callback.
  int bUseRhinoSpotlights = false; // I like to use a simple headlight
                                   // for a basic preview.

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  // lights
  if ( bUseRhinoSpotlights && glb_model ) {
    // Rhino spotlights (currently rotate along with geometry)
    myDisplayLighting( glb_model->m_view.m_vp, *glb_model );
  }
  else {
    // simple bright white headlight
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    GLfloat pos[4]  = { (GLfloat)0.0, (GLfloat)0.0, (GLfloat)1.0, (GLfloat)0.0 };
    glLightfv( GL_LIGHT0, GL_POSITION,  pos );
    GLfloat black[4] = { (GLfloat)0.0, (GLfloat)0.0, (GLfloat)0.0, (GLfloat)1.0 };
    GLfloat white[4] = { (GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0, (GLfloat)1.0 };
    glLightfv( GL_LIGHT0, GL_AMBIENT,  black );
    glLightfv( GL_LIGHT0, GL_DIFFUSE,  white );
    glLightfv( GL_LIGHT0, GL_SPECULAR, white );
    glEnable( GL_LIGHT0 );
    glPopMatrix();
  }

  // display list built with myBuildDisplayList()
  glCallList( glb_display_list_number ); 

  glFlush();

#if defined(MY_USE_WINDOWS_STUFF)
  // Windows decorations
  myDrawAxesSprite( glb_model->m_view.m_vp, wglGetCurrentDC() );
#endif
}

///////////////////////////////////////////////////////////////////////

void MY_GL_CALLBACK myNurbsErrorCallback( GLenum errCode )
{
  const GLubyte* s = gluErrorString( errCode );
  printf("GL NURBS ERROR: (%d) %s\n",errCode, s );
}

///////////////////////////////////////////////////////////////////////
