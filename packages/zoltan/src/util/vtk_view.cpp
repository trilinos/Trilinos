// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

//--------------------------------------------------------------------------
// This source file builds two applications: vtk_view and vtk_write.
// vtk_write is built when OUTPUT_TO_FILE is defined.
//
// Both applications read a Chaco or Exodus file and zdrive output
// files, and render the mesh colored by the Zoltan partitioning
// computed by zdrive.
//
// vtk_view brings up a window on your display in which you can rotate,
// zoom in, etc. the mesh.
//
// vtk_write writes one or more images to a file which you may then
// incorporate into a web page or report or mail to a friend.  It uses
// Mesa to do off screen rendering so it can be run in places where
// vtk_view can not be run.
//                                                                          
// Both can run with more or fewer processes than zdrive ran.  Choose   
// the size of your vtk_view/vtk_write application based only on the 
// computational demands of reading and rendering the Chaco or Nemesis files.
// 
// The inputs to vtk_view/vtk_write are:
//                                                                          
//   The Chaco or Nemesis files used by the zdrive test.                    
//                                                                          
//   The zdrive output files describing the calculated partitioning.        
//                                                                          
//   An input file similar to the input file used by zdrive, naming the     
//     Chaco or Nemesis files.  In fact, you can use the zdrive parameter
//     file since it named those input files.  When visualizing parallel
//     ExodusII/Nemesis files, it helps if you add the parameter "Zdrive Count" 
//     to tell vtk_view how many zdrive processes there were.  If you omit 
//     this, vtk_view will try to figure it out by searching for the
//     zdrive output files.                
//                                                                          
// vtk_view requires that you obtain and compile VTK (the Visualization     
// Toolkit, available for free from Kitware, Inc. at www.vtk.org), version
// 5.0 or later.  In the  Config.{arch} file that directed your compilation 
// of Zoltan, add directives to indicate where the VTK, GL and X libraries 
// are.  See "Config.generic" for an example of these directives.
//
// vtk_write must be compiled with a version of VTK libraries that were
// built with the MANGLED_MESA directive.  vtk_write must be linked with
// the Mesa libraries, which replace most of the Open GL libraries.
//
// Usage:  vtk_view parameterfileName
// Usage:  vtk_write parameterfileName
//
// The default parameter file name is "zdrive.inp" in the current working
// directory.  Parameter file options are described in the "usage" functions
// in this file.
//--------------------------------------------------------------------------

// system headers

#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include <map>
#include <iostream>
#include <cmath>

using namespace std;

//#define DEBUG_PARTITION_IDS

# ifdef MAXPATHLEN
// BSD
#define ZMAXPATH MAXPATHLEN
# else
#ifdef PATH_MAX
// Posix
#define ZMAXPATH PATH_MAX
#else
// Rational guess
#define ZMAXPATH 5024
#endif  //PATH_MAX
#endif //MAXPATHLEN

// VTK headers

#ifdef OUTPUT_TO_FILE
#include "vtkWindowToImageFilter.h"
#include "vtkTIFFWriter.h"
#include "vtkBMPWriter.h"
#include "vtkPNGWriter.h"
#include "vtkJPEGWriter.h"
#include "vtkPostScriptWriter.h"
#include "vtkXMesaRenderWindow.h"
#include "vtkMesaRenderer.h"
#include "vtkMesaActor.h"
#include "vtkMesaPolyDataMapper.h"
#include "vtkPKdTree.h"
#else
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
#endif

#include "vtkLookupTable.h"
#include "vtkCompositeRenderManager.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDistributedDataFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkMPIController.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIGroup.h"
#include "vtkCamera.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkIntArray.h"
#include "vtkCell.h"
#include "vtkPExodusReader.h"
#include "vtkPChacoReader.h"

// Zoltan headers

extern "C" {
#define MPICH_SKIP_MPICXX
#include "dr_input_const.h"
#include "../ch/ch_init_dist_const.h"
}

// MPI environment

vtkMPICommunicator *Comm;
vtkMPIController *Controller;
static int NumProcs, Proc;

#ifdef OUTPUT_TO_FILE
#define TIFF_FORMAT  1
#define BMP_FORMAT   2
#define PNG_FORMAT   3
#define PS_FORMAT    4
#define JPEG_FORMAT  5

static char *suffix[6] = {
"unknown",
"tif",
"bmp",
"png",
"ps",
"jpg"
};
#endif

// parameters from zdrive-style input file

#define LINELEN 80
 
static char *vis_opt_names[UNDEFINED_LIST_MAX]={
  "zdrive count",     // number of zdrive processes 
  "image height",    // pixel height of image
  "image width",     // pixel width of image
  "output partition number", // visualize this partition number

                      // Options ONLY for interactive window display:
  "omit caption",     // don't display default caption in window
  "omit scalar bar",  // don't display scalar bar in window
  "add caption",      // add user's caption to window
                      // add new interactive options here <<<<<<<

                      // Options ONLY for output image file :
  "output format",    // tiff, png, etc.
  "output name",      // base name for output image file
  "output frame start",  // first image (default is 0)
  "output frame stop",   // last image (default is 0, 359 is complete circle)
  "output frame stride", // degrees from one image to the next
  "output view up",      // the "up" direction in the image
                      // add new output image options here <<<<<<<
};

enum {option_zdrive_count,
      option_height,
      option_width,
      option_partition_number,

      option_omit_caption,
      option_omit_scalar_bar,
      option_add_caption,
                      // add new interactive options here <<<<<<<
      option_format,
      option_name,
      option_frame_start,
      option_frame_stop,
      option_frame_stride,
      option_view_up,
                      // add new output image options here <<<<<<<
      option_end};

static PARIO_INFO fname_opts;
static PROB_INFO prob_opts;
static UNDEFINED_INFO extra_options;
static char vis_opt_values[UNDEFINED_LIST_MAX][UNDEFINED_LENGTH_MAX];

// default option values

static int zdriveCount = 0; 
static int imageHeight = 300;
static int imageWidth= 300;
static int showPartitionNumber = -1;
static int notPartitionNumber = -1;

#ifndef OUTPUT_TO_FILE
static int omitCaption = 0;
static int omitScalarBar = 0;
static char addCaption[1024] = {'\0'};
#else
static int outputFormat = TIFF_FORMAT;
static char outputName[ZMAXPATH] = "outfile";
static int outputStart = 0;
static int outputStop = 0;
static int outputStride = 1;
static float outputViewUp[3] = {0, 1, 0};
#endif

// Two ways to specify distributed files to parallel VTK ExodusII reader:

static char filePattern[ZMAXPATH]; // "printf"-like pattern: "%s/basename.0016.%04d"
static char filePrefix[ZMAXPATH];  //   root directory for pattern
static int fileRange[2];      //   file numbers for pattern
static char **fileNames;  // OR list of all file names 
static char **fileNamesBase;
static int specialSingleFile = 0;
static int numNemesisFiles = 0; 
static int lastPartition = 0; 

static int read_broadcast_input_options(int &argc, char **argv);
static int read_mesh(vtkUnstructuredGrid *ug);
static int create_field_array(vtkUnstructuredGrid *ug, char *ca, char *pa, 
  double *range);
static int set_number_of_zdrive_processes(char *fname, char **disks);
static int set_nemesis_file_names_or_pattern(char *baseName,
  int numDisks, int diskListSize, int *diskList, int diskOffset, int useZeros,
  char *dirRoot, char *subDir);
static int check_partition_numbers(vtkUnstructuredGrid *ug);
static int assign_partition_numbers(vtkUnstructuredGrid *ug);
static void Run(vtkMultiProcessController *c, void *arg);
static int checkAllrc(int rc, vtkMPICommunicator *comm);
static int check_valid_range(int val, char *nm, int lo, int hi);
static void remove_trailing_junk(char *cbegin);
static void usage();
#ifdef DEBUG_PARTITION_IDS
static void debug_partition_ids(vtkPoints *pts, vtkIntArray *partids);
static void debug_partition_ids(vtkUnstructuredGrid *ug, vtkIntArray *partids);
#endif

#ifdef OUTPUT_TO_FILE
static void view_option_ignored(int option);
#else
static void write_option_ignored(int option);
static char *get_zdrive_output_file_name();
static char *captionText(char *p, char *c);
#endif

static int *realloc(int *buf, int newsize, int oldsize)
{
  if (newsize <= 0) return NULL;
  int *newbuf = new int [newsize];
  if (newbuf == NULL) return NULL;

  if (buf && (oldsize > 0)){
    memcpy(newbuf, buf, oldsize * sizeof(int));
    delete [] buf;
  }

  return newbuf;
}

int main(int argc, char **argv)
{
  // Establish that we are an MPI application

  Controller = vtkMPIController::New();
  Controller->Initialize(&argc, &argv);

  vtkMPIController::SetGlobalController(Controller);

  NumProcs = Controller->GetNumberOfProcesses();
  Proc = Controller->GetLocalProcessId();

  Comm = vtkMPICommunicator::SafeDownCast(Controller->GetCommunicator());

  // Read and broadcast the input parameters

  int rc = read_broadcast_input_options(argc, argv);

  if (rc){
    Controller->Delete();
    return 1;
  }

  if (fname_opts.file_type == NEMESIS_FILE){

    // Figure out the names and directories for the distributed files

    rc = set_nemesis_file_names_or_pattern(fname_opts.pexo_fname,
      fname_opts.num_dsk_ctrlrs,
      fname_opts.dsk_list_cnt, fname_opts.dsk_list,
      fname_opts.pdsk_add_fact, fname_opts.zeros, 
      fname_opts.pdsk_root, fname_opts.pdsk_subdir);


    if (NumProcs > 1){
      rc = checkAllrc(rc, Comm);
    }

    if (rc > 0){
      Controller->Finalize(); 
      Controller->Delete();
    
      return 1;
    }
  }
  else if ((fname_opts.file_type == CHACO_FILE) ||
           (fname_opts.file_type == NO_FILE_TRIANGLES)){

    // The chaco base name is in fname_opts.pexo_fname, but if we
    // are to read the zdrive output files, we need to know how many
    // zdrive processes there were.
  
    if (zdriveCount == 0){
      
      if (Proc == 0){
        rc = set_number_of_zdrive_processes(fname_opts.pexo_fname, NULL);
      }
  
      if (NumProcs > 1){
        int vals[2];
        vals[0] = rc;
        vals[1] = zdriveCount;
        Comm->Broadcast(vals, 2, 0);
 
        if (Proc > 0){
          rc = vals[0];
          zdriveCount = vals[1];
        }
      }

      if (rc){
        return 1; 
      }
    }
  }
  else{
    if (Proc == 0){
      cout << "Please specify either \"File Type = NemesisI\"";
      cout << " or \"File Type = Chaco\"" << endl;
      cout << " or \"File Type = random-triangles \"" << endl;
      cout << "in parameter file" << endl;
    }
    return 1;
  }

  // Run parallel method, VTK-style.

  Controller->SetSingleMethod(Run, NULL);
  Controller->SingleMethodExecute();

  Controller->Finalize();
  Controller->Delete();

  return 0;
}

static void Run(vtkMultiProcessController *c, void *arg)
{
  // Read in the Chaco or distributed Nemesis (exodusII) files.
  // Include at least one field array in the event we will not
  // have partition numbers to display.

  vtkMPIController *contr = vtkMPIController::SafeDownCast(c);

  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();

  int rc = read_mesh(ug);

  if (rc){
    ug->Delete();
    return;
  }

  // If we were processing zdrive output files, we'll visualize the
  // partition number assigned by zdrive.  If not, and there are
  // point or cell arrays other than global ID arrays, we'll visualize 
  // one of those.  If all else fails, we'll visualize the global element ID.

  char cellArray[128], pointArray[128];
  double range[2];

  rc = create_field_array(ug, cellArray, pointArray, range);

  if (rc){
    ug->Delete();
    return;
  }

  // Redistribute the cells for more load balanced rendering.  This
  // also creates ghost cells on each process if required for 
  // correct rendering.  For a single process application with 
  // distributed Exodus files, it removes duplicate points from
  // the vtkUnstructuredGrid output.

  vtkDistributedDataFilter *dd = vtkDistributedDataFilter::New();

  dd->SetInput(ug);
  dd->SetController(contr);
  dd->SetGlobalElementIdArrayName("GlobalElementId");
  dd->SetGlobalNodeIdArrayName("GlobalNodeId");

  // Create mapper to render the surface

  vtkDataSetSurfaceFilter *dss = vtkDataSetSurfaceFilter::New();
  dss->SetInput(dd->GetOutput());

#ifdef OUTPUT_TO_FILE
  vtkMesaPolyDataMapper *mapper = vtkMesaPolyDataMapper::New();
#else
  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
#endif

  mapper->SetInput(dss->GetOutput());
  mapper->SetColorModeToMapScalars();

  if (pointArray[0]){
    mapper->SetScalarModeToUsePointFieldData();
    mapper->SelectColorArray(pointArray);
  }
  else{
    mapper->SetScalarModeToUseCellFieldData();
    mapper->SelectColorArray(cellArray);
  }

  mapper->SetScalarRange(range[0], range[1]);

#if 0

  // Problem - with the default lookup table it is difficult to
  // distinguish the different regions, even with as few as
  // six partitions.  Here we try creating our own color map.
  // I haven't found one yet that give a satisfactory contrast.

  vtkIdType nvals = (vtkIdType)range[1] - (vtkIdType)range[0] + 1;
  vtkLookupTable *lut = vtkLookupTable::New();
  lut->SetTableRange(range);
  lut->SetNumberOfTableValues(nvals);
  double b0 = .1;                             // blue
  double bdiff = .8 / (nvals-1);
  double g0 = .0;                             // green
  double gdiff = 0.0; //.2 / (nvals-1);
  double r0 = .9;                             // red
  double rdiff = -.8 / (nvals-1);
  for (vtkIdType i=0; i<nvals; i++){
    lut->SetTableValue(i, r0 + i*rdiff, g0 + i*gdiff, b0 + i*bdiff);
  }
 
  mapper->SetLookupTable(lut);
#endif

  // Create actors (I don't know how to render text or scalar bar to a file)

#ifdef OUTPUT_TO_FILE

  vtkMesaActor *actor = vtkMesaActor::New();
  actor->SetMapper(mapper);

#else
  vtkScalarBarActor *sb = NULL;
  vtkTextActor *capActor = NULL;
  
  if (Proc == 0){
    if (!omitScalarBar){
      sb = vtkScalarBarActor::New();
      sb->SetOrientationToVertical();
      if (pointArray[0]){
        sb->SetTitle(pointArray);
      }
      else{
        sb->SetTitle(cellArray);
      }
      sb->SetLookupTable(mapper->GetLookupTable());
      sb->SetNumberOfLabels(4);
    }
  
    if (!omitCaption || addCaption[0]){
      capActor = vtkTextActor::New();
      char *info = captionText(pointArray, cellArray);

      if (info){
        capActor->SetInput(info);
        delete [] info;
      
        capActor->SetAlignmentPoint(0);
        capActor->GetTextProperty()->SetVerticalJustificationToBottom();
        capActor->GetTextProperty()->SetJustificationToLeft();
        capActor->GetTextProperty()->SetFontSize(16);
        capActor->GetTextProperty()->BoldOn();
        capActor->GetTextProperty()->SetLineSpacing(1.2);
      }
      else{
        capActor->Delete();
      }
    }
  }
      
  vtkActor *actor = vtkActor::New();
  actor->SetMapper(mapper);

#endif

  // Update the mapper now, because vtkDistributedDataFilter is a
  // parallel filter and when rendering, node 0 tries to execute
  // it early on (if it's not updated) and the others do not.

  mapper->SetPiece(Proc);
  mapper->SetNumberOfPieces(NumProcs);
  mapper->Update();

  // Set up the renderer

#ifdef OUTPUT_TO_FILE

  vtkXMesaRenderWindow *renWin = vtkXMesaRenderWindow::New();
  renWin->OffScreenRenderingOn();

  vtkMesaRenderer *renderer = vtkMesaRenderer::New();
  renderer->AddActor(actor);
  renWin->AddRenderer(renderer);

#else

  vtkRenderer *renderer = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();

  renderer->AddActor(actor);
  renWin->AddRenderer(renderer);

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

  if (Proc == 0){
    if (sb){
      renderer->AddActor(sb);
      sb->Delete();
    }
    if (capActor){
      renderer->AddViewProp(capActor);
      capActor->Delete();
    }
  }
#endif

  renderer->SetBackground(0,0,0);
  renWin->SetSize(imageWidth,imageHeight);

  vtkCompositeRenderManager *prm = vtkCompositeRenderManager::New();
  prm->SetRenderWindow(renWin);
  prm->SetController(contr);
  prm->InitializePieces();

#ifdef OUTPUT_TO_FILE

  // Create the images and write them out.

  if (Proc == 0){
    char fname[128];
    int skipFrames = outputStride - 1;
    int i, ii; 

    vtkWindowToImageFilter *wif = vtkWindowToImageFilter::New();
    wif->SetInput(renWin);

    vtkImageWriter *writer = NULL;

    if (outputFormat == TIFF_FORMAT){
      writer = vtkTIFFWriter::New();
    }
    else if (outputFormat == BMP_FORMAT){
      writer = vtkBMPWriter::New();
    }
    else if (outputFormat == PNG_FORMAT){
      writer = vtkPNGWriter::New();
    }
    else if (outputFormat == JPEG_FORMAT){
      writer = vtkJPEGWriter::New();
    }
    else if (outputFormat == PS_FORMAT){
      writer = vtkPostScriptWriter::New();
    }

    writer->SetInput(wif->GetOutput());

    if (outputStop == outputStart){
      sprintf(fname, "%s.%s", outputName, vis_opt_values[option_format]);
    }           

    vtkCamera *camera = renderer->GetActiveCamera();
    camera->SetViewUp(outputViewUp[0], outputViewUp[1], outputViewUp[2]);
    camera->UpdateViewport(renderer);
  
    for (i=0, ii=skipFrames; i<=outputStop; i++) {
      if (i >= outputStart) {
        if (ii == skipFrames) {

          prm->ResetCamera(renderer);
          prm->ResetCameraClippingRange(renderer);
          renWin->Render();

          wif->Modified();
          if (outputStop > outputStart){
            sprintf(fname, "%s.%03d.%s", outputName, i, 
                    vis_opt_values[option_format]);
          }

          writer->SetFileName(fname);
          writer->Write();
          cout << "Wrote: " << fname << endl;

          ii = 0;
        }
        else {
          ii++;
        }
      }
      camera->Azimuth(1);  // rotate around mesh 1 degree
    }
    writer->Delete();
    wif->Delete();
    prm->StopServices();
  }
  else{
    prm->StartServices();
  }

#else

  renWin->SetPosition(0, 360*Proc); 
prm->ResetCamera(renderer);
prm->ResetCameraClippingRange(renderer);
  prm->StartInteractor(); // now you can interact with window
  iren->Delete();

#endif

  mapper->Delete();
  actor->Delete();
  renderer->Delete();
  renWin->Delete();

  dd->Delete();
  ug->Delete();
  dss->Delete();

  prm->Delete();
}
static int read_mesh(vtkUnstructuredGrid *ug)
{
  if ((fname_opts.file_type == CHACO_FILE) || (fname_opts.file_type == NO_FILE_TRIANGLES)){
    // Process 0 reads in chaco file, request VTK reader to add
    // element ID array.

    vtkPChacoReader *rdr = vtkPChacoReader::New();
    rdr->SetBaseName(fname_opts.pexo_fname);

    rdr->GenerateGlobalElementIdArrayOn();  
    rdr->GenerateGlobalNodeIdArrayOn();  

    rdr->GetOutput()->SetUpdatePiece(Proc);
    rdr->GetOutput()->SetUpdateNumberOfPieces(NumProcs);

    rdr->UpdateInformation();  // get the file metadata

    int nvweights = rdr->GetNumberOfVertexWeights();
    if (nvweights > 0){
      rdr->GenerateVertexWeightArraysOn();
    }
    else{
      int neweights = rdr->GetNumberOfEdgeWeights();
      if (neweights > 0){
        rdr->GenerateEdgeWeightArraysOn();
       }
    }

    rdr->GetOutput()->SetUpdatePiece(Proc);
    rdr->GetOutput()->SetUpdateNumberOfPieces(NumProcs);

    rdr->Update(); // read in the files
    ug->DeepCopy(rdr->GetOutput());
    rdr->Delete();
  }
  else if (fname_opts.file_type == NEMESIS_FILE){

    // This is a parallel VTK Exodus II reader.  It will
    // distribute the input files across the processes.

    vtkPExodusReader *rdr = vtkPExodusReader::New();

    if (fileNames){
      rdr->SetFileNames(numNemesisFiles, (const char **)fileNames);
    }
    else{
      rdr->SetFilePrefix(filePrefix); 
      rdr->SetFilePattern(filePattern);
      rdr->SetFileRange(fileRange);
    }

    rdr->GenerateGlobalElementIdArrayOn();  
    rdr->GenerateGlobalNodeIdArrayOn();  
    rdr->GetOutput()->SetUpdatePiece(Proc);
    rdr->GetOutput()->SetUpdateNumberOfPieces(NumProcs);
    rdr->SetTimeStep(0);

    rdr->UpdateInformation();   // get the metadata

    // These counts don't include GlobalNodeIdArray or GlobalElementIdArray

    int nparrays = rdr->GetNumberOfPointArrays();
    if (nparrays > 0){
      rdr->SetPointArrayStatus(0, 1); 
    }
    else{
      int ncarrays = rdr->GetNumberOfCellArrays();
      if (ncarrays > 0){
        rdr->SetCellArrayStatus(0, 1); 
      }
    }

    rdr->GetOutput()->SetUpdatePiece(Proc);
    rdr->GetOutput()->SetUpdateNumberOfPieces(NumProcs);

    rdr->Update(); // read in the files
    ug->DeepCopy(rdr->GetOutput());
    rdr->Delete();
  }

  if (ug->GetNumberOfCells() < NumProcs) {
    if (Proc == 0) {
      if (ug->GetNumberOfCells() > 0){
        cout << "Input dataset is too small, only ";
        cout << ug->GetNumberOfCells() << " cells." << endl;
      }
      else{
        cout << "No input dataset." << endl;
      }
    }
    return 1;
  }
  return 0;
}
static int create_field_array(vtkUnstructuredGrid *ug, 
  char *ca, char *pa, double *range)
{
  int rc = 0;

  ca[0] = pa[0] = '\0';

  if (zdriveCount > 0){
    // Find the partition number for each "cell" in my subgrid.  This
    // requires reading the zdrive output files.  We name the
    // resulting element or point array "Partition".

    if ((fname_opts.file_type == CHACO_FILE) || (fname_opts.file_type == NO_FILE_TRIANGLES)){
      strcpy(pa, "Partition");
    }
    else{
      strcpy(ca, "Partition");
    }

    rc = assign_partition_numbers(ug);
  
    if (NumProcs > 1){
      rc = checkAllrc(rc, Comm);
    }
  
    if (rc > 0){
      return 1;
    }
  
    rc = check_partition_numbers(ug);
  
    if (rc){
      cout << Proc << " failed to obtain all partition numbers" << endl;
    }
  
    if (NumProcs > 1){
      rc = checkAllrc(rc, Comm);
    }
  
    if (rc > 0){
      return 1;
    }

    range[0] = 0;
    range[1] = lastPartition;
  }
  else{
    int nparrays = ug->GetPointData()->GetNumberOfArrays();

    if (nparrays > 1){
      const char *nm = ug->GetPointData()->GetArrayName(0);
      if (!strcmp(nm, "GlobalNodeId")){
        nm = ug->GetPointData()->GetArrayName(1);
      }
      strcpy(pa, nm);

      ug->GetPointData()->GetArray(pa)->GetRange(range);
    }
    else{
      int ncarrays = ug->GetCellData()->GetNumberOfArrays();

      if (ncarrays > 1){
        const char *nm = ug->GetCellData()->GetArrayName(0);
        if (!strcmp(nm, "GlobalElementId")){
          nm = ug->GetCellData()->GetArrayName(1);
        }
        strcpy(ca, nm);
      }
      else{
        strcpy(ca, "GlobalElementId");
      }

      ug->GetCellData()->GetArray(ca)->GetRange(range);
    }

    if (NumProcs > 1){
      double extreme;
      Comm->ReduceMin(range, &extreme, 1, 0);
      Comm->Broadcast(&extreme, 1, 0);

      range[0] = extreme;

      Comm->ReduceMax(range + 1, &extreme, 1, 0);
      Comm->Broadcast(&extreme, 1, 0);

      range[1] = extreme;
    }
  }

  return 0;
}

//----------------------------------------------------------------
// Functions to read parameter file
//----------------------------------------------------------------

int read_broadcast_input_options(int &argc, char **argv)
{
  // Options (same meaning as zdrive parameters):
  //   "File Type"     NemesisI or Chaco
  //   "File Name"     base name of Nemesis (base.p.n) or 
  //                     Chaco (base.coords, base.graph) files
  //   "Parallel Disk Info"
  //   "Parallel file location" 
  //
  // Special vtk_view options:
  //    See the "usage" functions in this source file.

  const char  *cmd_file;

  switch(argc)
  {
  case 1:
    cmd_file = "zdrive.inp";
    break;

  case 2:
    cmd_file = argv[1];

    if ( !strncmp(cmd_file,"-h", 2) || !strncmp(cmd_file,"-H", 2) ||
         !strcmp(cmd_file,"help") || !strcmp(cmd_file,"HELP")){

      if (Proc == 0){
        usage();
      }
      return 1;
    }

    break;

  default:
    if (Proc == 0) {
      usage();
    }
    return 1;
  }
  fname_opts.dsk_list_cnt         = -1;
  fname_opts.num_dsk_ctrlrs       = -1; 
  fname_opts.pdsk_add_fact        = -1; 
  fname_opts.zeros                = -1;
  fname_opts.file_type            = -1; 
  fname_opts.pdsk_root[0]         = '\0'; 
  fname_opts.pdsk_subdir[0]       = '\0';
  fname_opts.pexo_fname[0]        = '\0';

  for (int i=0; i<option_end; i++){
    vis_opt_values[i][0] = '\0'; // default is "not set"
  }

  int input_ok = 1;

  if (Proc == 0){
    input_ok = read_cmd_file(cmd_file, &prob_opts, &fname_opts, &extra_options);
  
    if (input_ok){
      if (extra_options.list_size > 0){
        for (int i=0; input_ok && (i<extra_options.list_size); i++){
          input_ok = 0;
          for (int j=0; j<option_end; j++){
            if (strstr(extra_options.line[i], vis_opt_names[j])){
              char *c = strchr(extra_options.line[i], '=');
              if (c){
                c++;
                while (*c && !isgraph(*c)) c++;
                if (*c){
                  strncpy(vis_opt_values[j], c, UNDEFINED_LENGTH_MAX);
                  input_ok = 1;
                }
              }
            }
          }
        } 
      }
    }

    if (fname_opts.file_type == HYPERGRAPH_FILE){
      // We don't have a way to visualize the hyperedges at this point
      // in time.  But some hypergraph files may have an associated
      // Chaco mesh that we can display, along with the partitioning
      // of vertices.
      //
      // Ideas for visualizing hypergraph partitionings:
      // We could read in the hypergraph file, create a vtkPolyData
      // containing the edges, and then use vtkGraphLayoutFilter to
      // get a vtkPolyData that we can render.  Another nice thing
      // would be a call back where every hit of a key causes the
      // next hyperedge to by highlighted.  The color of the vertices
      // would represent the partition.  For small hypergraphs, we
      // would be able to step through and see how the hyperedges are 
      // cut.  More difficult would be offsetting the edges from the 
      // vertices so we can see the separate overlapping hyperedges.
  
      fname_opts.file_type = CHACO_FILE;
  
      if (Proc == 0){
        cout << 
          "Warning: Hyperedge visualization not implemented yet." << endl;
        cout << 
          "We'll search for a Chaco file and display that if found." << endl;
      }
    }
  
    if (input_ok && 
        (fname_opts.file_type != NO_FILE_TRIANGLES ) &&  /* zdrive generated a chaco file */
        (fname_opts.file_type != CHACO_FILE) &&
        (fname_opts.file_type != NEMESIS_FILE)){
      cout << "vtk_view can only read chaco or exodusII/Nemesis files." << endl;
      cout << "Please specify either \"File Type = NemesisI\" or " << endl;
      cout << "\"File Type = Chaco\" in " << cmd_file << "." << endl;
      input_ok = 0;
    } 
  }

  if (NumProcs > 1){
    Comm->Broadcast(&input_ok, 1, 0);
  }
  
  if (!input_ok){
    if (Proc == 0){
      usage();
    }
    return 1;
  }
  
  if (NumProcs > 1){

    Comm->Broadcast((char *)&fname_opts, sizeof(fname_opts), 0);
  
    if (fname_opts.dsk_list_cnt > 0){
      if (Proc != 0){
        fname_opts.dsk_list = new int[fname_opts.dsk_list_cnt];
      }
      Comm->Broadcast(fname_opts.dsk_list, fname_opts.dsk_list_cnt, 0);
    }

    int *set = new int [option_end];

    if (Proc == 0){
      for (int i=0; i<option_end; i++){
        set[i] = (extra_options.line[i][0] != '\0');
      }
    }
    Comm->Broadcast(set, option_end, 0);

    for (int i=0; i<option_end; i++){
      if (!set[i]) continue;
      Comm->Broadcast(extra_options.line[i], UNDEFINED_LENGTH_MAX, 0);
    }
    delete [] set;
  }

  if (vis_opt_values[option_zdrive_count][0]){
    sscanf(vis_opt_values[option_zdrive_count], "%d", &zdriveCount);
    numNemesisFiles = zdriveCount;
  }
  if (vis_opt_values[option_height][0]){
    sscanf(vis_opt_values[option_height], "%d", &imageHeight);
  }
  if (vis_opt_values[option_width][0]){
    sscanf(vis_opt_values[option_width], "%d", &imageWidth);
  }
  if (vis_opt_values[option_partition_number][0]){
    sscanf(vis_opt_values[option_partition_number], "%d", 
      &showPartitionNumber);
    notPartitionNumber = (showPartitionNumber ? 0 : 1);
  }

  int rc = check_valid_range(imageHeight, vis_opt_names[option_height], 30, 2500);

  if (rc == 0){
    rc = check_valid_range(imageWidth, vis_opt_names[option_width], 30, 2500);
  }

  if (rc){
    return 1;
  }

#ifdef OUTPUT_TO_FILE
  if (vis_opt_values[option_format][0]){
    char *c = vis_opt_values[option_format];
    remove_trailing_junk(c);

    if (!strcmp(c, "tiff") || !strcmp(c, "tif")){
      outputFormat = TIFF_FORMAT;
    }
    else if (!strcmp(c, "bmp")){
      outputFormat = BMP_FORMAT;
    }
    else if (!strcmp(c, "png")){
      outputFormat = PNG_FORMAT;
    }
    else if (!strcmp(c, "ps")){
      outputFormat = PS_FORMAT;
    }
    else if (!strcmp(c, "jpeg") || !strcmp(c, "jpg")){
      outputFormat = JPEG_FORMAT;
    }
    else{
      if (Proc == 0) {
        cerr << "Unrecognized output format " << c << endl;
        cerr << "Choices are tiff, bmp, png, ps or jpeg" << endl;
      }

      return 1;
    } 
  }

  strcpy(vis_opt_values[option_format], suffix[outputFormat]);

  if (vis_opt_values[option_name][0]){
    remove_trailing_junk(vis_opt_values[option_name]);
    snprintf(outputName, ZMAXPATH-1, "%s", vis_opt_values[option_name]);
  }
  if (vis_opt_values[option_frame_start][0]){
    sscanf(vis_opt_values[option_frame_start], "%d", &outputStart);
  }
  if (vis_opt_values[option_frame_stop][0]){
    sscanf(vis_opt_values[option_frame_stop], "%d", &outputStop);
  }
  if (vis_opt_values[option_frame_stride][0]){
    sscanf(vis_opt_values[option_frame_stride], "%d", &outputStride);
  }
  if (vis_opt_values[option_view_up][0]){
    int nvals = sscanf(vis_opt_values[option_view_up], "%f %f %f", 
      outputViewUp, outputViewUp + 1, outputViewUp + 2);

    if (nvals != 3){
      if (Proc == 0){
        cerr << "\"" << vis_opt_names[option_view_up] << "\"" ;
        cout << " is a 3-d vector like \"0 1 1\"\n";
        cerr << "It points in the \"up\" direction of the output image." << endl;
      }
      return 1;
    }
  }

  // ignore options that are only used when displaying a window to the screen
  for (int i = option_omit_caption; i < option_format; i++){
    if (vis_opt_values[i][0]){
      view_option_ignored(i);
    }
  }

  // Sanity checks

  if ((outputStart < 0) || (outputStop > 360) || (outputStart > outputStop)){
    if (Proc == 0){
      cerr << "Error: \"" << vis_opt_names[option_frame_start] << "\"";
      cerr << " and/or \"" << vis_opt_names[option_frame_stop] << "\"" << endl;
      cerr << "Valid range is 0 through 360" << endl;
    }
    return 1;
  }
  if (outputStride < 1){
    outputStride = 1;
  }


#else
  if (vis_opt_values[option_omit_caption][0]){
    sscanf(vis_opt_values[option_omit_caption], "%d", &omitCaption);
  }
  if (vis_opt_values[option_omit_scalar_bar][0]){
    sscanf(vis_opt_values[option_omit_scalar_bar], "%d", &omitScalarBar);
  }
  if (vis_opt_values[option_add_caption][0]){
    remove_trailing_junk(vis_opt_values[option_add_caption]);
    strncpy(addCaption, vis_opt_values[option_add_caption], 1023);
  }
  // ignore options that are only used when writing to an image file
  for (int i = option_format; i < option_end; i++){
    if (vis_opt_values[i][0]){
      write_option_ignored(i);
    }
  }
#endif

  return 0;
}
static void remove_trailing_junk(char *cbegin)
{
  char *c = cbegin + strlen(cbegin) - 1;
  while (!isprint(*c)){
     *c-- = '\0';
  }
}
static int check_valid_range(int val, char *nm, int lo, int hi)
{
  if ((val < lo) || (val > hi)){
    if (Proc == 0){
      cerr << "Error: valid range for " << nm << " is ";
      cerr << lo << " through " << hi << endl;
    }
    return 1;
  }

  return 0;
}
#ifdef OUTPUT_TO_FILE
static void view_option_ignored(int option)
{
  if (Proc == 0){
    cout << "Warning: " << vis_opt_names[option];
    cout << " ignored.  It is only for use with vtk_view." << endl;
  }
}
#else
static void write_option_ignored(int option)
{
  if (Proc == 0){
    cout << "Warning: " << vis_opt_names[option];
    cout << " ignored.  It is only for use with vtk_write." << endl;
  }
}
#endif

static int set_nemesis_file_names_or_pattern(char *baseName,
  int numDisks, int diskListSize, int *diskList, 
  int diskOffset, int useZeros,
  char *dirRoot, char *subDir)
{
  char **disks=NULL;
  int nzeros, num;
  int retVal = 0;
  int len = 0;
  int rc = 0;

  fileNames = NULL;
  fileNamesBase = NULL;
  filePattern[0] = '\0';
  filePrefix[0] = '\0';
  fileRange[0] = fileRange[1] = -1;

  if ((numDisks > 0) || (diskListSize > 0)){  // get list of directory names 

    int ndisks = (numDisks > 0) ? numDisks : diskListSize;

    if (useZeros){
      nzeros = ((ndisks < 9) ? 2 : 1+(int)(std::log10((float)ndisks)));
    }
    else{
      nzeros = 0;
    }
    
    len = strlen(dirRoot) + 32;
    if (subDir)  len += strlen(subDir);

    disks = new char *[ndisks];

    for (int i=0; i<ndisks; i++){

      disks[i] = new char[len];
      int dnum = (diskList ? diskList[i] : i+ diskOffset);

      sprintf(disks[i], "%s%0*d/%s", dirRoot, nzeros, dnum, subDir);
    }
  }

  if (zdriveCount == 0){
    
    if (Proc == 0){
      // also sets numNemesisFiles, which is usually the same
      rc = set_number_of_zdrive_processes(baseName, disks);
    }

    if (NumProcs > 1){
      int counts[4];
      counts[0] = rc;
      counts[1] = zdriveCount;
      counts[2] = numNemesisFiles;
      counts[3] = specialSingleFile;
      Comm->Broadcast(counts, 4, 0);

      if (Proc > 0){
        rc = counts[0];
        zdriveCount = counts[1];
        numNemesisFiles = counts[2];
        specialSingleFile = counts[3];
      } 
    }

    if (rc){
      retVal = 1;
      goto done2;
    }
  }

  num = 1+(int)(std::log10((float)numNemesisFiles));

  if (specialSingleFile){
    fileNames = new char *[1];
    fileNamesBase = new char *[1];
    fileNames[0] = strdup(baseName);
    fileNamesBase[0] = strdup(baseName);
  }
  else if (numDisks > 1){

    // We need to provide the VTK reader with a list of all file names.
    
    fileNames = new char *[numNemesisFiles];
    fileNamesBase = new char *[numNemesisFiles];
    len += (sizeof(baseName) + 64);

    for (int i=0; i<numNemesisFiles; i++){
      fileNames[i] = new char[len];
      fileNamesBase[i] = new char[len];
      sprintf(fileNamesBase[i],"%s/%s", disks[i%numDisks], baseName);
      sprintf(fileNames[i],"%s.%d.%0*d", fileNamesBase[i], numNemesisFiles, num, i);
    }
  }
  else{

    // Sufficient to provide just the file name prefix, pattern and range.
    
    fileRange[0] = 0;
    fileRange[1] = numNemesisFiles - 1;

    fileNamesBase = new char *[1];

    if (numDisks == 1){
      strcpy(filePrefix, disks[0]);
      sprintf(filePattern, "%%s/%s.%d.%%0%dd", baseName, numNemesisFiles, num);
      fileNamesBase[0] = new char [ZMAXPATH];
      sprintf(fileNamesBase[0],"%s/%s", disks[0], baseName);
    }
    else{
      strcpy(filePrefix, baseName);
      sprintf(filePattern, "%%s.%d.%%0%dd", numNemesisFiles, num);
      fileNamesBase[0] = strdup(baseName);
    }
  }

done2:

  if (disks){
    for (int i=0; i<numDisks; i++){
      delete [] disks[i];
    }
    delete [] disks;
  }

  return retVal;
}

// Look for the zdrive output files.  The file name tells us how many
// zdrive output files there are.  If the input file is Exodus/Nemesis,
// this number is also the number of Exodus/Nemesis files.
//
// If there are no zdrive output files, we'll just display the input
// file geometry without partition numbers.  If the input file is
// Exodus/Nemesis however, we still need to find out how many input
// files there are.

static int set_number_of_zdrive_processes(char *fname, char **disks)
{
  char *dn = new char[ZMAXPATH];
  char *fn = new char[ZMAXPATH];
  int verbose = (Proc == 0);

  if (disks){
    strcpy(dn, disks[0]);
    strcpy(fn, fname);
  }
  else{
    int slash = (int)('/');
    strcpy(dn, fname);

    int len = strlen(dn)-1; // remove slashes at the end of file name

    while ((len >= 0) && (dn[len] == slash)){
      dn[len] = '\0';
      len--;
    }

    if (len < 0){
      if (verbose){
        cout << "invalid file name " << fname << endl;
      }
      delete [] dn;
      delete [] fn;
      return 1;
    }

    // find the rightmost slash, which separates directory from file name
    
    char *c = strrchr(dn, slash);
    if (c==NULL){
      strcpy(dn, ".");
      strcpy(fn, fname);
    }
    else{
      *c = '\0';
      strcpy(fn, c+1);
    }
  }
  int len = strlen(fn);

  DIR *dir = opendir(dn);

  if (dir == NULL){
    if (verbose){
      cout << "zdrive output file directory: " << dn;
      cout << " " << strerror(errno) << endl;
    }
    delete [] dn;
    delete [] fn;
    return 1;
  }

  int fileno;
  int count = 0;
  struct dirent *de;

  while ((de = readdir(dir)) != NULL){
    if (strncmp(fn, de->d_name, len) == 0){
      char *c = de->d_name + len;

      if (*c == '.'){
        int nmatches = sscanf(c, ".out.%d.%d", &count, &fileno);

        if (nmatches == 2){
          break;
        }
      }
    }
  }

  if (count > 0){
    zdriveCount = count;
    numNemesisFiles = count;
  }
  else if (fname_opts.file_type == NEMESIS_FILE){
    // we must have a count for the number of input files
    rewinddir(dir);
    while ((de = readdir(dir)) != NULL){
      if (strncmp(fn, de->d_name, len) == 0){
        char *c = de->d_name + len;
  
        if (*c == '.'){
          int nmatches = sscanf(c, ".%d.%d", &count, &fileno);
  
          if (nmatches == 2){
            break;
          }
        }
        else if (!*c){
          /* 
           * The name given was not the base name for a distributed
           * file, it was the whole name of a single file.
           */
          count = 1;
          specialSingleFile = 1;
          break;
        }
      }
    }

    zdriveCount = 0;
    numNemesisFiles = count;
  }

  closedir(dir);

  int fail = ((fname_opts.file_type == NEMESIS_FILE) && (numNemesisFiles == 0));

  if (verbose){
    if (fail){
      cerr << "Unable to locate specified nemesis file(s)"<< endl;
    }
    else if (zdriveCount == 0){
      cerr << "Unable to locate zdrive output file(s) ";
      cerr << dn << "/" << fn << ".out.{numfiles}.{fileno}" << endl;
      cerr << "We'll just show you the geometry of the input files." << endl;
    }
  }

  delete [] dn;
  delete [] fn;

  return fail;
}
//----------------------------------------------------------------
// Functions to read the zdrive output files, and to create an
// element array in the distributed unstructured grid containing
// the partition number assigned to each element by Zoltan.
//----------------------------------------------------------------

static map<int, int> partitionMap;   // for searching
static int *gids=NULL;               // for message passing
static int *pids=NULL;
static int listSize=0;
static int listMax=0;
static int increment=10000;

static void clear_partition_ids()
{
  partitionMap.erase(partitionMap.begin(), partitionMap.end());
  if (gids) delete [] gids;
  if (pids) delete [] pids;
}
static void update_partition_id_map()
{
  partitionMap.erase(partitionMap.begin(), partitionMap.end());

  for (int i=0; i<listSize; i++){
    partitionMap.insert(pair<int, int>(gids[i], pids[i]));

    if (pids[i] > lastPartition){
      lastPartition = pids[i];
    }
  }
}
static void update_partition_ids(int *eltIds, vtkIntArray *partids)
{
  int *ids = partids->GetPointer(0);

  for (int i=0; i<partids->GetNumberOfTuples(); i++){

    if (ids[i] < 0){
      map<int,int>::iterator it = partitionMap.find(eltIds[i]);

      if (it->first == eltIds[i]){

        if (showPartitionNumber >= 0){
          if (it->second == showPartitionNumber){
            ids[i] = showPartitionNumber;
          }
          else{
            ids[i] = notPartitionNumber;
          }
        }
        else{
          ids[i] = it->second;
        }
      }
    }
  }
}
#ifdef DEBUG_PARTITION_IDS
static void debug_partition_ids(vtkUnstructuredGrid *ug, vtkIntArray *partids)
{
  double pcoords[3], weights[64], center[3];

  for (vtkIdType i=0; i<ug->GetNumberOfCells(); i++){
    vtkCell *c = ug->GetCell(i);
    int subId = c->GetParametricCenter(pcoords);
    c->EvaluateLocation(subId, pcoords, center, weights);
    
    cout << center[0] << " " << center[1] << " " << center[2] << ", ";
    cout << partids->GetTuple(i) << endl;
  }
}
static void debug_partition_ids(vtkPoints *pts, vtkIntArray *partids)
{
  for (vtkIdType i=0; i<pts->GetNumberOfPoints(); i++){
    double *p = pts->GetPoint(i);
    cout << p[0] << " " << p[1] << " " << p[2] << ", ";
    cout << partids->GetValue(i) << endl;
  }
}
#endif
static int read_zdrive_output(int from, int to)
{
char nm[ZMAXPATH], line[1024];
char *c = NULL;
FILE *fp;
int g, p, nvals;

  listSize = 0;

  int num = 1+(int)(std::log10((float)zdriveCount));

  if (!fileNames){
    if (fileNamesBase){
      c = fileNamesBase[0];
    }
    else{
      c = fname_opts.pexo_fname;
    }
  }

  for (int i=from; i<=to; i++){

    if (fileNames){
      c = fileNamesBase[i];
    }

    sprintf(nm, "%s.out.%d.%0*d",c,zdriveCount,num,i);

    fp = fopen(nm, "r");

    if (!fp){
      cout << Proc << " unable to open " << nm << endl;
      return 1;
    }

    while (fgets(line, 1024, fp) != NULL){

      if (!isdigit(line[0])) continue;
      nvals = sscanf(line, "%d %d", &g, &p);

      if (nvals != 2){
        fclose(fp);
        cout << Proc << " bad file format " << nm << endl;
        return 1;
      }

      if (listSize == listMax){
        gids = realloc(gids, listMax + increment, listMax);
        pids = realloc(pids, listMax + increment, listMax);
        listMax += increment;
      }

      gids[listSize] = g;
      pids[listSize] = p;

      listSize++;
    }
    fclose(fp);
  }
 
  update_partition_id_map();

  return 0;
}
static int pass_along_zdrive_output(int sndr, int recvr)
{
int msgsize;
vtkMPICommunicator::Request request;
int *buf;

  Comm->NoBlockReceive(&msgsize, 1, sndr, 0x1111, request);
  Comm->Send(&listSize, 1, recvr, 0x1111);

  request.Wait();

  buf = new int[msgsize];

  if (msgsize > listMax){
    gids = realloc(gids, msgsize, listMax);
    pids = realloc(pids, msgsize, listMax);
    listMax = msgsize;
  }

  Comm->NoBlockReceive(buf, msgsize, sndr, 0x1112, request);
  Comm->Send(gids, listSize, recvr, 0x1112);

  request.Wait();

  memcpy(gids, buf, sizeof(int) * msgsize);

  Comm->NoBlockReceive(buf, msgsize, sndr, 0x1113, request);
  Comm->Send(pids, listSize, recvr, 0x1113);

  request.Wait();

  memcpy(pids, buf, sizeof(int) * msgsize);

  listSize = msgsize;

  delete [] buf;

  update_partition_id_map();

  return 0; // how do we detect an error here?
}
static int check_partition_numbers(vtkUnstructuredGrid *ug)
{
  vtkIntArray *ids = NULL;

  if ((fname_opts.file_type == CHACO_FILE) || (fname_opts.file_type == NO_FILE_TRIANGLES)){
    ids = vtkIntArray::SafeDownCast(ug->GetPointData()->GetArray("Partition"));
  }
  else{
    ids = vtkIntArray::SafeDownCast(ug->GetCellData()->GetArray("Partition"));
  }

  if (!ids){
    return 1;
  }

  int *vals = ids->GetPointer(0);

  for (int i=0; i<ids->GetNumberOfTuples(); i++){
    if (vals[i] < 0){
      return 1;
    }
  }
  return 0;
}

int assign_partition_numbers(vtkUnstructuredGrid *ug)
{
  int nprocs, myrank, me, rc, pointArray;
  int from, to;
  int sndr, recvr;
  int nfiles, extra, noextra;
  int *procs, *hasCells;

  if ((fname_opts.file_type == CHACO_FILE) || (fname_opts.file_type == NO_FILE_TRIANGLES)){
    // zdrive treats the Chaco file vertices as the "elements" to be partitioned.
    pointArray = 1;
  }
  else{
    pointArray = 0;
  }

  vtkMPICommunicator *subComm = vtkMPICommunicator::New();
  vtkMPIGroup *group = vtkMPIGroup::New();

  int ncells = 0;

  if (pointArray){
    ncells = ug->GetNumberOfPoints();
  }
  else{
    ncells = ug->GetNumberOfCells();
  }

  // Add a new element array to the unstructured grid.  This
  // will contain the element's partition number.

  vtkIntArray *partids = vtkIntArray::New();
  partids->SetNumberOfValues(ncells);
  partids->SetName("Partition");

  for (int i=0; i<ncells; i++){
    partids->SetValue(i, -1);
  }

  // Get the global element IDs (Exodus) or point IDs (CHACO)

  vtkIntArray *ia = NULL;

  if (pointArray){
    ia = vtkIntArray::SafeDownCast(ug->GetPointData()->GetArray("GlobalNodeId"));
  }
  else{
    ia = vtkIntArray::SafeDownCast(ug->GetCellData()->GetArray("GlobalElementId"));
  }

  int *eltIds = ia->GetPointer(0);

  if (NumProcs == 1){
    nprocs = 1;
    myrank = 0;
  }
  else{
    // It's possible that the VTK reader output empty grids to
    // some of the processes.  They will not participate.

    hasCells = new int[NumProcs];
    procs    = new int[NumProcs];
    me = ((ncells > 0) ? 1 : 0);
  
    Comm->AllGather(&me, hasCells, 1);
  
    nprocs = 0;
    myrank = -1;
    group->Initialize(Controller);
  
    for (int i=0; i<NumProcs; i++){
      if (hasCells[i]){
        procs[nprocs] = i;
        group->AddProcessId(i);
        if (i == Proc){
          myrank = nprocs;
        }
        nprocs++;
      }
    }

    subComm->Initialize(Comm, group);

    if (myrank < 0){
      goto done;
    }
  }

  if (nprocs == 1){
    rc = read_zdrive_output(0, zdriveCount-1);
    if (rc){
      partids->Delete();
      partids = NULL;
    }
    else{
      update_partition_ids(eltIds, partids);
    }

    goto done;
  }

  // Divide the zdrive output files between the participating
  // processes.  Each participating process reads in it's share
  // of the files.

  nfiles = zdriveCount / nprocs;
  extra = zdriveCount - (nfiles * nprocs);
  noextra = nprocs - extra;

  if (myrank < noextra){
    from = nfiles * myrank;
    to = from  + nfiles - 1;
  }
  else{
    from = (nfiles * noextra) + ((myrank - noextra) * (nfiles + 1));
    to = from + nfiles;
  }


  rc = read_zdrive_output(from, to);

  if (NumProcs > 1){
    rc = checkAllrc(rc, subComm);
  }

  if (rc > 0){
    partids->Delete();
    partids = NULL;
    goto done;
  }

  update_partition_ids(eltIds, partids);

  sndr = (myrank + nprocs - 1) % nprocs;
  recvr = (myrank + 1) % nprocs;

  for (int i=0; i<nprocs-1; i++){
    rc = pass_along_zdrive_output(procs[sndr], procs[recvr]);

    rc = checkAllrc(rc, subComm);
    
    if (rc > 0){
      partids->Delete();
      partids = NULL;
      goto done;
    }
    update_partition_ids(eltIds, partids);
  }

done:

#ifdef DEBUG_PARTITION_IDS
  if (partids){
    for (int i=0; i<nprocs; i++){
      if (i == myrank){
        if (pointArray){
          cout << "Rank " << i << " partition IDs for points " << endl;
          debug_partition_ids(ug->GetPoints(), partids);
        }
        else{
          cout << "Rank " << i << " partition IDs for cell centers " << endl;
          debug_partition_ids(ug, partids);
        }
      }
      if ((nprocs > 1) && (myrank >= 0)){
        checkAllrc(1, subComm);   // barrier
      }
    }
  }
#endif

  subComm->Delete();
  group->Delete();

  clear_partition_ids();

  if (partids){
    if (pointArray){
      ug->GetPointData()->AddArray(partids);
    }
    else{
      ug->GetCellData()->AddArray(partids);
    }
    partids->Delete();
    return 0;
  }

  return 1;
}
static int checkAllrc(int rc, vtkMPICommunicator *comm)
{
  int allrc;
  comm->ReduceSum(&rc, &allrc, 1, 0);
  comm->Broadcast(&allrc, 1, 0);

  return allrc;
}
#ifndef OUTPUT_TO_FILE
static char *captionText(char *pnm, char *cnm)
{
  if (omitCaption && (addCaption[0] == '\0')){
    return NULL;
  }

  char *buf1 = new char [LINELEN  + 1];
  char *buf = new char [LINELEN * 3 + 1];
  char *c = buf;
  int lenleft = LINELEN * 3;
  int linelen = 0, used = 0;

  if (addCaption[0]){
    strncpy(buf, addCaption, lenleft);
    used = strlen(buf);
    lenleft -= used;
    c += used;
  }

  if (omitCaption){
    return buf;
  }

  used = snprintf(buf1, LINELEN,
     "(key \"t\" toggles interaction mode, \"r\" resets, \"q\" quits)\n");

  if (used <= lenleft){
    strcpy(c, buf1);
    lenleft -= used;
    c += used;
  }

  if (zdriveCount == 0){
    used = snprintf(buf1, LINELEN, "%s: %s\n",
      (fname_opts.file_type != NEMESIS_FILE ? "Chaco" : "Exodus/Nemesis"),
      fname_opts.pexo_fname);

    if (used <= lenleft){
      strcpy(c, buf1);
      lenleft -= used;
      c += used;

      used = snprintf(buf1, LINELEN, "displaying %s array: %s\n",
        pnm[0] ? "point" : "cell",
        pnm[0] ? pnm : cnm);

      if (used <= lenleft){
        strcpy(c, buf1);
        lenleft -= used;
        c += used;
      }
    }
  }
  else {

    // We had zdrive output files.  Display pertinent info about zdrive run

    int dtype = fname_opts.init_dist_type;
    if ((fname_opts.file_type != NEMESIS_FILE) && (dtype != INITIAL_FILE)) {
      snprintf(buf1, LINELEN, "Chaco: %s, %s initial distribution, %s",
        fname_opts.pexo_fname, 
        (dtype == INITIAL_LINEAR ? 
           "linear" : 
           (dtype == INITIAL_CYCLIC ? "cyclic" : "owner")),
        prob_opts.method);
    }
    else{
      snprintf(buf1, LINELEN, "%s: %s, %s",
        (fname_opts.file_type != NEMESIS_FILE ? "Chaco" : "Exodus/Nemesis"),
        fname_opts.pexo_fname, 
        prob_opts.method);
    }
  
    used = snprintf(c, LINELEN+1, "%s\n", buf1);
  
    c += used;
    lenleft -= used;
    linelen = 0;
  
    if (prob_opts.num_params > 0) {
      for (int i=0; i<prob_opts.num_params; i++) {
        Parameter_Pair p = prob_opts.params[i];
        used = snprintf(buf1, LINELEN, "%s(%s) ", p.Name, p.Val);
  
        if (used >= lenleft){  
          snprintf(c, lenleft, "...");
          break;  // that's all we have room for
        }
  
        if (linelen + used >= LINELEN){
          *c++ = '\n';
          linelen = 0;
          lenleft--;
        }
     
        strcpy(c, buf1);
        c += used;
        lenleft -= used;
        linelen += used;
      }
    }

    char *nm = get_zdrive_output_file_name();
    used = strlen(nm) + 1;
  
    if (lenleft > used){
      sprintf(c, "\n%s", nm);
      lenleft -= used;
      c += used;
    }
  
    delete [] nm;
  }

  return buf;
}
static char *get_zdrive_output_file_name()
{
  char *nm = new char [ZMAXPATH];

  sprintf(nm, "%s.out.%d.{n}",
    fileNames ? fileNamesBase[0] : fname_opts.pexo_fname, zdriveCount);

  return nm;
}
#endif

static void all_usage()
{
  cout << "zdrive count =" << endl;
  cout << "  The number of zdrive processes.  If you leave it blank," << endl;
  cout << "  we'll try to determine it by looking for zdrive output files." << endl;
  cout << endl;

  cout << "image height =" << endl;
  cout << "  The pixel height of the produced image." << endl;
  cout << "  Default is 300." << endl;
  cout << endl;

  cout << "image width =" << endl;
  cout << "  The pixel width of the produced image." << endl;
  cout << "  Default is 300." << endl;
  cout << endl;
}
#ifndef OUTPUT_TO_FILE
static void view_usage()
{
  cout << "omit caption = 1" << endl;
  cout << "  Don't display the default caption." << endl;
  cout << endl;

  cout << "omit scalar bar = 1" << endl;
  cout << "  Don't display the scalar bar." << endl;
  cout << endl;

  cout << "add caption = my caption text" << endl;
  cout << "  Add the specified caption to the image." << endl;
  cout << endl;
}
#else
static void write_usage()
{
  cout << "output format = " << endl;
  cout << "  The format for the image file (tiff, png, jpeg, ps, bmp)." << endl;
  cout << "  The default format is TIFF." << endl;
  cout << endl;

  cout << "output name = " << endl;
  cout << "  The basename for the output image file." << endl;
  cout << "  The default basename is \"outfile\"." << endl;
  cout << endl;

  cout << "output frame start = " << endl;
  cout << "  The camera is set up pointing toward the center of the mesh," << endl;
  cout << "  aligned with the Z-axis,  pointing to the negative Z-axis.  The \"up\"" << endl;
  cout << "  direction is the positive Y-axis.  (We use a right hand coordinate" << endl;
  cout << "  system.)  That is frame 0.  The camera can rotate around the " << endl;
  cout << "  mesh one degree at a time.  So frame numbers range from 0 to 359." << endl;
  cout << "  This is the number of the first frame.  The default is 0." << endl;
  cout << endl;

  cout << "output frame stop = " << endl;
  cout << "  This is the number of the last frame.  The default is 0." << endl;
  cout << "  (By default, we only write out one image.)" << endl;
  cout << endl;

  cout << "output frame stride = " << endl;
  cout << "  This the difference in degrees from one frame to the next." << endl;
  cout << "  The default is 1 degree." << endl;
  cout << endl;

  cout << "output view up = x y z" << endl;
  cout << "  These three numbers represent the direction that is \"up\" in" << endl;
  cout << "  the image.  The default is \"0 1 0\", the positive Y-axis." << endl;
};
#endif
static void usage()
{
#ifdef OUTPUT_TO_FILE
  cout << "vtk_write ";
#else
  cout << "vtk_view ";
#endif

  cout << "{parameter file name}" << endl;
  cout << endl;
  cout << "Parameters may include:" << endl;

  all_usage();

#ifdef OUTPUT_TO_FILE
  write_usage();
#else
  view_usage();
#endif
}

  
//----------------------------------------------------------------
// We link in one source file from zdrive.  It references some
// of the globals defined in another file.  We define the globals
// here to permit use of the source file we want.
//----------------------------------------------------------------

int Debug_Driver = 1;
int Number_Iterations = 1;
int Driver_Action = 1;  
int Chaco_In_Assign_Inv = 0;
struct Test_Flags Test;
struct Output_Flags Output;
