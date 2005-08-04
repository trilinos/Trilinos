/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2000,2001,2002, Sandia National Laboratories.               *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/

//--------------------------------------------------------------------------
// vtk_view is a parallel MPI program which may be used to visualize the    
// results of zdrive (Zoltan's test program) when it has calculated a       
// geometric partitioning.                                                   
//                                                                          
// The inputs to vtk_view are:                                              
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
// vtk_view can run with more or fewer processes than zdrive ran.  Choose   
// the size of your vtk_view application based only on the computational    
// demands of reading in and rendering the chaco or Nemesis files.
//                                                                          
// vtk_view requires that you obtain and compile VTK (the Visualization     
// Toolkit, available for free from Kitware, Inc. at www.vtk.org), version
// 5.0 or later.  In the  Config.{arch} file that directed your compilation 
// of Zoltan, add directives to indicate where the VTK, GL and X libraries 
// are.  See "Config.linux" for an example of these directives.
//
// Usage:  vtk_view parameterfileName
//
// The default parameter file name is "zdrive.inp" in the current working
// directory.
//--------------------------------------------------------------------------

// system headers

#include <dirent.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>

#include <map>
#include <iostream>

using namespace std;

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

#include "vtkCompositeRenderManager.h"
#include "vtkUnstructuredGrid.h"
#include "vtkDistributedDataFilter.h"
#include "vtkDataSetSurfaceFilter.h"
#include "vtkMultiProcessController.h"
#include "vtkMPIController.h"
#include "vtkMPICommunicator.h"
#include "vtkMPIGroup.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkActor.h"
#include "vtkScalarBarActor.h"
#include "vtkPolyDataMapper.h"
#include "vtkCamera.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkIntArray.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
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

// parameters from zdrive-style input file

#define LINELEN 80
#define MAXVAL 64
 
static char *vis_opt_names[UNDEFINED_LIST_MAX]={
  "zdrive count",     // number of zdrive processes
  "no rendering",     // don't bring up render window NOT IMPLEMENTED
  "output format",    // output this graphics format to file NOT IMPLEMENTED
  "output frames",    // output number of frames, default is 1 NOT IMPLEMENTED
  "output resolution" // resolution of output file NOT IMPLEMENTED
};

enum {option_zdrive_count,
      option_norender, 
      option_format, 
      option_frames, 
      option_resolution, 
      option_end};

static PARIO_INFO fname_opts;
static PROB_INFO prob_opts;
static UNDEFINED_INFO extra_options;
static char vis_opt_values[UNDEFINED_LIST_MAX][MAXVAL];

// Two ways to specify distributed files to parallel VTK ExodusII reader:

static char filePattern[ZMAXPATH]; // "printf"-like pattern: "%s/basename.0016.%04d"
static char filePrefix[ZMAXPATH];  //   root directory for pattern
static int fileRange[2];      //   file numbers for pattern
static char **fileNames;  // OR list of all file names 
static char **fileNamesBase;

// Number of processes that ran zdrive.  This is also the number of
// zdrive output files, and in the case of Exodus files, the number
// of input Exodus files.  If this is not provided in the input
// parameter file, vtk_view will try to figure it by looking in the 
// dirctory where the input files are.  If there are no zdrive output 
// files, vtk_view will just display the geometry.

static int numZdriveProcs = 0; 
static int numNemesisFiles = 0; 
static int lastPartition = 0; 

static int read_broadcast_input_options(int &argc, char **argv);
static int set_number_of_zdrive_processes(char *fname, char **disks);
static int set_nemesis_file_names_or_pattern(char *baseName,
  int numDisks, int diskListSize, int *diskList, int diskOffset, int useZeros,
  char *dirRoot, char *subDir);
static int check_partition_numbers(vtkUnstructuredGrid *ug);
static int assign_partition_numbers(vtkUnstructuredGrid *ug);
static void Run(vtkMultiProcessController *contr, void *arg);
static int checkAllrc(int rc, vtkMPICommunicator *comm);
static char *captionText(char *p, char *c);
static char *get_zdrive_output_file_name();

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

  vtkMultiProcessController *contr = vtkMultiProcessController::New();
  contr->Initialize(&argc, &argv);

  vtkMultiProcessController::SetGlobalController(contr);

  NumProcs = contr->GetNumberOfProcesses();
  Proc = contr->GetLocalProcessId();

  if (!contr->IsA("vtkMPIController"))
    {
    if (Proc == 0)
      {
      cout << "vtk_view requires MPI" << endl;
      }
    contr->Delete();
    return 1;
    }

  Comm = vtkMPICommunicator::SafeDownCast(contr->GetCommunicator());
  Controller = vtkMPIController::SafeDownCast(contr);

  // Read and broadcast the input parameters

  int rc = read_broadcast_input_options(argc, argv);

  if (rc){
    contr->Delete();
    return 1;
  }

  if (vis_opt_values[option_zdrive_count][0]){
    sscanf(vis_opt_values[option_zdrive_count], "%d", &numZdriveProcs);
    numNemesisFiles = numZdriveProcs;
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
      contr->Finalize(); 
      contr->Delete();
    
      return 1;
    }
  }
  else if (fname_opts.file_type == CHACO_FILE){

    // The chaco base name is in fname_opts.pexo_fname, but if we
    // are to read the zdrive output files, we need to know how many
    // zdrive processes there were.
  
    if (numZdriveProcs == 0){
      
      if (Proc == 0){
        rc = set_number_of_zdrive_processes(fname_opts.pexo_fname, NULL);
      }
  
      if (NumProcs > 1){
        int vals[2];
        vals[0] = rc;
        vals[1] = numZdriveProcs;
        Comm->Broadcast(vals, 2, 0);
 
        if (Proc > 0){
          rc = vals[0];
          numZdriveProcs = vals[1];
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
      cout << "in parameter file" << endl;
    }
    return 1;
  }

  // Run parallel method, VTK-style.

  contr->SetSingleMethod(Run, NULL);
  contr->SingleMethodExecute();

  contr->Finalize();
  contr->Delete();

  return 0;
}

static void Run(vtkMultiProcessController *contr, void *arg)
{
  // Read in the chaco or distributed Nemesis (exodusII) files.
  // Include at least one field array in the event we will not
  // have partition numbers to display.

  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();

  if (fname_opts.file_type == CHACO_FILE){
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
      cout << "Input dataset is too small, only ";
      cout << ug->GetNumberOfCells() << " cells." << endl;
    }
    ug->Delete();
    return;
  }

  // If we were processing zdrive output files, we'll visualize the
  // partition number assigned by zdrive.  If not, and there are
  // point or cell arrays other than global ID arrays, we'll visualize 
  // one of those.  If all else fails, we'll visualize the global element ID.

  char cellArray[128], pointArray[128];
  cellArray[0] = '\0';
  pointArray[0] = '\0';
  double range[2];

  if (numZdriveProcs > 0){
    // Find the partition number for each "cell" in my subgrid.  This
    // requires reading the zdrive output files.  We name the
    // resulting element or point array "Partition".

    if (fname_opts.file_type == CHACO_FILE){
      strcpy(pointArray, "Partition");
    }
    else{
      strcpy(cellArray, "Partition");
    }

    int rc = assign_partition_numbers(ug);
  
    if (NumProcs > 1){
      rc = checkAllrc(rc, Comm);
    }
  
    if (rc > 0){
      ug->Delete();
      return;
    }
  
    rc = check_partition_numbers(ug);
  
    if (rc){
      cout << Proc << " failed to obtain all partition numbers" << endl;
    }
  
    if (NumProcs > 1){
      rc = checkAllrc(rc, Comm);
    }
  
    if (rc > 0){
      ug->Delete();
      return;
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
      strcpy(pointArray, nm);

      ug->GetPointData()->GetArray(pointArray)->GetRange(range);
    }
    else{
      int ncarrays = ug->GetCellData()->GetNumberOfArrays();

      if (ncarrays > 1){
        const char *nm = ug->GetCellData()->GetArrayName(0);
        if (!strcmp(nm, "GlobalElementId")){
          nm = ug->GetCellData()->GetArrayName(1);
        }
        strcpy(cellArray, nm);
      }
      else{
        strcpy(cellArray, "GlobalElementId");
      }

      ug->GetCellData()->GetArray(cellArray)->GetRange(range);
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

  // Redistribute the cells for more load balanced rendering.  This
  // also creates ghost cells on each process if required for 
  // correct rendering.

  vtkDistributedDataFilter *dd = vtkDistributedDataFilter::New();

  dd->SetInput(ug);
  dd->SetController(contr);
  dd->SetGlobalElementIdArrayName("GlobalElementId");
  dd->SetGlobalNodeIdArrayName("GlobalNodeId");

  // Render the surface

  vtkDataSetSurfaceFilter *dss = vtkDataSetSurfaceFilter::New();
  dss->SetInput(dd->GetOutput());

  vtkPolyDataMapper *mapper = vtkPolyDataMapper::New();
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

  vtkScalarBarActor *sb = NULL;

  if (Proc == 0){
  
    sb = vtkScalarBarActor::New();
    sb->SetOrientationToVertical();
    if (pointArray[0]){
      sb->SetTitle(pointArray);
    }
    else{
      sb->SetTitle(cellArray);
    }
    sb->SetLookupTable(mapper->GetLookupTable());
//    sb->SetMaximumNumberOfColors(lastPartition + 1);
    sb->SetNumberOfLabels(4);
  }

  // Helpful text

  vtkTextActor *capActor = vtkTextActor::New();
  char *info = captionText(pointArray, cellArray);
  capActor->SetInput(info);
  delete [] info;

  capActor->SetAlignmentPoint(0);
  capActor->GetTextProperty()->SetVerticalJustificationToBottom();
  capActor->GetTextProperty()->SetJustificationToLeft();
  capActor->GetTextProperty()->SetFontSize(16);
  capActor->GetTextProperty()->BoldOn();
  capActor->GetTextProperty()->SetLineSpacing(1.2);

  vtkActor *actor = vtkActor::New();
  actor->SetMapper(mapper);

  // Set up the renderer

  vtkCompositeRenderManager *prm = vtkCompositeRenderManager::New();
  vtkRenderer *renderer = prm->MakeRenderer();

  renderer->AddActor(actor);

  if (Proc == 0){
    if (sb){
      renderer->AddActor(sb);
    }
    renderer->AddViewProp(capActor);
  }

  vtkRenderWindow *renWin = prm->MakeRenderWindow();
  renWin->AddRenderer(renderer);

  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);

  renderer->SetBackground(0,0,0);
  renWin->SetSize(400,400);
  renWin->SetPosition(0, 360*Proc);

  prm->SetRenderWindow(renWin);
  prm->SetController(contr);

  prm->InitializeOffScreen();   // Mesa GL only

  // We must update the whole pipeline here, otherwise node 0
  // goes into GetActiveCamera which updates the pipeline, putting
  // it into vtkDistributedDataFilter::Execute() which then hangs.
  // If it executes here, dd will be up-to-date won't have to 
  // execute in GetActiveCamera.

  mapper->SetPiece(Proc);
  mapper->SetNumberOfPieces(NumProcs);
  mapper->Update();

  prm->StartInteractor(); // now you can interact with window
  iren->Delete();

  mapper->Delete();
  actor->Delete();
  if (sb){
    sb->Delete();
  }
  capActor->Delete();
  renderer->Delete();
  renWin->Delete();

  dd->Delete();
  ug->Delete();
  dss->Delete();

  prm->Delete();
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
  //    "Zdrive Count" The number of zdrive processes.  vtk_view will try
  //       to figure it out if you don't supply it.  We don't assume it's
  //       equal to the number of vtk_view processes.
  // NOT IMPLEMENTED YET VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
  //    "No Rendering"  By default vtk_view creates a window on your
  //       display to render the partitioning created by zdrive.  This
  //       tells vtk_view to not do this
  //    "Output Format"  Normally vtk_view does not write an image
  //       file depicting the partitioning.  Set this option to
  //       "jpeg", blah, blah for an image file
  //    "Output Frames"  By default vtk_view outputs one frame to the
  //       file if you specified an "Output Format".  For an animation
  //       with several frames, state the number of frames you would 
  //       like here.
  //    "Output Resolution"  By default, if you have specified an
  //       "Output Format", vtk_view will write 300x300 images.  If you
  //       would like higher resolution images, state so here.

  const char  *cmd_file;

  switch(argc)
  {
  case 1:
    cmd_file = "zdrive.inp";
    break;

  case 2:
    cmd_file = argv[1];
    break;

  default:
    cerr << "MAIN: ERROR in command line " ;

    if (Proc == 0)
    {
      cerr << "usage:" << endl;
      cerr << "\t" << argv[0] << " [command file]";
    }
    cerr << endl;

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
                while (*c && isblank(*c)) c++;
                if (*c){
                  strncpy(vis_opt_values[j], c, MAXVAL-1);
                  input_ok = 1;
                }
              }
            }
          }
        } 
      }
    }

    if (input_ok && 
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

  return 0;
}

static int set_nemesis_file_names_or_pattern(char *baseName,
  int numDisks, int diskListSize, int *diskList, 
  int diskOffset, int useZeros,
  char *dirRoot, char *subDir)
{
  char **disks=NULL;
  int nzeros, len, rc, num;
  int retVal = 0;

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

  if (numZdriveProcs == 0){
    
    if (Proc == 0){
      // also sets numNemesisFiles, which is usually the same
      rc = set_number_of_zdrive_processes(baseName, disks);
    }

    if (NumProcs > 1){
      int counts[3];
      counts[0] = rc;
      counts[1] = numZdriveProcs;
      counts[2] = numNemesisFiles;
      Comm->Broadcast(counts, 3, 0);

      if (Proc > 0){
        rc = counts[0];
        numZdriveProcs = counts[1];
        numNemesisFiles = counts[2];
      } 
    }

    if ((numNemesisFiles == 0) && (Proc == 0)){
      cerr << "Error: Unable to locate input Exodus/Nemesis files" << endl;
    }
    if ((numNemesisFiles == 0) || rc){
      retVal = 1;
      goto done2;
    }
  }

  num = 1+(int)(std::log10((float)numNemesisFiles));

  if (numDisks > 1){

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

    if (numDisks == 1){
      strcpy(filePrefix, disks[0]);
      sprintf(filePattern, "%%s/%s.%d.%%0%dd", baseName, numNemesisFiles, num);
    }
    else{
      strcpy(filePrefix, baseName);
      sprintf(filePattern, "%%s.%d.%%0%dd", numNemesisFiles, num);
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
    numZdriveProcs = count;
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
      }
    }

    numZdriveProcs = 0;
    numNemesisFiles = count;
  }

  closedir(dir);

  if ((numZdriveProcs == 0) && verbose){
    cerr << "Unable to locate zdrive output file(s) ";
    cerr << dn << "/" << fn << ".out.{numfiles}.{fileno}" << endl;
    cerr << "We'll just show you the geometry of the input files." << endl;
  }

  delete [] dn;
  delete [] fn;

  return 0;
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
        ids[i] = it->second;
      }
    }
  }
}
static char *get_zdrive_output_file_name()
{
  char *nm = new char [ZMAXPATH];

  sprintf(nm, "%s.out.%d.0",
    fileNames ? fileNamesBase[0] : fname_opts.pexo_fname, numZdriveProcs);

  return nm;
}
static int read_zdrive_output(int from, int to)
{
char nm[ZMAXPATH], line[1024], *c;
FILE *fp;
int g, p, nvals;

  listSize = 0;

  int num = 1+(int)(std::log10((float)numZdriveProcs));

  if (!fileNames){
    c = fname_opts.pexo_fname;
  }

  for (int i=from; i<=to; i++){

    if (fileNames){
      c = fileNamesBase[i];
    }

    sprintf(nm, "%s.out.%d.%0*d",c,numZdriveProcs,num,i);

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

  if (fname_opts.file_type == CHACO_FILE){
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

  if (fname_opts.file_type == CHACO_FILE){
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
    rc = read_zdrive_output(0, numZdriveProcs-1);
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

  nfiles = numZdriveProcs / nprocs;
  extra = numZdriveProcs - (nfiles * nprocs);
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
static char *captionText(char *pnm, char *cnm)
{
  char *buf1 = new char [LINELEN  + 1];
  char *buf = new char [LINELEN * 3 + 1];
  char *c = buf;
  int lenleft = LINELEN * 3;
  int linelen = 0;

  strcpy(buf, "(key \"t\" toggles interaction mode, \"r\" resets, \"q\" quits)\n");

  int used = strlen(buf);
  lenleft -= used;
  c += used;

  if (numZdriveProcs == 0){
    used = snprintf(buf1, LINELEN, "%s: %s\n",
      (fname_opts.file_type == CHACO_FILE ? "Chaco" : "Exodus/Nemesis"),
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
    if ((fname_opts.file_type == CHACO_FILE) && (dtype != INITIAL_FILE)) {
      snprintf(buf1, LINELEN, "Chaco: %s, %s initial distribution, %s",
        fname_opts.pexo_fname, 
        (dtype == INITIAL_LINEAR ? 
           "linear" : 
           (dtype == INITIAL_CYCLIC ? "cyclic" : "owner")),
        prob_opts.method);
    }
    else{
      snprintf(buf1, LINELEN, "%s: %s, %s",
        (fname_opts.file_type == CHACO_FILE ? "Chaco" : "Exodus/Nemesis"),
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
