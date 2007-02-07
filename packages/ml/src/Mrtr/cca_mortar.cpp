#ifdef TRILINOS_PACKAGE
#ifdef MORTAR


#include <ctime>
#include <cstdlib>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Teuchos_ParameterList.hpp"

// include the mortar stuff
#include "mrtr_functions.H"
#include "mrtr_segment.H"
#include "mrtr_segment_linear1D.H"
#include "mrtr_segment_bilinearquad.H"
#include "mrtr_node.H"
#include "mrtr_interface.H"
#include "mrtr_manager.H"

// user defined headers
#include "cca_mortar.H"


// the Epetra_Comm
#ifdef PARALLEL
static class Epetra_MpiComm* comm;
#else
static class Epetra_SerialComm* comm;
#endif

// the mortar manager class
static class MOERTEL::Manager* mrtr_manager;
class Epetra_CrsMatrix* saddleproblem;

static const int outlevel = 10;

/*----------------------------------------------------------------------*
 |  routine to create mortar interfaces                      m.gee 06/05|
 *----------------------------------------------------------------------*/
#ifdef PARALLEL
int create_mortar(FIELD *actfield, PARTITION *actpart,
                  int disnum, DESIGN *design, int init, MPI_Comm mpicomm)
#else
int create_mortar(FIELD *actfield, PARTITION *actpart,
                  int disnum, DESIGN *design, int init)
#endif
{
  // Create a communicator for Epetra objects
#ifdef PARALLEL
  comm = new Epetra_MpiComm(mpicomm);
  int MyPID = comm->MyPID();
#else
  comm = new Epetra_SerialComm();
  int MyPID = comm->MyPID();
#endif

  // output level (0-10)

  //-------------------------------------------------------------------
  // get discretization

  DISCRET* actdis = &(actfield->dis[disnum]);

  //-------------------------------------------------------------------
  // create the mortar manager

  mrtr_manager = new MOERTEL::Manager(*comm,outlevel);

  //=================================== handle the line interfaces (2D)
  
  //-------------------------------------------------------------------
  // count number of 2D mortar interfaces ninter

  int* ids    = NULL;
  int  ninter = cca_mrtr_2D_numinterfaces(&ids,design);
  if (ninter) 
    mrtr_manager->SetDimension(MOERTEL::Manager::manager_2D);
  //-------------------------------------------------------------------
  // loop over number of mortar interfaces 

  MOERTEL::Interface* interface = NULL;
  for (int k=0; k<ninter; ++k)
  {
    //-----------------------------------------------------------------
    // create an interface class

    int id    = ids[k];
    interface = new MOERTEL::Interface(id,true,mrtr_manager->Comm(),outlevel);

    //-----------------------------------------------------------------
    // find dlines for given interface id

    DLINE* dline1 = NULL;
    DLINE* dline2 = NULL;
    bool ok = cca_mrtr_2D_finddlines(id,&dline1,&dline2,design);
    
    
    //-----------------------------------------------------------------
    // find all glines that are on dline1/dline2

    int ngline1=0;
    int ngline2=0;
    GLINE** gline1 = NULL; // needs to be deleted
    GLINE** gline2 = NULL; // needs to be deleted
    ngline1 = cca_mrtr_2D_find_glines_on_dline(&gline1,dline1,actdis);    
    ngline2 = cca_mrtr_2D_find_glines_on_dline(&gline2,dline2,actdis);    
    
    // loop all lines and create segment classes from them

    // side 0 of interface
    for (int i=0; i<ngline1; ++i)
    {
      // only my segments
      if (gline1[i]->proc != MyPID) continue;
      if (gline1[i]->ngsurf != 1)
        cout << "***WRN*** gline " << i << " has more then one gsurf!\n";
      int* nodeIds = NULL; // new int[gline1[i]->ngnode];
      
      MOERTEL::Segment::SegmentType typ = MOERTEL::Segment::seg_none;
      cca_mrtr_2D_prepare_gline_data(gline1[i],&nodeIds,&typ);  

      MOERTEL::Segment* seg;
      if (typ==MOERTEL::Segment::seg_Linear1D)
        seg = new MOERTEL::Segment_Linear1D(gline1[i]->Id,gline1[i]->ngnode,nodeIds,outlevel);
      else dserror("Unknown type of 1D segment");
      
      delete [] nodeIds;
      
      interface->AddSegment(*seg,0);
      delete seg; seg = NULL;
    }

    // side 1 of interface
    for (int i=0; i<ngline2; ++i)
    {
      // only my segments
      if (gline2[i]->proc != MyPID) continue;
      if (gline2[i]->ngsurf != 1)
        cout << "***WRN*** gline " << i << " has more then one gsurf!\n";

      int* nodeIds = NULL;
      MOERTEL::Segment::SegmentType typ = MOERTEL::Segment::seg_none;
      cca_mrtr_2D_prepare_gline_data(gline2[i],&nodeIds,&typ); 

      MOERTEL::Segment* seg;
      if (typ==MOERTEL::Segment::seg_Linear1D)
        seg = new MOERTEL::Segment_Linear1D(gline2[i]->Id,gline2[i]->ngnode,nodeIds,outlevel);
      else dserror("Unknown type of 1D segment");
      
      delete [] nodeIds;

      interface->AddSegment(*seg,1);
      delete seg; seg = NULL;
    }
    
    //-----------------------------------------------------------------
    // find all gnodes that are on dline1/dline2

    int ngnode1 = 0;
    int ngnode2 = 0;
    GNODE** gnode1 = NULL; // needs to be deleted
    GNODE** gnode2 = NULL; // needs to be deleted
    bool*   boundary1 = NULL;
    bool*   boundary2 = NULL;
    ngnode1 = cca_mrtr_2D_find_gnodes_on_dline(&gnode1,&boundary1,dline1,actdis);
    ngnode2 = cca_mrtr_2D_find_gnodes_on_dline(&gnode2,&boundary2,dline2,actdis);
    
    // side 0 of interface
    for (int i=0; i<ngnode1; ++i)
    {
      // only my own nodes
      if (gnode1[i]->node->proc != MyPID) continue;
      MOERTEL::Node* node = 
        new MOERTEL::Node(gnode1[i]->Id,gnode1[i]->node->x,
                       2/*gnode1[i]->node->numdf*/,gnode1[i]->node->dof,boundary1[i],outlevel);
      interface->AddNode(*node,0);
      delete node; node = NULL;
    }
    
    // side 1 of interface
    for (int i=0; i<ngnode2; ++i)
    {
      // only my own nodes
      if (gnode2[i]->node->proc != MyPID) continue;
      MOERTEL::Node* node = 
        new MOERTEL::Node(gnode2[i]->Id,gnode2[i]->node->x,
                       2/*gnode2[i]->node->numdf*/,gnode2[i]->node->dof,boundary2[i],outlevel);
      interface->AddNode(*node,1);
      delete node; node = NULL;
    }
    
    //-----------------------------------------------------------------
    // manually choose mortar (master side)
    // mortar side is either 0 or 1 or -2 for automatic
    interface->SetMortarSide(-2);

    //-----------------------------------------------------------------
    // set linear shape functions on both sides, 
    // set additional dual shape functions on non-mortar side 
    // Currently ONLY linear functions!!!!
    if (interface->MortarSide() != -2)
    {
      MOERTEL::Function_Linear1D* func = new MOERTEL::Function_Linear1D(outlevel);
      interface->SetFunctionAllSegmentsSide(0,0,func);
      interface->SetFunctionAllSegmentsSide(1,0,func);
      delete func; func = NULL;
    }
    // mortar side is not yet chosen, we cannot set functions.
    // So we just set what kind of functions we want to have set
    // The Mortar::Manager is setting the functions later accordingly
    else
    {
      interface->SetFunctionTypes(MOERTEL::Function::func_Linear1D,      // the isoparametric function
                                  MOERTEL::Function::func_DualLinear1D); // the LM space
                                  //MOERTEL::Function::func_Linear1D); // the LM space
    }

    //-----------------------------------------------------------------
    // get the mortar side and set dual linear shape function to non mortar side

    int side = interface->MortarSide();
    if (side==1 || side==0)
    {
      side     = interface->OtherSide(side);
      //MOERTEL::Function_Linear1D* func = new MOERTEL::Function_Linear1D(outlevel);
      MOERTEL::Function_DualLinear1D* func = new MOERTEL::Function_DualLinear1D(outlevel);
      interface->SetFunctionAllSegmentsSide(side,1,func);
      delete func; func = NULL;
    }
    
    //-----------------------------------------------------------------
    // set type of projection to be used on this interface

    interface->SetProjectionType(MOERTEL::Interface::proj_continousnormalfield);    
    //interface->SetProjectionType(MOERTEL::Interface::proj_orthogonal);    
    
    //-----------------------------------------------------------------
    // Complete interface 

    ok = interface->Complete();
    if (!ok)
    {
      cout << "***ERR*** interface->Complete() returned false\n";
      exit(EXIT_FAILURE);
    }
    
    //-----------------------------------------------------------------
    // Add interface to Mortar Manager 
    mrtr_manager->AddInterface(*interface); // deep copy

    //-----------------------------------------------------------------
    // tidy up 

    if (gline1) delete [] gline1;
    if (gline2) delete [] gline2;
    if (gnode1) delete [] gnode1;
    if (gnode2) delete [] gnode2;
    if (boundary1) delete [] boundary1;
    if (boundary2) delete [] boundary2;
    if (interface) delete interface; interface = NULL;
    
  //-------------------------------------------------------------------
  } // for (int k=0; k<ninter; ++k) loop over mortar interfaces
  delete [] ids;

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  //-------------------------------------------------------------------

  //-------------------------------------------------------------------
  // count number of 3D mortar interfaces 
  ninter = cca_mrtr_3D_numinterfaces(&ids,design);
  if (ninter) 
    mrtr_manager->SetDimension(MOERTEL::Manager::manager_3D);
  //-------------------------------------------------------------------
  // loop over number of mortar interfaces 
  for (int k=0; k<ninter; ++k)
  {
    
    //-----------------------------------------------------------------
    // create an interface class
    
    int id = ids[k];
    interface = new MOERTEL::Interface(id,false,mrtr_manager->Comm(),outlevel);
  
    //-----------------------------------------------------------------
    // find dsurfs for given interface id

    DSURF* dsurf1 = NULL;
    DSURF* dsurf2 = NULL;
    bool ok = cca_mrtr_3D_finddsurfs(id,&dsurf1,&dsurf2,design);
  
    //-----------------------------------------------------------------
    // find all gsurfs on dsurf1/dsurf2

    int ngsurf1=0;
    int ngsurf2=0;
    GSURF** gsurf1 = NULL;  
    GSURF** gsurf2 = NULL;  
    ngsurf1 = cca_mrtr_3D_find_gsurfs_on_dsurf(&gsurf1,dsurf1,actdis);  
    ngsurf2 = cca_mrtr_3D_find_gsurfs_on_dsurf(&gsurf2,dsurf2,actdis);  
  
    // loop all gsurfs and create segment classes from them
    
    // side 0 of interface
    for (int i=0; i<ngsurf1; ++i)
    {
      // only my own segments
      if (gsurf1[i]->ngvol>1)
        dserror("ERR: gsurf attached to more then one gvol");
      if (gsurf1[i]->gvol[0]->element->proc != MyPID) continue;
      
      int* nodeIds = NULL;
      MOERTEL::Segment::SegmentType typ = MOERTEL::Segment::seg_none;
      if (!cca_mrtr_3D_prepare_gsurf_data(gsurf1[i],&nodeIds,&typ))
        dserror("ERR: cannot find nodal topology");
      
      MOERTEL::Segment* seg;
      if (typ==MOERTEL::Segment::seg_BiLinearQuad)
        seg = new MOERTEL::Segment_BiLinearQuad(gsurf1[i]->Id,gsurf1[i]->ngnode,nodeIds,outlevel);  
      else
        dserror("Unknown type of 2D element"); 
      delete [] nodeIds;
      interface->AddSegment(*seg,0); 
      delete seg; seg = NULL;      
    }
    
    // side 1 of interface
    for (int i=0; i<ngsurf2; ++i)
    {
      // only my own segments
      if (gsurf2[i]->ngvol>1)
        dserror("ERR: gsurf attached to more then one gvol");
      if (gsurf2[i]->gvol[0]->element->proc != MyPID) continue;
      
      int* nodeIds = NULL;
      MOERTEL::Segment::SegmentType typ = MOERTEL::Segment::seg_none;
      if (!cca_mrtr_3D_prepare_gsurf_data(gsurf2[i],&nodeIds,&typ))
        dserror("ERR: cannot find nodal topology");
      
      MOERTEL::Segment* seg;
      if (typ==MOERTEL::Segment::seg_BiLinearQuad)
        seg = new MOERTEL::Segment_BiLinearQuad(gsurf2[i]->Id,gsurf2[i]->ngnode,nodeIds,outlevel);
      else
        dserror("Unknown type of 2D element"); 
      delete [] nodeIds;
      interface->AddSegment(*seg,1); 
      delete seg; seg = NULL;      
    }
  
    //-----------------------------------------------------------------
    // find all gnodes on dsurf1/dsurf2
    
    int ngnode1 = 0;
    int ngnode2 = 0;
    GNODE** gnode1 = NULL;
    GNODE** gnode2 = NULL;
    bool*   boundary1 = NULL;
    bool*   boundary2 = NULL;
    
    ngnode1 = cca_mrtr_3D_find_gnodes_on_dsurf(&gnode1,&boundary1,dsurf1,actdis);
    ngnode2 = cca_mrtr_3D_find_gnodes_on_dsurf(&gnode2,&boundary2,dsurf2,actdis);

    // side 0 of interface
    for (int i=0; i<ngnode1; ++i)
    {
      // only my own nodes
      if (gnode1[i]->node->proc != MyPID) continue;
      MOERTEL::Node* node = 
        new MOERTEL::Node(gnode1[i]->node->Id,gnode1[i]->node->x,
                       gnode1[i]->node->numdf,gnode1[i]->node->dof,
                       boundary1[i],outlevel);
      interface->AddNode(*node,0);
      delete node; node = NULL;
    }
  
    // side 1 of interface
    for (int i=0; i<ngnode2; ++i)
    {
      // only my own nodes
      if (gnode2[i]->node->proc != MyPID) continue;
      MOERTEL::Node* node = 
        new MOERTEL::Node(gnode2[i]->node->Id,gnode2[i]->node->x,
                       gnode2[i]->node->numdf,gnode2[i]->node->dof,
                       boundary2[i],outlevel);
                       
      interface->AddNode(*node,1);
      delete node; node = NULL;
    }
  
    //-----------------------------------------------------------------
    // manually choose mortar (master side)
    // mortar side is either 0 or 1 or -2 for automatic
    interface->SetMortarSide(0);

    //-----------------------------------------------------------------
    // set linear shape functions on both sides, 
    // set additional dual shape functions on non-mortar side 
    // Currently ONLY linear functions!!!!
    if (interface->MortarSide() != -2)
    {
      MOERTEL::Function_LinearTri* func = new MOERTEL::Function_LinearTri(outlevel);
      interface->SetFunctionAllSegmentsSide(0,0,func);
      interface->SetFunctionAllSegmentsSide(1,0,func);
      delete func; func = NULL;
    }
    // mortar side is not yet chosen, we cannot set functions.
    // So we just set what kind of functions we want to have set
    // The Mortar::Manager is setting the functions later accordingly
    else
    {
      interface->SetFunctionTypes(MOERTEL::Function::func_LinearTri,      // the isoparametric function
                                  MOERTEL::Function::func_DualLinearTri); // the LM space
                                  //MOERTEL::Function::func_LinearTri); // the LM space
                                  //MOERTEL::Function::func_ConstantTri); // the LM space
    }
      
    //-----------------------------------------------------------------
    // get the mortar side and set dual linear shape function to non mortar side

    int side = interface->MortarSide();
    if (side==1 || side==0)
    {
      side     = interface->OtherSide(side);
      MOERTEL::Function_DualLinearTri* func = new MOERTEL::Function_DualLinearTri(outlevel);
      //MOERTEL::Function_LinearTri* func = new MOERTEL::Function_LinearTri(outlevel);
      interface->SetFunctionAllSegmentsSide(side,1,func);
      delete func; func = NULL;
    }
    
    //-----------------------------------------------------------------
    // set type of projection to be used on this interface

    interface->SetProjectionType(MOERTEL::Interface::proj_continousnormalfield);    
    
    //-----------------------------------------------------------------
    // Complete interface 

    ok = interface->Complete();
    if (!ok)
    {
      cout << "***ERR*** interface->Complete() returned false\n";
      exit(EXIT_FAILURE);
    }
    
    //-----------------------------------------------------------------
    // Add interface to Mortar Manager 
    mrtr_manager->AddInterface(*interface); // deep copy

    //-----------------------------------------------------------------
    // tidy up 

    if (gsurf1) delete [] gsurf1;
    if (gsurf2) delete [] gsurf2;
    if (gnode1) delete [] gnode1;
    if (gnode2) delete [] gnode2;
    if (boundary1) delete [] boundary1;
    if (boundary2) delete [] boundary2;
    if (interface) delete interface; interface = NULL;
    
  //-------------------------------------------------------------------
  }// for (int k=0; k<ninter; ++k) loop over mortar interfaces
  if (ids) delete [] ids;
  
#if 0
  //-----------------------------------------------------------------
  // print all interfaces
  fflush(stdout);
  comm->Barrier();
  cout << *mrtr_manager;
#endif  

  return (1);
}




/*----------------------------------------------------------------------*
 |  routine to compute mortar matrices                       m.gee 06/05|
 *----------------------------------------------------------------------*/
int compute_mortar(SPARSE_TYP* arraytyp, SPARSE_ARRAY* array, DESIGN *design)
{
  //-------------------------------------------------------------------
  // count number of 2D mortar interfaces ninter
  //-------------------------------------------------------------------
  int* ids    = NULL;
  int  ninter = cca_mrtr_2D_numinterfaces(&ids,design);
  if (ninter==0)
    if (ids) delete [] ids;

  ninter += cca_mrtr_3D_numinterfaces(&ids,design);
  if (ninter==0)
  {
    if (ids) delete [] ids;
    return 1;
  }
  
  //-------------------------------------------------------------------
  // check type of matrix, currently only msr/spooles supported
  int *update = NULL;
  int numeq = 0;
  int numeq_total = 0;
  if (*arraytyp == msr)
  {
    update = array->msr->update.a.iv;
    numeq  = array->msr->numeq;
    numeq_total  = array->msr->numeq_total;
  }
  else if (*arraytyp == spoolmatrix)
  {
    update = array->spo->update.a.iv;
    numeq  = array->spo->numeq;
    numeq_total  = array->spo->numeq_total;
  }
  else
  {
    cout << "Mortar works with Aztec msr or Spooles matrix only!\n";
    exit(EXIT_FAILURE);
  }
  if (!update)
  {
    cout << "update=NULL\n";
    exit(EXIT_FAILURE);
  }
  if (!mrtr_manager)
  {
    cout << "mrtr_manager=NULL\n";
    exit(EXIT_FAILURE);
  }
  
  //-------------------------------------------------------------------
  // create a map for the current input matrix  
  Epetra_Map* matrixmap = new Epetra_Map(numeq_total,numeq,update,0,*comm);

  //-------------------------------------------------------------------
  // store this map in the mortar manager
  mrtr_manager->SetProblemMap(matrixmap);
  
  //-------------------------------------------------------------------
  // integrate the whole mortar interfaces
  mrtr_manager->Mortar_Integrate();  
      
  //-------------------------------------------------------------------
  // create an Epetra_CrsMatrix from the inputmatrix
  Epetra_CrsMatrix* inputmatrix = new Epetra_CrsMatrix(Copy,*matrixmap,60);
  
  //-------------------------------------------------------------------
  // fill the inputmatrix from ccarat depending on msr or spooles matrix
  if (*arraytyp == msr)
  {
    update = array->msr->update.a.iv;
    numeq  = array->msr->numeq;
    int* bindx = array->msr->bindx.a.iv;
    double* val = array->msr->val.a.dv;
    for (int i=0; i<numeq; ++i)
    {
      int grow = update[i];
      // diagonal
      inputmatrix->InsertGlobalValues(grow,1,&(val[i]),&grow);
      // off-diagonal
      inputmatrix->InsertGlobalValues(grow,bindx[i+1]-bindx[i],&(val[bindx[i]]),&(bindx[bindx[i]]));
    }
  }
  else if (*arraytyp == spoolmatrix)
  {
    int* irn = array->spo->irn_loc.a.iv;
    int* jcn = array->spo->jcn_loc.a.iv;
    double* val = array->spo->A_loc.a.dv;
    int nnz = array->spo->A_loc.fdim;
    for (int i=0; i<nnz; ++i)
      inputmatrix->InsertGlobalValues(irn[i],1,&(val[i]),&jcn[i]);
  }
  inputmatrix->FillComplete();
  inputmatrix->OptimizeStorage();

  //-------------------------------------------------------------------
  // store the inputmatrix in the Mortar manager
  mrtr_manager->SetInputMatrix(inputmatrix,false);

#if 0
  //-------------------------------------------------------------------
  // produce the saddle point problem and use spooles
  // (Need to adjust defines in solver_spooles.c as well)
  if (saddleproblem) 
    delete saddleproblem; 
  saddleproblem = mrtr_manager->MakeSaddleProblem();

  //-------------------------------------------------------------------
  // copy the saddle point problem matrix back to ccarat in msr or spooles
  if (*arraytyp == msr)
  {
    cout << "***ERR** msr matrix not yet impl\n\n"; fflush(stdout);
  }
  else if (*arraytyp == spoolmatrix)
  {
    array->spo->numeq       = saddleproblem->NumMyRows();
    array->spo->numeq_total = saddleproblem->NumGlobalRows();

    amdel(&(array->spo->update));
    update = (int*)amdef("update",&(array->spo->update),array->spo->numeq,1,"IV");
    for (int i=0; i<array->spo->numeq; ++i) 
    {
      update[i] = saddleproblem->GRID(i);
      if (update[i] < 0) {
        cout << "update[i]<0!\n";
        exit(EXIT_FAILURE); 
      }
    }
    array->spo->nnz = saddleproblem->NumMyNonzeros();
    amdel(&(array->spo->rowptr));
    amdel(&(array->spo->irn_loc));
    amdel(&(array->spo->jcn_loc));
    amdel(&(array->spo->A_loc));
    int* rowptr   = (int*)   amdef("rowptr" ,&(array->spo->rowptr) ,array->spo->numeq+1,1,"IV");
    int* irn_loc  = (int*)   amdef("irn_loc",&(array->spo->irn_loc),array->spo->nnz,1,"IV");
    int* jcn_loc  = (int*)   amdef("jcn_loc",&(array->spo->jcn_loc),array->spo->nnz,1,"IV");
    double* A_loc = (double*)amdef("A_loc"  ,&(array->spo->A_loc)  ,array->spo->nnz,1,"DV");
  
    //saddleproblem->ColMap();
    int rowcount = 0;
    int nnzcount = 0;
    //FILE *out = fopen("matrix.txt","w");
    //fprintf(out,"nrow %d nnz %d\n",saddleproblem->NumMyRows(),saddleproblem->NumMyNonzeros());
    for (int i=0; i<array->spo->numeq; ++i)
    {
      
      int Row = saddleproblem->GRID(i);
      if (Row<0) {
        cout << "Row<0!\n";
        exit(EXIT_FAILURE); 
      }
      int NumEntries; double* Values; int* Indices;
      int err = saddleproblem->ExtractMyRowView(i,NumEntries,Values,Indices);
      if (err) {
        cout << "ExtractMyRowView returned " << err << endl;
        exit(EXIT_FAILURE); 
      }
      for (int j=0; j<NumEntries; ++j)
      {
        //fprintf(out,"%10d %10d %20.10e\n",Row+1,saddleproblem->GCID(Indices[j])+1,Values[j]);
        irn_loc[nnzcount] = Row;
        jcn_loc[nnzcount] = saddleproblem->GCID(Indices[j]);
        if (jcn_loc[nnzcount]<0) {
          cout << "saddleproblem->GCID returned " << jcn_loc[nnzcount] << endl;
          exit(EXIT_FAILURE); 
        }
        A_loc[nnzcount]   = Values[j];
        ++nnzcount;
      }
      rowptr[i] = rowcount;
      rowcount += NumEntries;
    } // for (int i=0; i<array->spo->numeq; ++i)
    rowptr[array->spo->numeq] = rowcount;
  }
#endif  
  
  //-------------------------------------------------------------------
  // tidy up
  delete matrixmap;
  return(1);
}

/*----------------------------------------------------------------------*
 |  routine to solve mortar matrices                         m.gee 12/05|
 *----------------------------------------------------------------------*/
int solve_mortar(struct _DIST_VECTOR *sol, struct _DIST_VECTOR *rhs)
{
  //-------------------------------------------------------------------
  // create Epetra_Vectors from rhs and sol
  Epetra_Map* rowmap = mrtr_manager->ProblemMap();  
  Epetra_Vector esol(*rowmap,true);
  esol.Random();
  Epetra_Vector erhs(*rowmap,false);
  for (int i=0; i<erhs.MyLength(); ++i)
    erhs[i] = rhs->vec.a.dv[i];
  
  //-------------------------------------------------------------------
  // create a Teuchos Parameter List holding solver arguments
  Teuchos::ParameterList params;
  //params.set("System","SaddleSystem");
  params.set("System","SPDSystem");

  // choose solver package
  //params.set("Solver","Amesos");
  params.set("Solver","ML/Aztec");

  // argument sublist for amesos
  Teuchos::ParameterList& amesosparams = params.sublist("Amesos");
  amesosparams.set("Solver","Amesos_Klu");
  amesosparams.set("PrintTiming",true);
  amesosparams.set("PrintStatus",true);
  amesosparams.set("UseTranspose",true);
  
  // argument sublist for aztec
  Teuchos::ParameterList& aztecparams = params.sublist("Aztec");
  aztecparams.set("AZ_solver","AZ_cg");
  aztecparams.set("AZ_precond","AZ_user_precond");
  aztecparams.set("AZ_max_iter",1200);
  aztecparams.set("AZ_output",100);
  aztecparams.set("AZ_tol",1.0e-7);
  aztecparams.set("AZ_scaling","AZ_none");

  // argument sublist for ml
  Teuchos::ParameterList& mlparams = params.sublist("ML");
  ML_Epetra::SetDefaults("SA",mlparams);
  mlparams.set("ML output",10);
  mlparams.set("print unused",1/*-2*/);
  mlparams.set("increasing or decreasing","increasing");
  mlparams.set("PDE equations",3);
  mlparams.set("max levels",3);
  mlparams.set("aggregation: type","Uncoupled");
  mlparams.set("aggregation: damping factor",1.33);
  mlparams.set("coarse: max size",15);
  mlparams.set("coarse: type","Amesos-KLU"); // Amesos-KLU MLS
  mlparams.set("smoother: type","ML MLS"); /* MLS ML symmetric Gauss-Seidel ML Gauss-Seidel*/
  mlparams.set("smoother: MLS polynomial order",3);
  mlparams.set("smoother: damping factor",0.67);
  mlparams.set("smoother: sweeps",1);
  mlparams.set("smoother: pre or post","both");
  mlparams.set("null space: type","pre-computed");
  mlparams.set("null space: add default vectors",false);
  int dimnullspace = 6;
  int nummyrows = mrtr_manager->ProblemMap()->NumMyElements();
  int numglobalrows = mrtr_manager->ProblemMap()->NumGlobalElements();
  int* update       = mrtr_manager->ProblemMap()->MyGlobalElements();
  int dimnsp        = dimnullspace*nummyrows;
  double* nsp   = new double[dimnsp];
  int numdf = 0;
  solver_ml_nullspace(nsp,dimnsp,update,1,NULL,nummyrows,numglobalrows,&numdf);
  mlparams.set("null space: dimension",dimnullspace);
  mlparams.set("null space: vectors",nsp);

  //-------------------------------------------------------------------
  // solve
  bool ok = mrtr_manager->Solve(params,esol,erhs); 
  if (!ok) cout << "***ERR*** MOERTEL::Manager::Solve() returned false\n";

  mrtr_manager->ResetSolver();

  //-------------------------------------------------------------------
  // copy result back
  for (int i=0; i<esol.MyLength(); ++i)
    sol->vec.a.dv[i] = esol[i];

  delete mrtr_manager; mrtr_manager = 0;
  delete [] nsp; nsp = NULL;
  
  return (int)ok;
}
#endif // MORTAR
#endif // TRILINOS_PACKAGE
