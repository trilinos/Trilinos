// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// File:      driver.cc                                                     //
// Project:   Local HSFC Ordering                                           //
// Author:    Michael Wolf                                                  //
// Date Started:  3/11/2010                                                 //
//                                                                          //
// Description:                                                             //
//              File tests local HSFC ordering for simple test problem      //
//              0 objects of process 0                                      //
//                                                                          //
// $Id: driver.cc 19 2009-12-18 17:52:36Z mmwolf $                          //
//////////////////////////////////////////////////////////////////////////////
#include "Zoltan_config.h"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include <iostream>
#include <math.h>
#include <vector>
#include <set>
#include <iterator>
#include <list>
#include <stdlib.h>


extern "C" {
#include <zoltan.h>
#include <zoltan_types.h>
}


#define SUCCESS 0
#define FAIL 1

int locNumObj = 16;   //number of objects locally

unsigned int answers[4][16] = {{12,  8,  9, 13, 14, 15, 11, 10,  6,  7,  3,  2,  1,  5,  4,  0},
                               {28, 24, 25, 29, 30, 31, 27, 26, 22, 23, 19, 18, 17, 21, 20, 16},
		               {44, 40, 41, 45, 46, 47, 43, 42, 38, 39, 35, 34, 33, 37, 36, 32},
		               {60, 56, 57, 61, 62, 63, 59, 58, 54, 55, 51, 50, 49, 53, 52, 48}};

int order(int worldsize, int myrank);
int checkResults(int myrank, ZOLTAN_ID_PTR permGIDs);

//////////////////////////////////////////////////////////////
// Zoltan query functions
//////////////////////////////////////////////////////////////
int zoltNumObjs(void *data,int *ierr);

void zoltGetObjs(void *data, int num_gid_entries, 
                      int num_lid_entries, ZOLTAN_ID_PTR global_ids, 
                      ZOLTAN_ID_PTR local_ids, int wgt_dim, 
		      float *obj_wgts, int *ierr);
int zoltNumGeom(void *data,int *ierr);
void zoltGeom(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, 
              ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr);
//////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int main(int argc, char *argv[])
{
  float ver;

#ifdef HAVE_MPI
  MPI_Init(&argc,&argv);
#endif

  Zoltan_Initialize(argc,argv,&ver); //initialize Zoltan

  int worldsize, myrank;

#ifdef HAVE_MPI
  MPI_Comm_size(MPI_COMM_WORLD,&worldsize);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  worldsize = 1;
  myrank = 0;
#endif

  if(worldsize!=1 && worldsize!=2 && worldsize!=3 && worldsize!=5)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif

    if(myrank==0)
    {
      std::cout << "WARNING for testPartitionMVInput() with number of processors = " << worldsize << std::endl
	        << "        Number of processes not supported .... (by default) returning success" << std::endl;
    }

    return SUCCESS;
  }

  int retFlag = 0;

  retFlag = order(worldsize,myrank);
 
#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  if (retFlag==0)
  {
    return FAIL;
  }

  return SUCCESS;
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
  int order(int worldsize, int myrank)
{
  struct Zoltan_Struct *zz;

  //  sprintf(paramPartitions,"%d",numprocs);

  zz = Zoltan_Create(MPI_COMM_WORLD);


  // General Zoltan Parameters
  Zoltan_Set_Param(zz,"ORDER_METHOD","LOCAL_HSFC");
  //Zoltan_Set_Param(zz,"ORDER_START_INDEX","1");


  // register query functions
  Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) zoltNumObjs, 0);
  Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) zoltGetObjs, 0);
  Zoltan_Set_Fn(zz, ZOLTAN_NUM_GEOM_FN_TYPE, (void (*)()) zoltNumGeom, 0);
  Zoltan_Set_Fn(zz, ZOLTAN_GEOM_FN_TYPE, (void (*)()) zoltGeom, 0);

  int numGidEntries=1;

  ZOLTAN_ID_PTR global_ids=0;
  ZOLTAN_ID_PTR permGIDs = 0;


  if(myrank != 0)
  {
    global_ids = new ZOLTAN_ID_TYPE[locNumObj];
    permGIDs = new ZOLTAN_ID_TYPE[locNumObj];

 
    for(int i=0;i<locNumObj; i++)
    {
      global_ids[i]=16*myrank + (i+1);
    }
    Zoltan_Order (zz, numGidEntries, locNumObj, global_ids, permGIDs);
  }
  else
  {
    global_ids = 0;
    permGIDs = 0;

     Zoltan_Order (zz, numGidEntries, 0, global_ids, permGIDs);

  }


  int successFlag = checkResults(myrank, permGIDs);

  int returnFlag;

#ifdef HAVE_MPI
  MPI_Allreduce(&successFlag,&returnFlag,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
#else
  returnFlag = successFlag;
#endif

  if(myrank==0)
  {
    if(returnFlag == 1)
    {
      std::cout   << "Success in testing local HSFC ordering on simple problem, numProcs = "
                  << worldsize << std::endl;
    }
    else
    {
      std::cout   << "Failure in testing local HSFC ordering on simple problem, numProcs = "
                  << worldsize << std::endl;
    }
  }

  delete [] global_ids;
  delete [] permGIDs;
  Zoltan_Destroy(&zz);

  return returnFlag;
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
int checkResults(int myrank, ZOLTAN_ID_PTR permGIDs)
{
  if(myrank==0)
  {
    std::cout << "HERE" << permGIDs << std::endl;
    return 1;
  }

   for(int i=0;i<locNumObj; i++)
   {
     if(permGIDs[i] != answers[myrank-1][i])
     {
       return 0;
     }
   }  
   return 1;
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int zoltNumObjs(void *data,int *ierr)
{
  int myrank;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=0;
#endif


  *ierr=0;

  if(myrank==0)
  {
    return 0;
  }

  return locNumObj;
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void zoltGetObjs(void *data, int num_gid_entries, 
                               int num_lid_entries, ZOLTAN_ID_PTR global_ids, 
                               ZOLTAN_ID_PTR local_ids, int wgt_dim, 
                               float *obj_wgts, int *ierr)
{
  int myrank;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank=0;
#endif

  if(myrank==0)
  {
    return;
  }

  for(int i=0;i<locNumObj;i++)
  {
    global_ids[i]=16*(myrank-1) + (i+1);
    local_ids[i]=i+1;
  }

  *ierr=0;

  return;
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
int zoltNumGeom(void *data,int *ierr)
{
  return 2;
}
//////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
void zoltGeom(void *data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_id, 
              ZOLTAN_ID_PTR local_id, double *geom_vec, int *ierr)
{

  int myrank;

#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
#else
  myrank = 0;
#endif

  if(myrank==0)
  {
    return;
  }

  geom_vec[0] = (double) ( 4*(((*global_id)-1)/32) + (((*global_id)-1)%4)+1);
  geom_vec[1] = (double)  ((((*global_id)-1)%32)/4 + 1);

  *ierr=0;

  return;
}
//////////////////////////////////////////////////////////////

