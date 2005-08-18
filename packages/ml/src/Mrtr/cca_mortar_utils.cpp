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

// user defined headers
#include "cca_mortar.H"

/*----------------------------------------------------------------------*
 |  counter number of 2D mortar interfaces                   m.gee 06/05|
 *----------------------------------------------------------------------*/
int cca_mrtr_2D_numinterfaces(int **iids, DESIGN* design)
{
  int size   = 50;
  int ninter = 0;
  *iids = new int[size];
  int* ids = *iids;
  for (int i=0; i<design->ndline; ++i)
  {
    if (design->dline[i].hasmortar==0) continue;
    int id = design->dline[i].mrtr_Id;
    bool foundit = false;
    for (int j=0; j<ninter; ++j)
      if (ids[j]==id)
      {
        foundit = true;
        break;
      }
    if (foundit) continue;
    if (ninter==size)
    {
      size += 50;
      int* tmp = new int[size];
      for (int j=0; j<=ninter; ++j)
        tmp[j] = ids[j];
      delete [] ids;
      ids   = tmp;
      *iids = tmp;
    }
    ids[ninter] = id;
    ++ninter;
  }
  return (ninter);
}

/*----------------------------------------------------------------------*
 |  find dlines for given interface id                       m.gee 06/05|
 *----------------------------------------------------------------------*/
bool cca_mrtr_2D_finddlines(int id, DLINE** dline1, DLINE** dline2, DESIGN* design)
{
  // find the first dline of this pair
  for (int i=0; i<design->ndline; ++i)
  {
    *dline1 = NULL;
    *dline2 = NULL;
    if (design->dline[i].hasmortar==0) continue;
    if (design->dline[i].mrtr_Id != id) continue;
    
    // the first design line of this interface
    *dline1 = &(design->dline[i]);
    
    // find the second line of this pair
    for (int j=0; j<design->ndline; ++j)
    {
      if (design->dline[j].hasmortar==0) continue;
      if (design->dline[j].mrtr_Id != id) continue;
      if (&(design->dline[j]) == *dline1) continue;
      
      // the second design line of this pair
      *dline2 = &(design->dline[j]);
      break;
    }
    if (!(*dline1) || !(*dline2)) dserror("Cannot find dlines of mortar interface");
    else break;
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  find glines for given dline                              m.gee 06/05|
 *----------------------------------------------------------------------*/
int  cca_mrtr_2D_find_glines_on_dline(GLINE*** gline, DLINE* dline, 
                                      DISCRET* actdis)
{
  int size   = 50;
  int ngline = 0;
  *gline = new (GLINE*)[size];
  for (int i=0; i<actdis->ngline; ++i)
  {
    if (actdis->gline[i].dline != dline) continue;
    if (ngline==size)
    {
      size += 50;
      GLINE** tmp = new (GLINE*)[size];
      for (int j=0; j<ngline; ++j)
        tmp[j] = (*gline)[j];
      delete [] (*gline);
      *gline = tmp;
    }
    (*gline)[ngline] = &(actdis->gline[i]);
    ++ngline;
  }
  return (ngline);
}

/*----------------------------------------------------------------------*
 |  prepare data of a gline                                  m.gee 06/05|
 *----------------------------------------------------------------------*/
int cca_mrtr_2D_prepare_gline_data(GLINE* gline, int** nodeIds,
                                   MRTR::Segment::SegmentType* typ)
{
  // do typ of segment
  if (gline->ngnode==2)      *typ = MRTR::Segment::seg_Linear1D;
  else if (gline->ngnode==3) *typ = MRTR::Segment::seg_Quadratic1D;
  else cout << "***WRN*** gline " << gline->Id << " does not have 2 or 3 nodes\n";
  
  // The nodeIds of the segment have to be given in math positive order.
  // This is important as the mortar code will assume the outward normal
  // going through each segment in math positive direction 
  ELEMENT* ele = gline->gsurf[0]->element;
  int eleids[100];
  int ids[100];
  int nodelocalids[100];
  int count=0;

  (*nodeIds) = new int [gline->ngnode];
  
  for (int i=0; i<ele->numnp; ++i)
    eleids[i] = ele->node[i]->gnode->Id;
  for (int j=0; j<gline->ngnode; ++j)
    ids[j] = gline->gnode[j]->Id;
  
  for (int i=0; i<ele->numnp; ++i)
  {
    for (int j=0; j<gline->ngnode; ++j)
      if (eleids[i] == ids[j])
      {
        nodelocalids[count] = i;
        (*nodeIds)[count] = ids[j];
        ++count;
        break;
      }
  }
  if (count != gline->ngnode)
    cout << "***ERR*** Something wrong in finding nodeids for gline\n";

  // we have to check the nodelocalids for the combination 0 3
  // because in this case the order is actually 3 0
  if (nodelocalids[0]==0 && nodelocalids[1]==3)
  {
    int tmp = (*nodeIds)[0];
    (*nodeIds)[0] = (*nodeIds)[1];
    (*nodeIds)[1] = tmp;
  }
  return (count);
}

/*----------------------------------------------------------------------*
 |  find gnodes for given dline                              m.gee 06/05|
 *----------------------------------------------------------------------*/
int cca_mrtr_2D_find_gnodes_on_dline(GNODE*** gnode, DLINE* dline, DISCRET* actdis)
{
  int size   = 50;
  int ngnode = 0;
  *gnode = new (GNODE*)[size];

  // find the gnodes that are on this line
  for (int i=0; i<actdis->ngnode; ++i)
  {
    if (actdis->gnode[i].d.dline != dline) continue;
    if (ngnode==size)
    {
      size += 50;
      GNODE** tmp = new (GNODE*)[size];
      for (int j=0; j<ngnode; ++j)
        tmp[j] = (*gnode)[j];
      delete [] (*gnode);
      *gnode = tmp;
    }
    (*gnode)[ngnode] = &(actdis->gnode[i]);
    ++ngnode;
  }
  
  // add the gnodes that are on the dnodes of this line
  if (dline->ndnode != 2) 
    cout << "***WRN*** cca_mrtr_2D_find_gnodes_on_dline:\n"
         << "***WRN*** more then 2 dnodes on this dline\n";
  for (int i=0; i<dline->ndnode; ++i)
  {
    DNODE* dnode = dline->dnode[i];
    for (int j=0; j<actdis->ngnode; ++j)
    {
      if (actdis->gnode[j].d.dnode != dnode) continue;
      if (ngnode==size)
      {
        size += 50;
        GNODE** tmp = new (GNODE*)[size];
        for (int j=0; j<ngnode; ++j)
          tmp[j] = (*gnode)[j];
        delete [] (*gnode);
        *gnode = tmp;
      }
      (*gnode)[ngnode] = &(actdis->gnode[j]);
      ++ngnode;
      break;
    }
  }
    
  return (ngnode);
}





#endif // MORTAR
#endif // TRILINOS_PACKAGE
