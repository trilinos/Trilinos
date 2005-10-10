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
 |  counter number of 3D mortar interfaces                   m.gee 10/05|
 *----------------------------------------------------------------------*/
int cca_mrtr_3D_numinterfaces(int **iids, DESIGN* design)
{
  int size   = 50;
  int ninter = 0;
  *iids = new int[size];
  int* ids = *iids;
  for (int i=0; i<design->ndsurf; ++i)
  {
    if (design->dsurf[i].hasmortar==0) continue;
    int id = design->dsurf[i].mrtr_Id;
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
 |  find dsurfs for given interface id                       m.gee 10/05|
 *----------------------------------------------------------------------*/
bool cca_mrtr_3D_finddsurfs(int id, DSURF** dsurf1, DSURF** dsurf2, DESIGN* design)
{
  // find the first dline of this pair
  for (int i=0; i<design->ndsurf; ++i)
  {
    *dsurf1 = NULL;
    *dsurf2 = NULL;
    if (design->dsurf[i].hasmortar==0) continue;
    if (design->dsurf[i].mrtr_Id != id) continue;
    
    // the first design line of this interface
    *dsurf1 = &(design->dsurf[i]);
    
    // find the second line of this pair
    for (int j=0; j<design->ndsurf; ++j)
    {
      if (design->dsurf[j].hasmortar==0) continue;
      if (design->dsurf[j].mrtr_Id != id) continue;
      if (&(design->dsurf[j]) == *dsurf1) continue;
      
      // the second design line of this pair
      *dsurf2 = &(design->dsurf[j]);
      break;
    }
    if (!(*dsurf2) || !(*dsurf2)) dserror("Cannot find dlines of mortar interface");
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
 |  find gsurfs for given dsurf                              m.gee 10/05|
 *----------------------------------------------------------------------*/
int  cca_mrtr_3D_find_gsurfs_on_dsurf(GSURF*** gsurf, DSURF* dsurf, 
                                      DISCRET* actdis)
{
  int size   = 50;
  int ngsurf = 0;
  *gsurf = new (GSURF*)[size];
  for (int i=0; i<actdis->ngsurf; ++i)
  {
    if (actdis->gsurf[i].dsurf != dsurf) continue;
    if (ngsurf==size)
    {
      size += 50;
      GSURF** tmp = new (GSURF*)[size];
      for (int j=0; j<ngsurf; ++j)
        tmp[j] = (*gsurf)[j];
      delete [] (*gsurf);
      *gsurf = tmp;
    }
    (*gsurf)[ngsurf] = &(actdis->gsurf[i]);
    ++ngsurf;
  }
  return (ngsurf);
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
 |  prepare data of a gsurf                                  m.gee 10/05|
 *----------------------------------------------------------------------*/
static bool findnode(int ngnode, int* gsurfids, int eleid);
int cca_mrtr_3D_prepare_gsurf_data(GSURF* gsurf, int** nodeIds,
                                   MRTR::Segment::SegmentType* typ)
{
  // do typ of segment
  if (gsurf->ngnode==4)      *typ = MRTR::Segment::seg_BiLinearQuad;
  else if (gsurf->ngnode==3) *typ = MRTR::Segment::seg_BiLinearTri;
  else cout << "***WRN*** gsurf " << gsurf->Id << " does not have 4 or 3 nodes\n";
  
  // The nodeIds of the segment have to be given in math positive order.
  // This is important as the mortar code will assume the outward normal
  // going through each segment in math positive direction 
  ELEMENT* ele = gsurf->gvol[0]->element;

  int eleids[27];
  for (int i=0; i<ele->numnp; ++i)
    eleids[i] = ele->node[i]->Id;
    
  int gsurfids[9];
  for (int i=0; i<gsurf->ngnode; ++i)
    gsurfids[i] = gsurf->gnode[i]->node->Id;

  bool foundit = false;
  *nodeIds = new int[4];

  //--------------------------------- front surface
  // node 0 is on gsurf
  if (findnode(gsurf->ngnode,gsurfids,eleids[0]))
  {
    // node 1 is on gsurf
    if (findnode(gsurf->ngnode,gsurfids,eleids[1]))
    {
      // node 4 is on gsurf
      if (findnode(gsurf->ngnode,gsurfids,eleids[5]))
      {
        // node 5 is on gsurf
        if (findnode(gsurf->ngnode,gsurfids,eleids[4]))
        {
          foundit = true;
          (*nodeIds)[0] = eleids[0];
          (*nodeIds)[1] = eleids[1];
          (*nodeIds)[2] = eleids[5];
          (*nodeIds)[3] = eleids[4];
          return (1);
        }
      }
    }
  }
  
  //--------------------------------- bottom surface
  // node 0 is on gsurf
  if (findnode(gsurf->ngnode,gsurfids,eleids[0]))
  {
    // node 1 is on gsurf
    if (findnode(gsurf->ngnode,gsurfids,eleids[3]))
    {
      // node 4 is on gsurf
      if (findnode(gsurf->ngnode,gsurfids,eleids[2]))
      {
        // node 5 is on gsurf
        if (findnode(gsurf->ngnode,gsurfids,eleids[1]))
        {
          foundit = true;
          (*nodeIds)[0] = eleids[0];
          (*nodeIds)[1] = eleids[3];
          (*nodeIds)[2] = eleids[2];
          (*nodeIds)[3] = eleids[1];
          return (1);
        }
      }
    }
  }
  
  //--------------------------------- right surface
  // node 0 is on gsurf
  if (findnode(gsurf->ngnode,gsurfids,eleids[1]))
  {
    // node 1 is on gsurf
    if (findnode(gsurf->ngnode,gsurfids,eleids[2]))
    {
      // node 4 is on gsurf
      if (findnode(gsurf->ngnode,gsurfids,eleids[6]))
      {
        // node 5 is on gsurf
        if (findnode(gsurf->ngnode,gsurfids,eleids[5]))
        {
          foundit = true;
          (*nodeIds)[0] = eleids[1];
          (*nodeIds)[1] = eleids[2];
          (*nodeIds)[2] = eleids[6];
          (*nodeIds)[3] = eleids[5];
          return (1);
        }
      }
    }
  }
  
  //--------------------------------- left surface
  // node 0 is on gsurf
  if (findnode(gsurf->ngnode,gsurfids,eleids[0]))
  {
    // node 1 is on gsurf
    if (findnode(gsurf->ngnode,gsurfids,eleids[4]))
    {
      // node 4 is on gsurf
      if (findnode(gsurf->ngnode,gsurfids,eleids[7]))
      {
        // node 5 is on gsurf
        if (findnode(gsurf->ngnode,gsurfids,eleids[3]))
        {
          foundit = true;
          (*nodeIds)[0] = eleids[0];
          (*nodeIds)[1] = eleids[4];
          (*nodeIds)[2] = eleids[7];
          (*nodeIds)[3] = eleids[3];
          return (1);
        }
      }
    }
  }
  
  //--------------------------------- back surface
  // node 0 is on gsurf
  if (findnode(gsurf->ngnode,gsurfids,eleids[2]))
  {
    // node 1 is on gsurf
    if (findnode(gsurf->ngnode,gsurfids,eleids[3]))
    {
      // node 4 is on gsurf
      if (findnode(gsurf->ngnode,gsurfids,eleids[7]))
      {
        // node 5 is on gsurf
        if (findnode(gsurf->ngnode,gsurfids,eleids[6]))
        {
          foundit = true;
          (*nodeIds)[0] = eleids[2];
          (*nodeIds)[1] = eleids[3];
          (*nodeIds)[2] = eleids[7];
          (*nodeIds)[3] = eleids[6];
          return (1);
        }
      }
    }
  }
  
  //--------------------------------- top surface
  // node 0 is on gsurf
  if (findnode(gsurf->ngnode,gsurfids,eleids[4]))
  {
    // node 1 is on gsurf
    if (findnode(gsurf->ngnode,gsurfids,eleids[5]))
    {
      // node 4 is on gsurf
      if (findnode(gsurf->ngnode,gsurfids,eleids[6]))
      {
        // node 5 is on gsurf
        if (findnode(gsurf->ngnode,gsurfids,eleids[7]))
        {
          foundit = true;
          (*nodeIds)[0] = eleids[4];
          (*nodeIds)[1] = eleids[5];
          (*nodeIds)[2] = eleids[6];
          (*nodeIds)[3] = eleids[7];
          return (1);
        }
      }
    }
  }
  return (0);
}


bool findnode(int ngnode, int* gsurfids, int eleid)
{
  for (int i=0; i<ngnode; ++i)
    if (gsurfids[i]==eleid)
      return true;
  return false;
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


/*----------------------------------------------------------------------*
 |  find gnodes for given dsurf                              m.gee 10/05|
 *----------------------------------------------------------------------*/
int cca_mrtr_3D_find_gnodes_on_dsurf(GNODE*** gnode, DSURF* dsurf, DISCRET* actdis)
{
  int size   = 50;
  int ngnode = 0;
  *gnode = new (GNODE*)[size];

  // find the gnodes that are on this line
  for (int i=0; i<actdis->ngnode; ++i)
  {
    if (actdis->gnode[i].d.dsurf != dsurf) continue;
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
  
  // add the gnodes that are on the dlines of this dsurf
  for (int k=0; k<dsurf->ndline; ++k)
  {
    DLINE* actdline = dsurf->dline[k];
    for (int i=0; i<actdis->ngnode; ++i)
    {
      if (actdis->gnode[i].d.dline != actdline) continue;
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
  }
  
  // add the gnodes that are on dnodes of dlines of dsurf
  // build stack of dnodes
  DNODE* dnodes[100];
  int ndnode=0;
  for (int k=0; k<dsurf->ndline; ++k)
  {
    DLINE* actdline = dsurf->dline[k];
    for (int i=0; i<actdline->ndnode; ++i)
    {
      DNODE* actdnode = actdline->dnode[i];
      // look whether it's on stack already
      bool foundit = false;
      for (int j=0; j<ndnode; ++j)
        if (dnodes[j]==actdnode)
        {
          foundit = true;
          break;
        }
        if (foundit) continue;
        dnodes[ndnode] = actdnode;
        ++ndnode;
        if (ndnode==100) dserror("Overflow appears");
    }
  }
  // add gnodes on dnodes
  for (int k=0; k<ndnode; ++k)
  {
    DNODE* actdnode = dnodes[k];
    for (int i=0; i<actdis->ngnode; ++i)
    {
      if (actdis->gnode[i].d.dnode != actdnode) continue;
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
  }
  
    
  return (ngnode);
}





#endif // MORTAR
#endif // TRILINOS_PACKAGE
