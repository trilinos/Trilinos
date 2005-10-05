#ifdef MORTAR
#ifdef TRILINOS_PACKAGE
/*----------------------------------------------------------------------*
 |  routine to control mortar interfaces                     m.gee 06/05|
 *----------------------------------------------------------------------*/
#ifndef CCA_MORTAR_H
#define CCA_MORTAR_H

#include "mrtr_segment.H"

extern "C" 
{
#include "../headers/standardtypes.h"
#ifdef PARALLEL
int create_mortar(FIELD *actfield, PARTITION *actpart,
                  int disnum, DESIGN *design, int init, MPI_Comm mpicomm);
#else
int create_mortar(FIELD *actfield, PARTITION *actpart,
                  int disnum, DESIGN *design, int init);
#endif 
int compute_mortar(SPARSE_TYP* arraytyp, SPARSE_ARRAY* array, DESIGN *design);
}

int  cca_mrtr_2D_numinterfaces(int **iids, DESIGN* design);
bool cca_mrtr_2D_finddlines(int id, DLINE** dline1, DLINE** dline2, DESIGN* design);
int  cca_mrtr_2D_find_glines_on_dline(GLINE*** gline, DLINE* dline, DISCRET* actdis);
int  cca_mrtr_2D_find_gnodes_on_dline(GNODE*** gnode, DLINE* dline, DISCRET* actdis);
int cca_mrtr_2D_prepare_gline_data(GLINE* gline, int** nodeIds,
                                   MRTR::Segment::SegmentType* typ);

int  cca_mrtr_3D_numinterfaces(int **iids, DESIGN* design);
bool cca_mrtr_3D_finddsurfs(int id, DSURF** dsurf1, DSURF** dsurf2, DESIGN* design);
int  cca_mrtr_3D_find_gsurfs_on_dsurf(GSURF*** gsurf, DSURF* dsurf, 
                                      DISCRET* actdis);
int cca_mrtr_3D_prepare_gsurf_data(GSURF* gline, int** nodeIds,
                                   MRTR::Segment::SegmentType* typ);
int cca_mrtr_3D_find_gnodes_on_dsurf(GNODE*** gnode, DSURF* dsurf, DISCRET* actdis);
#endif
#endif // TRILINOS_PACKAGE
#endif // MORTAR
