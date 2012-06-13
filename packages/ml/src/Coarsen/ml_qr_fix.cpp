/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */        
/* ******************************************************************** */
#include "ml_include.h"
#include "ml_qr_fix.h"
#include "ml_qr_fix.hpp"
#include <stdio.h>

/* -mb: make this global var for now - decide where to put in ML */
static ML_qr_fix *QRFixStructure = NULL;

int ML_qr_fix_Create(const int nCoarseNod, int const nullspaceDim)
{
   QRFixStructure = new ML_qr_fix();
   QRFixStructure->numAggsWithDeadDofs = 0;
   QRFixStructure->coarseDOFState = std::vector<std::vector<bool> >(nCoarseNod, std::vector<bool>(nullspaceDim,true));
   QRFixStructure->aggTooSmall = std::vector<bool>(nCoarseNod,false);

   return(0);
}

int ML_qr_fix_Destroy(void)
{
  if (QRFixStructure == NULL)
    return(0);
  delete QRFixStructure;
  QRFixStructure = NULL;
  return(0);
}

int ML_qr_fix_NumDeadNodDof(void) { 
    if (QRFixStructure == NULL) return 0;
    return QRFixStructure->numAggsWithDeadDofs; 
}

void ML_qr_fix_setNumDeadNod(int num)
{
     if (QRFixStructure == NULL) return;
     QRFixStructure->numAggsWithDeadDofs = num;
}

//   nodeNum -- local aggregate number
//   aggSize -- number of DOFs in this aggregate
//   nullspaceDim -- dimension of nullspace
void ML_qr_fix_markDOFsAsDead(const int nodeNum, const int aggSize, const int nullspaceDim)
{
  if (QRFixStructure == NULL) return;
  // Note:  false = alive
  //        true  = dead
  for (int dof=0; dof < aggSize; ++dof) ((QRFixStructure->coarseDOFState)[nodeNum])[dof] = false;
  for (int dof=aggSize; dof < nullspaceDim; ++dof) ((QRFixStructure->coarseDOFState)[nodeNum])[dof] = true;
  if (nullspaceDim > aggSize)
    QRFixStructure->aggTooSmall[nodeNum] = true;
}

int ML_qr_fix_isDOFDead(const int nodeNum, const int coarseDOF)
{
  if (QRFixStructure == NULL) return 0;
  if ( ((QRFixStructure->coarseDOFState)[nodeNum])[coarseDOF] == true ) return 1;
  else                                                           return 0;
}

int ML_qr_fix_nodeHasDeadDOFs(const int nodeNum)
{
  if (QRFixStructure == NULL) return 0;
  if ( QRFixStructure->aggTooSmall[nodeNum] == true ) return 1;
  else                                             return 0;
}
