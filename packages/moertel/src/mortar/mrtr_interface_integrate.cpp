/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include <ctime>
#include <vector>

#include "mrtr_interface.H"
#include "mrtr_utils.H"
#include "mrtr_pnode.H"
#include "mrtr_segment.H"
#include "mrtr_integrator.H"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  assemble values from integration                                    |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Mortar_Assemble(Epetra_CrsMatrix& D, 
                                       Epetra_CrsMatrix& M)
{ 

  //-------------------------------------------------------------------
  // interface needs to be complete
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Assemble:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // send all procs not member of this interface's intra-comm out of here
  if (!lComm()) return true;

  //-------------------------------------------------------------------
  // interface needs to have a mortar side assigned
  if (MortarSide()==-1)
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Assemble:\n"
           << "***ERR*** mortar side was not assigned on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // interface need to be integrated
  if (!IsIntegrated())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Assemble:\n"
           << "***ERR*** interface " << Id_ << " not integrated\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // call assembly of 2D and 3D problems
  return Assemble_3D(D,M);
}

#if 0 // old version
/*----------------------------------------------------------------------*
 |  make mortar integration of this interface (2D problem)              |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Mortar_Integrate(Epetra_CrsMatrix& D, 
                                          Epetra_CrsMatrix& M)
{ 
  bool ok = false;
  
  //-------------------------------------------------------------------
  // time this process
  Epetra_Time time(*lComm());
  time.ResetStartTime();

  //-------------------------------------------------------------------
  if (!IsOneDimensional())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** This is not a 2D problem, we're in the wrong method here!!!\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // interface needs to be complete
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // send all procs not member of this interface's intra-comm out of here
  if (!lComm()) return true;

  //-------------------------------------------------------------------
  // interface needs to have a mortar side assigned
  if (MortarSide()==-1)
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** mortar side was not assigned on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // interface segments need to have at least one function on the mortar side
  // and two functions on the slave side
  int mside = MortarSide();
  int sside = OtherSide(mside);
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator scurr;
  for (scurr=seg_[mside].begin(); scurr!=seg_[mside].end(); ++scurr)
    if (scurr->second->Nfunctions() < 1)
    {
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", mortar side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 1 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  for (scurr=seg_[sside].begin(); scurr!=seg_[sside].end(); ++scurr)
    if (scurr->second->Nfunctions() < 2)
    {
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", slave side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 2 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    
  //-------------------------------------------------------------------
  // do the integration of the master and slave side
  ok = Integrate_2D(M,D);
  if (!ok) return false;

  //-------------------------------------------------------------------
  // set the flag that this interface has been successfully integrated
  isIntegrated_ = true;
  
  //-------------------------------------------------------------------
  // time this process
  if (OutLevel()>5)
  {
    cout << "MOERTEL::Interface " << Id() << ": Integration on proc " << gComm().MyPID()
         << " finished in " << time.ElapsedTime() << " sec\n"; fflush(stdout);
  }
  return true;
}
#endif

/*----------------------------------------------------------------------*
 |  make mortar integration of this interface (2D problem)              |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Mortar_Integrate_2D(
                           Teuchos::RCP<Teuchos::ParameterList> intparams)
{ 
  bool ok = false;
  intparams_ = intparams;
  
  //-------------------------------------------------------------------
  // time this process
  Epetra_Time time(*lComm());
  time.ResetStartTime();

  //-------------------------------------------------------------------
  if (!IsOneDimensional())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** This is not a 2D problem, we're in the wrong method here!!!\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // interface needs to be complete
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // send all procs not member of this interface's intra-comm out of here
  if (!lComm()) return true;

  //-------------------------------------------------------------------
  // interface needs to have a mortar side assigned
  if (MortarSide()==-1)
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** mortar side was not assigned on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // interface segments need to have at least one function on the mortar side
  // and two functions on the slave side
  int mside = MortarSide();
  int sside = OtherSide(mside);
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator scurr;
  for (scurr=seg_[mside].begin(); scurr!=seg_[mside].end(); ++scurr)
    if (scurr->second->Nfunctions() < 1)
    {
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", mortar side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 1 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  for (scurr=seg_[sside].begin(); scurr!=seg_[sside].end(); ++scurr)
    if (scurr->second->Nfunctions() < 2)
    {
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", slave side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 2 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    
  //-------------------------------------------------------------------
  // do the integration of the master and slave side
  ok = Integrate_2D();
  if (!ok) return false;

  //-------------------------------------------------------------------
  // set the flag that this interface has been successfully integrated
  isIntegrated_ = true;
  
  //-------------------------------------------------------------------
  // time this process
  if (OutLevel()>5)
  {
    cout << "MOERTEL::Interface " << Id() << ": Integration on proc " << gComm().MyPID()
         << " finished in " << time.ElapsedTime() << " sec\n"; fflush(stdout);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |  make mortar integration of master/slave side in 2D (1D interface)   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_2D()
{ 
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Integrate_2D:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  
  // loop over all segments of slave side
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator scurr;
  for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)
  {
    // the segment to be integrated
	Teuchos::RCP<MOERTEL::Segment> actsseg = scurr->second;

#if 0
    cout << "\nActive sseg id " << actsseg->Id() << "\n\n";
#endif

    // check whether I own at least one of the nodes on this slave segment
    int nnode = actsseg->Nnode();
    MOERTEL::Node** nodes = actsseg->Nodes();
    bool foundone = false;
    for (int i=0; i<nnode; ++i)
      if (NodePID(nodes[i]->Id()) == lComm()->MyPID())
      {
        foundone = true;
        break;
      }
    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // loop over all segments on the master side
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator mcurr;
    for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)    
    {
	  Teuchos::RCP<MOERTEL::Segment> actmseg = mcurr->second;
      
#if 0
    cout << "Active mseg id " << actmseg->Id() << endl;
#endif
      // if there is an overlap, integrate the pair
      // (whether there is an overlap or not will be checked inside)
      Integrate_2D_Section(*actsseg,*actmseg);
      
    } // for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)  
  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}


#if 0 // old version
/*----------------------------------------------------------------------*
 |  make mortar integration of master/slave side in 2D (1D interface)   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_2D(Epetra_CrsMatrix& M,
                                      Epetra_CrsMatrix& D)
{ 
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Integrate_2D:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  
  // loop over all segments of slave side
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator scurr;
  for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)
  {
    // the segment to be integrated
	Teuchos::RCP<MOERTEL::Segment> actsseg = scurr->second;

#if 0
    cout << "\nActive sseg id " << actsseg->Id() << "\n\n";
#endif

    // check whether I own at least one of the nodes on this slave segment
    int nnode = actsseg->Nnode();
    MOERTEL::Node** nodes = actsseg->Nodes();
    bool foundone = false;
    for (int i=0; i<nnode; ++i)
      if (NodePID(nodes[i]->Id()) == lComm()->MyPID())
      {
        foundone = true;
        break;
      }
    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // loop over all segments on the master side
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator mcurr;
    for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)    
    {
	  Teuchos::RCP<MOERTEL::Segment> actmseg = mcurr->second;
      
#if 0
    cout << "Active mseg id " << actmseg->Id() << endl;
#endif
      // if there is an overlap, integrate the pair
      // (whether there is an overlap or not will be checked inside)
      Integrate_2D_Section(*actsseg,*actmseg,M,D);
      
    } // for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)  
  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}
#endif

/*----------------------------------------------------------------------*
 | integrate the master/slave side's contribution from the overlap      |
 | of 2 segments (2D version) IF there is an overlap                    |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_2D_Section(MOERTEL::Segment& sseg, 
                                              MOERTEL::Segment& mseg)
{ 
  // if one of the segments is quadratic, we have to do something here
  if (sseg.Type()!=MOERTEL::Segment::seg_Linear1D || mseg.Type()!=MOERTEL::Segment::seg_Linear1D)
  {
	std::stringstream oss;
		oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
         << "***ERR*** Integration of other than linear segments not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw ReportError(oss);
  }

#if 0
  cout << "\n\nSlave Segment:\n";
  cout << sseg;
  cout << "Master Segment:\n";
  cout << mseg;
#endif

  // there is several cases on how these 2 segments can overlap
  // handle all of them, including the ones that they don't overlap 
  // at all
  
  // Do a coarse check to see if the segments are even close to each other.

  bool overlap = false;

  overlap = QuickOverlapTest_2D(mseg, sseg);

  if (!overlap)
    return true;

  // get slave and master's projections of the end points
  MOERTEL::Node** snodes = sseg.Nodes();
  MOERTEL::Node** mnodes = mseg.Nodes();
  
  // determine the overlap of the 2 segments if there is any
  MOERTEL::Projector projector(IsOneDimensional(),OutLevel());
  // project master nodes onto slave segment
  std::vector<double> mxi(mseg.Nnode());
  std::vector<double> mgap(mseg.Nnode());
  for (int i=0; i<mseg.Nnode(); ++i)
    projector.ProjectNodetoSegment_SegmentNormal(*mnodes[i],sseg,&mxi[i],mgap[i]);
  //cout << mxi[0] << " " << mxi[1] << endl;

  // project slave nodes onto master segment
  std::vector<double> sxi(sseg.Nnode());
  std::vector<double> sgap(sseg.Nnode());
  for (int i=0; i<sseg.Nnode(); ++i)
    projector.ProjectNodetoSegment_NodalNormal(*snodes[i],mseg,&sxi[i],sgap[i]);
  //cout << sxi[0] << " " << sxi[1] << endl;
  
  // Depending on mxi and sxi decide on the overlap
  bool snode0 = false;
  bool snode1 = false;
  bool mnode0 = false;
  bool mnode1 = false;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_spnode0 = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_spnode1 = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_mpnode0 = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_mpnode1 = Teuchos::null;
  double xi[2]; xi[0] = xi[1] = 0.0;
  if ( -1.05 <= mxi[0] && mxi[0] <= 1.05) 
  {
    mnode0 = true;
    xi[0] = mxi[0];
    is_mpnode0 = Teuchos::rcp(new MOERTEL::ProjectedNode(*mnodes[0],xi,&sseg));
	mnodes[0]->SetGap(mgap[0]);
  }
  if ( -1.05 <= mxi[1] && mxi[1] <= 1.05) 
  {
    mnode1 = true;
    xi[0] = mxi[1];
    is_mpnode1 = Teuchos::rcp(new MOERTEL::ProjectedNode(*mnodes[1],xi,&sseg));
	mnodes[1]->SetGap(mgap[1]);
  }
  if ( -1.05 <= sxi[0] && sxi[0] <= 1.05) 
  {
    snode0 = true;
    xi[0] = sxi[0];
    is_spnode0 = Teuchos::rcp(new MOERTEL::ProjectedNode(*snodes[0],xi,&mseg));
	snodes[0]->SetGap(sgap[0]);
  }
  if ( -1.05 <= sxi[1] && sxi[1] <= 1.05) 
  {
    snode1 = true;
    xi[0] = sxi[1];
    is_spnode1 = Teuchos::rcp(new MOERTEL::ProjectedNode(*snodes[1],xi,&mseg));
	snodes[1]->SetGap(sgap[1]);
  }
  //cout << mnode0 << "  " << mnode1 << "  " << snode0 << "  " << snode1 << endl;
  
  // Make decision upon overlap
  overlap = false;
  Teuchos::RCP<MOERTEL::ProjectedNode> nstart = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> nend   = Teuchos::null;
  double sxia=999.0,sxib=999.0;
  double mxia=999.0,mxib=999.0;
  
  // no overlap
  if (!snode0 && !snode1 && !mnode0 && !mnode1);
  // no overlap
  else if (snode0 && !snode1 && !mnode0 && !mnode1)
  {
    if (sxi[0]>-0.95)
      cout << "***WRN*** Significant overlap ignored\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  // no overlap
  else if (!snode0 && !snode1 && mnode0 && !mnode1)
  {
    if (mxi[0]>-0.95)
      cout << "MOERTEL: ***WRN*** Significant overlap ignored\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  else if (!snode0 && !snode1 && !mnode0 && mnode1)
  {
    if (mxi[1]<0.95)
      cout << "MOERTEL: ***WRN*** Significant overlap ignored\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  else if (!snode0 && snode1 && !mnode0 && !mnode1)
  {
    if (sxi[1]<0.95)
      cout << "***WRN*** Significant overlap ignored\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  else if (mnode0 && mnode1)
  {
    overlap = true;
    nstart = is_mpnode0;
    nend   = is_mpnode1;
    sxia = nend->Xi()[0];
    sxib = nstart->Xi()[0];
    mxia = -1.0;
    mxib = 1.0;
  }
  else if (snode0 && snode1)
  {
    overlap = true;
    nstart = is_spnode0;
    nend   = is_spnode1;
    sxia = -1.0;
    sxib =  1.0;
    mxia = nend->Xi()[0];
    mxib = nstart->Xi()[0];
  }
  else if (snode0 && !snode1 && mnode0 && !mnode1)
  {
    overlap = true;
    nstart = is_spnode0;
    nend   = is_mpnode0;
    sxia = -1.0;
    sxib = nend->Xi()[0];
    mxia = -1.0;
    mxib = nstart->Xi()[0];
  }
  else if (snode1 && !snode0 && mnode1 && !mnode0)
  {
    overlap = true;
    nstart = is_mpnode1;
    nend   = is_spnode1;
    sxia = nstart->Xi()[0];
    sxib = 1.0;
    mxia = nend->Xi()[0];
    mxib = 1.0;
  }
  else
  {
	
	std::stringstream oss;
		oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
         << "***ERR*** Unknown overlap case found\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw ReportError(oss);
  }
  if (!overlap)
    return true;

#if 0  
  cout << "slave  xi range " << sxia << " - " << sxib << endl;
  cout << "master xi range " << mxia << " - " << mxib << endl;
#endif

  // create an integrator instance of some given order
  MOERTEL::Integrator integrator(5,IsOneDimensional(),OutLevel());

  // do the integration of the master side
  Epetra_SerialDenseMatrix* Mdense = 
                            integrator.Integrate(sseg,sxia,sxib,mseg,mxia,mxib);
  
  // do the integration of the slave side
  Epetra_SerialDenseMatrix* Ddense = integrator.Integrate(sseg,sxia,sxib);

  // Assemble contributions Mdense into nodes (scalar only)
  integrator.Assemble(*this,sseg,mseg,*Mdense);

  // Assemble contributions Ddense into nodes (scalar only)
  integrator.Assemble(*this,sseg,*Ddense);

#if 0 // modification for curved interfaces from paper by B. Wohlmuth

  if (sseg.Type() == MOERTEL::Segment::seg_Linear1D && 
      mseg.Type() == MOERTEL::Segment::seg_Linear1D)
  if (sseg.FunctionType(1) == MOERTEL::Function::func_DualLinear1D)
  if (mnodes[0]->Ndof() == mnodes[1]->Ndof() &&
      mnodes[0]->Ndof() == 2)
  {
    Epetra_SerialDenseMatrix* Mmod = NULL;
    
    // get the normal at slave nodes
    const double* n0 = snodes[0]->N();
    const double* n1 = snodes[1]->N();

    // build the tangential orthogonal to the normal
    double t[2][2];
    t[0][0] = -n0[1]; t[1][0] = -n1[1];
    t[0][1] =  n0[0]; t[1][1] =  n1[0];
    double n[2][2];
    n[0][0] =  n0[0]; n[1][0] =  n1[0]; 
    n[0][1] =  n0[1]; n[1][1] =  n1[1]; 
    
    // build delta values of normal and tangential
    double dn[2]; double dt[2];
    dn[0] = n0[0] - n1[0];  
    dn[1] = n0[1] - n1[1];  
    dt[0] = t[0][0] - t[1][0];
    dt[1] = t[0][1] - t[1][1];
    
    // build norm of dn. If it's zero, don't do anything
    bool doit = true;
//    double delta = dn[0]*dn[0]+dn[1]*dn[1];
//    if (abs(delta)>1.0e-11) doit = true;

    if (doit)
    {
      // do the integration of the modification of the master side
      // integral ( -0.5 * psi_12 * phi_k ) k=1,...,nnode_master 
      Epetra_SerialDenseMatrix* Mmod_scalar =
                        integrator.Integrate_2D_Mmod(sseg,sxia,sxib,mseg,mxia,mxib);

      // create an Epetra_SerialDenseMatrix of dimension (nsnode x nlmdof , nmnode x nmdof)
      int nsnode = sseg.Nnode();
      int nsdof  = snodes[0]->Ndof();
      int nmnode = mseg.Nnode();
      int nmdof  = mnodes[0]->Ndof();
      Mmod =  new Epetra_SerialDenseMatrix(nsnode*nsdof,nmnode*nmdof);

      // add modification values to Mmod
      for (int snode=0; snode<nsnode; ++snode)
        for (int sdof=0; sdof<nsdof; ++sdof)
        {
          double nt[2];
          nt[0] = n[snode][sdof] * dn[0] + t[snode][sdof] * dt[0];
          nt[1] = n[snode][sdof] * dn[1] + t[snode][sdof] * dt[1];
          for (int mnode=0; mnode<nmnode; ++mnode)
            for (int mdof=0; mdof<nmdof; ++mdof)
            {
              double val = nt[mdof] * (*Mmod_scalar)(mnode,0);
              (*Mmod)(snode*nsdof+sdof,mnode*nmdof+mdof) = val;
            }
        } // for (int sdof=0; sdof<nsdof; ++sdof)

#if 0  // verification of the expression by expressions given in paper
      Epetra_SerialDenseMatrix* Mmod2 = new Epetra_SerialDenseMatrix(nsnode*nsdof,nmnode*nmdof);
      // n1 dot n2
      double n1n2 = 0.0;
      for (int i=0; i<2; ++i) n1n2 += n[0][i]*n[1][i];
      // third row of n1 x n2
      double n1xn2 = n[0][0]*n[1][1] - n[0][1]*n[1][0];
      
      // slave 0 sdof 0 master 0 mdof 0 
      (*Mmod2)(0,0) = (*Mmod_scalar)(0,0) * (1.0-n1n2);
      // slave 0 sdof 0 master 0 mdof 1
      (*Mmod2)(0,1) = - (*Mmod_scalar)(0,0) * n1xn2;
      // slave 0 sdof 0 master 1 mdof 0
      (*Mmod2)(0,2) = (*Mmod_scalar)(1,0) * (1.0-n1n2);
      // slave 0 sdof 0 master 1 mdof 1
      (*Mmod2)(0,3) = - (*Mmod_scalar)(1,0) * n1xn2;
      // slave 0 sdof 1 master 0 mdof 0 
      (*Mmod2)(1,0) = (*Mmod_scalar)(0,0) * n1xn2;
      // slave 0 sdof 1 master 0 mdof 1
      (*Mmod2)(1,1) = (*Mmod_scalar)(0,0) * (1.0-n1n2);
      // slave 0 sdof 1 master 1 mdof 0
      (*Mmod2)(1,2) = (*Mmod_scalar)(1,0) * n1xn2;
      // slave 0 sdof 1 master 1 mdof 1
      (*Mmod2)(1,3) = (*Mmod_scalar)(1,0) * (1.0-n1n2);
      // slave 1 sdof 0 master 0 mdof 0
      (*Mmod2)(2,0) = (*Mmod_scalar)(0,0) * (n1n2-1.0);
      // slave 1 sdof 0 master 0 mdof 1
      (*Mmod2)(2,1) = - (*Mmod_scalar)(0,0) * n1xn2;
      // slave 1 sdof 0 master 1 mdof 0
      (*Mmod2)(2,2) = (*Mmod_scalar)(1,0) * (n1n2-1.0);
      // slave 1 sdof 0 master 1 mdof 1
      (*Mmod2)(2,3) = - (*Mmod_scalar)(1,0) * n1xn2;
      // slave 1 sdof 1 master 0 mdof 0
      (*Mmod2)(3,0) = (*Mmod_scalar)(0,0) * n1xn2;
      // slave 1 sdof 1 master 0 mdof 1
      (*Mmod2)(3,1) = (*Mmod_scalar)(0,0) * (n1n2-1.0);
      // slave 1 sdof 1 master 1 mdof 0
      (*Mmod2)(3,2) = (*Mmod_scalar)(1,0) * n1xn2;
      // slave 1 sdof 1 master 1 mdof 1
      (*Mmod2)(3,3) = (*Mmod_scalar)(1,0) * (n1n2-1.0);
      //cout << *Mmod2;
      //delete Mmod2; Mmod2 = NULL;
#endif

      //  assemble -Mmod into M
      integrator.Assemble_2D_Mod(*this,sseg,mseg,*Mmod);
      
      // tidy up 
      if (Mmod)        delete Mmod;        Mmod = NULL;
      if (Mmod_scalar) delete Mmod_scalar; Mmod_scalar = NULL;
    } // if (doit)
  } // if modification
  

#endif

  return true;
}

bool MOERTEL::Interface::QuickOverlapTest_2D(MOERTEL::Segment& sseg, MOERTEL::Segment& mseg)
{

  MOERTEL::Node** snode = sseg.Nodes();
  MOERTEL::Node** mnode = mseg.Nodes();
  const int nsnode = sseg.Nnode();
  const int nmnode = mseg.Nnode();

  double mcen[3], scen[3], mrad[3], srad[3], vec[3], mdiam, sdiam, length;

  mcen[0] = mcen[1] = mcen[2] = 0;
  scen[0] = scen[1] = scen[2] = 0;
  mdiam = sdiam = 0;

  for (int i=0; i<nmnode; ++i){
	mcen[0] += mnode[i]->X()[0];
	mcen[1] += mnode[i]->X()[1];
	mcen[2] += mnode[i]->X()[2];
  }
  mcen[0] /= (double)nmnode;
  mcen[1] /= (double)nmnode;
  mcen[2] /= (double)nmnode;

  for (int i=0; i<nsnode; ++i){
	scen[0] += snode[i]->X()[0];
	scen[1] += snode[i]->X()[1];
	scen[2] += snode[i]->X()[2];
  }
  scen[0] /= (double)nsnode;
  scen[1] /= (double)nsnode;
  scen[2] /= (double)nsnode;

  for (int i=0; i<nmnode; ++i){
	mrad[0] = mnode[i]->X()[0] - mcen[0];
	mrad[1] = mnode[i]->X()[1] - mcen[1];
	mrad[2] = mnode[i]->X()[2] - mcen[2];
	length = MOERTEL::length(mrad,3);
	if (mdiam < length) mdiam = length;
  }

  for (int i=0; i<nsnode; ++i){
	srad[0] = snode[i]->X()[0] - scen[0];
	srad[1] = snode[i]->X()[1] - scen[1];
	srad[2] = snode[i]->X()[2] - scen[2];
	length = MOERTEL::length(srad,3);
	if (sdiam < length) sdiam = length;
  }

  vec[0] = mcen[0] - scen[0];
  vec[1] = mcen[1] - scen[1];
  vec[2] = mcen[2] - scen[2];
  length = MOERTEL::length(vec,3);

  // GAH EPSILON - max distance between mseg and sseg for contact purposes
  
  double maxdia = 2.5;

  if (length > maxdia * (sdiam + mdiam)){

    // std::cerr << " test NOT passed\n";
    return false;
  }
  
  return true;

}

#if 0 // old version
/*----------------------------------------------------------------------*
 | integrate the master/slave side's contribution from the overlap      |
 | of 2 segments (2D version) IF there is an overlap                    |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_2D_Section(MOERTEL::Segment& sseg, 
                                              MOERTEL::Segment& mseg,
                                              Epetra_CrsMatrix& M,
                                              Epetra_CrsMatrix& D)
{ 
  // if one of the segments is quadratic, we have to do something here
  if (sseg.Type()!=MOERTEL::Segment::seg_Linear1D || mseg.Type()!=MOERTEL::Segment::seg_Linear1D)
  {
	std::stringstream oss;
    oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
         << "***ERR*** Integration of other then linear segments not yet implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    throw ReportError(oss);
  }

#if 0
  cout << "\n\nSlave Segment:\n";
  cout << sseg;
  cout << "Master Segment:\n";
  cout << mseg;
#endif

  // there is several cases on how these 2 segments can overlap
  // handle all of them, including the ones that they don't overlap 
  // at all
  
  // get slave and master's projections of the end points
  MOERTEL::Node** snodes = sseg.Nodes();
  MOERTEL::Node** mnodes = mseg.Nodes();
  
#if 0
  cout << "snodes[0]\n" << *snodes[0];
  cout << "snodes[1]\n" << *snodes[1];
  cout << "mnodes[0]\n" << *mnodes[0];
  cout << "mnodes[1]\n" << *mnodes[1];
#endif  

  bool snode0 = false;
  bool snode1 = false;
  bool mnode0 = false;
  bool mnode1 = false;
  int foundcase =  0;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_spnode0 = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_spnode1 = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_mpnode0 = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> is_mpnode1 = Teuchos::null;
  
  // projection along continous normal field results in projection points
  // that are unique
  if (GetProjectionType() == proj_continousnormalfield)
  {  
    if (snodes[0]->GetProjectedNode() != Teuchos::null)
      if (snodes[0]->GetProjectedNode()->Segment())
        if (snodes[0]->GetProjectedNode()->Segment()->Id() == mseg.Id())
        {
          snode0     = true;
          is_spnode0 = snodes[0]->GetProjectedNode();
        }
    if (snodes[1]->GetProjectedNode() != Teuchos::null)
      if (snodes[1]->GetProjectedNode()->Segment())
        if (snodes[1]->GetProjectedNode()->Segment()->Id() == mseg.Id())
        {
          snode1     = true;
          is_spnode1 = snodes[1]->GetProjectedNode(); 
        }
    if (mnodes[0]->GetProjectedNode() != Teuchos::null)
      if (mnodes[0]->GetProjectedNode()->Segment())
        if (mnodes[0]->GetProjectedNode()->Segment()->Id() == sseg.Id())
        {
          mnode0     = true;
          is_mpnode0 = mnodes[0]->GetProjectedNode();
        }
    if (mnodes[1]->GetProjectedNode() != Teuchos::null)
      if (mnodes[1]->GetProjectedNode()->Segment())
        if (mnodes[1]->GetProjectedNode()->Segment()->Id() == sseg.Id())
        {
          mnode1     = true;
          is_mpnode1 = mnodes[1]->GetProjectedNode();
        }
  }
  // projection orthogonal to some slave segment results in multiple projection
  // point for the slave side. Here, we pick the one that has been projected
  // orthogonal to the current slave segment
  else if (GetProjectionType() == proj_orthogonal)
  {
    int nspnode0;
	Teuchos::RCP<MOERTEL::ProjectedNode>* spnode0 = snodes[0]->GetProjectedNode(nspnode0);
    if (spnode0)
      for (int i=0; i<nspnode0; ++i)
        if (spnode0[i]->Segment())
          if (spnode0[i]->Segment()->Id() == mseg.Id())
            if (spnode0[i]->OrthoSegment() == sseg.Id())
            {
#if 0
              cout << " snode id: " << spnode0[i]->Id()
                   << " projects on mseg: " << mseg.Id()
                   << " orth to sseg: " << spnode0[i]->OrthoSegment() << endl;
#endif
              snode0     = true;
              is_spnode0 = spnode0[i];  
              break;
            }
    
    int nspnode1;
	Teuchos::RCP<MOERTEL::ProjectedNode>* spnode1 = snodes[1]->GetProjectedNode(nspnode1);
    if (spnode1)
      for (int i=0; i<nspnode1; ++i)
        if (spnode1[i]->Segment())
          if (spnode1[i]->Segment()->Id() == mseg.Id())
            if (spnode1[i]->OrthoSegment() == sseg.Id())
            {
#if 0
              cout << " snode id: " << spnode1[i]->Id()
                   << " projects on mseg: " << mseg.Id()
                   << " orth to sseg: " << spnode1[i]->OrthoSegment() << endl;
#endif
              snode1 = true;  
              is_spnode1 = spnode1[i];  
              break;
            }

    if (mnodes[0]->GetProjectedNode() != Teuchos::null)
      if (mnodes[0]->GetProjectedNode()->Segment())
        if (mnodes[0]->GetProjectedNode()->Segment()->Id() == sseg.Id())
        {
          mnode0     = true;
          is_mpnode0 = mnodes[0]->GetProjectedNode(); 
        }
    if (mnodes[1]->GetProjectedNode() != Teuchos::null)
      if (mnodes[1]->GetProjectedNode()->Segment())
        if (mnodes[1]->GetProjectedNode()->Segment()->Id() == sseg.Id())
        {
          mnode1 = true;
          is_mpnode1 = mnodes[1]->GetProjectedNode(); 
        }
  }
  
        
  Teuchos::RCP<MOERTEL::ProjectedNode> nstart = Teuchos::null;
  Teuchos::RCP<MOERTEL::ProjectedNode> nend   = Teuchos::null;

  // the xi range to integrate
  double sxia=999.0,sxib=999.0;
  double mxia=999.0,mxib=999.0;
  
  // case 1: snodes don't project into master element and
  //         mnodes don't project into slave element
  if (!snode0 && !snode1 && !mnode0 && !mnode1)
  {
    //cout << "Case 1: no overlap\n";
    ++foundcase;
  }
  
  // case 2: snode0 projects into master element
  //         snode1 not projects into master element
  //         mnodes not project into slave element
  // Note: this case is due to tolerance in projection
  if (snode0 && !snode1 && !mnode0 && !mnode1)
    ++foundcase;
  
  // case 3: mnode0 projects into slave element
  //         mnode1 not projects into slave element
  //         snodes don't project into master element
  // Note: this case is due to tolerance in projection
  if (!snode0 && !snode1 && mnode0 && !mnode1)
    ++foundcase;
  
  // case 4: mnode0 does not project into slave element
  //         mnode1 projects into slave element
  //         snodes don't project into master element
  // Note: this case is due to tolerance in projection
  if (!snode0 && !snode1 && !mnode0 && mnode1)
  {
    bool ok = false;
    // to do nothing, mnode1 has to project really low in the slave segment
    nstart = is_mpnode1;
    sxia = nstart->Xi()[0];
    if (sxia>0.95)
    {
      ++foundcase;
      nstart = Teuchos::null;
      nend   = Teuchos::null;
    }
    else
    {
      ok = true;
      sxib = 1.0;
      // for the range of the master element, we need to check whether
      // 1.) mnodes[0] projects into a neighbor of sseg
      // 2.) mnodes[0]'s projection hast to be low in xi
      nend = mnodes[0]->GetProjectedNode();
      if (nend != Teuchos::null) ok = true;
      else      ok = false;
      if (ok)
      {
        if (!nend->Segment()) ok = true;
        else
        {
          int nseg             = snodes[1]->Nseg();
          MOERTEL::Segment** segs = snodes[1]->Segments();
          int segid = nend->Segment()->Id();
          for (int i=0; i<nseg; ++i)
            if (segid==segs[i]->Id()) { ok = true; break;}
            else ok = false;
          if (ok)
          {
            double xi = nend->Xi()[0];
            if (xi>-0.95)
              ok = false;
          }
        }
      }
      if (ok)
      {
        mxia = -1.0;
        mxib = 1.0;
        ++foundcase;
      }
      else // do nothing?
      {
        nstart = Teuchos::null;
        nend   = Teuchos::null;
        ++foundcase;
      }
    }
  }
  
  // case 5: mnodes do not project into slave element
  //        snode0 does not project into master element
  //        snode1 does project into master element
  // Note: this case might happen when mnode1 and snode0
  //       project exactly on an opposite node and are assigned
  //       an other element then this one
  if (!snode0 && snode1 && !mnode0 && !mnode1)
  {
    //cout << "Case 5: weirdo projection\n";
    bool ok = true;
    // Have to check whether snodes[0] has a projection 
    // (into a neighbor master segment) and whether that projection point is
    // low in xi range (should be -1.0)
    nstart = snodes[0]->GetProjectedNode(); // check whether a projection exists 
    if (nstart == Teuchos::null) ok = false;
    if (ok) // projection nstart has to be in neighbour master element
    {
      if (!nstart->Segment()) ok = true; // nstart is virtual
      else
      {
        int nseg             = mnodes[1]->Nseg();
        MOERTEL::Segment** segs = mnodes[1]->Segments();
        int segid = nstart->Segment()->Id();
        for (int i=0; i<nseg; ++i)
        if (segid == segs[i]->Id()) { ok = true; break; }
        else ok = false;
      }
    }
    if (ok) sxia = nstart->Xi()[0]; 
    if (ok && sxia > -1.1 && sxia < -0.9) ok = true; // check whether projection is good
    else                                  ok = false;  
    if (ok)
    {    
      nend =  is_spnode1; 
      sxia = -1.0;
      sxib =  1.0;
      mxia =  is_spnode1->Xi()[0];
      mxib =  1.0;
      ++foundcase;
    }
    else // do nothing?
    {
      ++ foundcase;
      nstart = Teuchos::null;
      nend = Teuchos::null; 
    }
  }

  // case 6: both master node project into slave segment
  if (mnode0 && mnode1)
  {
    ++foundcase;
    nstart = is_mpnode0;
    nend   = is_mpnode1;
    sxia = nend->Xi()[0];
    sxib = nstart->Xi()[0];
    mxia = -1.0;
    mxib = 1.0;
  }
  
  // case 7: both slave nodes project into master segment
  if (snode0 && snode1)
  {
    ++foundcase;
    nstart = is_spnode0;
    nend   = is_spnode1;
    sxia = -1.0;
    sxib =  1.0;
    mxia = nend->Xi()[0];
    mxib = nstart->Xi()[0];
  }

  // case 8: first slave node in master segment and first master node in slave segment
  if (snode0 && !snode1 && mnode0 && !mnode1)
  {
    ++foundcase;
    nstart = is_spnode0;
    nend   = is_mpnode0;
    sxia = -1.0;
    sxib = nend->Xi()[0];
    mxia = -1.0;
    mxib = nstart->Xi()[0];
  }

  // case 9: last slave node in master segment and last master node in slave segment
  if (snode1 && !snode0 && mnode1 && !mnode0)
  {
    ++foundcase;
    nstart = is_mpnode1;
    nend   = is_spnode1;
    sxia = nstart->Xi()[0];
    sxib = 1.0;
    mxia = nend->Xi()[0];
    mxib = 1.0;
  }

  if (foundcase != 1)
  {

	std::stringstream oss;
		oss << "***ERR*** MOERTEL::Interface::Integrate_2D_Section:\n"
         << "***ERR*** # cases that apply here: " << foundcase << "\n"
			<< "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n"
			<< "Slave :" << sseg << "\n " << *ssnodes[0] << "   " << *ssnodes[1] << "\n"
			<< "Master:" << mseg << "\n" << *mmnodes[0] << "   " << *mmnodes[1] << "\n"
			<< "snode0: " << snode0 << " snode1: " << snode1 <<
		" mnode0: " << mnode0 << " mnode1: " << mnode1 << "\n";
    throw ReportError(oss);
  }
  
  // there might be no overlap
  if (nstart==Teuchos::null && nend==Teuchos::null)
    return true;

#if 0  
  cout << "slave  xi range " << sxia << " - " << sxib << endl;
  cout << "master xi range " << mxia << " - " << mxib << endl;
#endif
  
  // FIXME: need to get the number of multipliers attached to the slave segment 
  //        when using discontinous lambda, lambdas are attached to segment!
  
  // create an integrator instance of some given order
  MOERTEL::Integrator integrator(5,IsOneDimensional(),OutLevel());
  
  // do the integration of the master side
  Epetra_SerialDenseMatrix* Mdense = 
                            integrator.Integrate(sseg,sxia,sxib,mseg,mxia,mxib);
  
  // do the integration of the slave side
  Epetra_SerialDenseMatrix* Ddense = integrator.Integrate(sseg,sxia,sxib);
  
     // put results -Mdense into Epetra_CrsMatrix M
   // note the sign change for M here
  integrator.Assemble(*this,sseg,mseg,M,*Mdense);

   // put results Ddense into Epetra_CrsMatrix D
  integrator.Assemble(*this,sseg,D,*Ddense);

#if 1 // modification for curved interfaces from paper by B. Wohlmuth
  // do this modification for
  // linear elements
  // vector valued PDE (ndof=2, e.g. elasticity)
  // |delta n| != 0
  if (sseg.Type() == MOERTEL::Segment::seg_Linear1D && 
      mseg.Type() == MOERTEL::Segment::seg_Linear1D)
  if (sseg.FunctionType(1) == MOERTEL::Function::func_DualLinear1D)
  if (snodes[0]->Nlmdof() == snodes[1]->Nlmdof() &&
      mnodes[0]->Ndof() == mnodes[1]->Ndof() &&
      snodes[0]->Nlmdof() == mnodes[0]->Ndof())
  {
    Epetra_SerialDenseMatrix* Mmod = NULL;
    
    // get the normal at slave nodes
    const double* n0 = snodes[0]->N();
    const double* n1 = snodes[1]->N();

    // build the tangential orthogonal to the normal
    double t[2][2];
    t[0][0] = -n0[1]; t[1][0] = -n1[1];
    t[0][1] =  n0[0]; t[1][1] =  n1[0];
    double n[2][2];
    n[0][0] =  n0[0]; n[1][0] =  n1[0]; 
    n[0][1] =  n0[1]; n[1][1] =  n1[1]; 
    
    // build delta values of normal and tangential
    double dn[2]; double dt[2];
    dn[0] = n0[0] - n1[0];  
    dn[1] = n0[1] - n1[1];  
    dt[0] = t[0][0] - t[1][0];
    dt[1] = t[0][1] - t[1][1];
    
    // build norm of dn. If it's zero, don't do anything
    bool doit = true;
//    double delta = dn[0]*dn[0]+dn[1]*dn[1];
//    if (abs(delta)>1.0e-11) doit = true;

    if (doit)
    {
      // do the integration of the modification of the master side
      // integral ( -0.5 * psi_12 * phi_k ) k=1,...,nnode_master 
      Epetra_SerialDenseMatrix* Mmod_scalar =
                        integrator.Integrate_2D_Mmod(sseg,sxia,sxib,mseg,mxia,mxib);

      // create an Epetra_SerialDenseMatrix of dimension (nsnode x nlmdof , nmnode x nmdof)
      int nsnode = sseg.Nnode();
      int nsdof  = snodes[0]->Nlmdof();
      int nmnode = mseg.Nnode();
      int nmdof  = mnodes[0]->Ndof();
      Mmod =  new Epetra_SerialDenseMatrix(nsnode*nsdof,nmnode*nmdof);

      // add modification values to Mmod
      for (int snode=0; snode<nsnode; ++snode)
        for (int sdof=0; sdof<nsdof; ++sdof)
        {
          double nt[2];
          nt[0] = n[snode][sdof] * dn[0] + t[snode][sdof] * dt[0];
          nt[1] = n[snode][sdof] * dn[1] + t[snode][sdof] * dt[1];
          for (int mnode=0; mnode<nmnode; ++mnode)
            for (int mdof=0; mdof<nmdof; ++mdof)
            {
              double val = nt[mdof] * (*Mmod_scalar)(mnode,0);
              (*Mmod)(snode*nsdof+sdof,mnode*nmdof+mdof) = val;
            }
        } // for (int sdof=0; sdof<nsdof; ++sdof)

#if 0  // verification of the expression by expressions given in paper
      Epetra_SerialDenseMatrix* Mmod2 = new Epetra_SerialDenseMatrix(nsnode*nsdof,nmnode*nmdof);
      // n1 dot n2
      double n1n2 = 0.0;
      for (int i=0; i<2; ++i) n1n2 += n[0][i]*n[1][i];
      // third row of n1 x n2
      double n1xn2 = n[0][0]*n[1][1] - n[0][1]*n[1][0];
      
      // slave 0 sdof 0 master 0 mdof 0 
      (*Mmod2)(0,0) = (*Mmod_scalar)(0,0) * (1.0-n1n2);
      // slave 0 sdof 0 master 0 mdof 1
      (*Mmod2)(0,1) = - (*Mmod_scalar)(0,0) * n1xn2;
      // slave 0 sdof 0 master 1 mdof 0
      (*Mmod2)(0,2) = (*Mmod_scalar)(1,0) * (1.0-n1n2);
      // slave 0 sdof 0 master 1 mdof 1
      (*Mmod2)(0,3) = - (*Mmod_scalar)(1,0) * n1xn2;
      // slave 0 sdof 1 master 0 mdof 0 
      (*Mmod2)(1,0) = (*Mmod_scalar)(0,0) * n1xn2;
      // slave 0 sdof 1 master 0 mdof 1
      (*Mmod2)(1,1) = (*Mmod_scalar)(0,0) * (1.0-n1n2);
      // slave 0 sdof 1 master 1 mdof 0
      (*Mmod2)(1,2) = (*Mmod_scalar)(1,0) * n1xn2;
      // slave 0 sdof 1 master 1 mdof 1
      (*Mmod2)(1,3) = (*Mmod_scalar)(1,0) * (1.0-n1n2);
      // slave 1 sdof 0 master 0 mdof 0
      (*Mmod2)(2,0) = (*Mmod_scalar)(0,0) * (n1n2-1.0);
      // slave 1 sdof 0 master 0 mdof 1
      (*Mmod2)(2,1) = - (*Mmod_scalar)(0,0) * n1xn2;
      // slave 1 sdof 0 master 1 mdof 0
      (*Mmod2)(2,2) = (*Mmod_scalar)(1,0) * (n1n2-1.0);
      // slave 1 sdof 0 master 1 mdof 1
      (*Mmod2)(2,3) = - (*Mmod_scalar)(1,0) * n1xn2;
      // slave 1 sdof 1 master 0 mdof 0
      (*Mmod2)(3,0) = (*Mmod_scalar)(0,0) * n1xn2;
      // slave 1 sdof 1 master 0 mdof 1
      (*Mmod2)(3,1) = (*Mmod_scalar)(0,0) * (n1n2-1.0);
      // slave 1 sdof 1 master 1 mdof 0
      (*Mmod2)(3,2) = (*Mmod_scalar)(1,0) * n1xn2;
      // slave 1 sdof 1 master 1 mdof 1
      (*Mmod2)(3,3) = (*Mmod_scalar)(1,0) * (n1n2-1.0);
      //cout << *Mmod2;
      //delete Mmod2; Mmod2 = NULL;
#endif

      //  assemble -Mmod into M
      integrator.Assemble_2D_Mod(*this,sseg,mseg,M,*Mmod);
      
      // tidy up 
      if (Mmod)        delete Mmod;        Mmod = NULL;
      if (Mmod_scalar) delete Mmod_scalar; Mmod_scalar = NULL;
    } // if (doit)
  } // if a lot of stuff
#endif

  
  // tidy up 
  if (Mdense) delete Mdense; Mdense = NULL;
  if (Ddense) delete Ddense; Ddense = NULL;

  return true;
}
#endif
