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
#include "mrtr_projector.H"
#include "mrtr_overlap.H"

#include "Epetra_SerialDenseMatrix.h"
#include "Epetra_Time.h"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"


const double CONSTRAINT_MATRIX_ZERO = 1.0e-11;

/*----------------------------------------------------------------------*
 |  make mortar integration of this interface (2D/3D problem)           |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Mortar_Integrate(
    Teuchos::RCP<Teuchos::ParameterList> intparams) {
  bool ok = false;
  intparams_ = intparams;
  
  //-------------------------------------------------------------------
  // time this process
  Epetra_Time time(*lComm());
  time.ResetStartTime();

  //-------------------------------------------------------------------
	if(IsOneDimensional()) {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** This is not a 3D problem, we're in the wrong method here!!!\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";

    return false;
  }

  //-------------------------------------------------------------------
  // interface needs to be complete
	if(!IsComplete()) {
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
	if(MortarSide() == -1) {
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
		if(scurr->second->Nfunctions() < 1) {
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", mortar side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 1 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }

  for (scurr=seg_[sside].begin(); scurr!=seg_[sside].end(); ++scurr)
		if(scurr->second->Nfunctions() < 2) {
      cout << "***ERR*** MOERTEL::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", slave side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 2 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    
  //-------------------------------------------------------------------
  // do the integration of the master and slave side
  ok = Integrate_3D();

  if (!ok) return false;

  //-------------------------------------------------------------------
  // set the flag that this interface has been successfully integrated
  isIntegrated_ = true;
  
  //-------------------------------------------------------------------
  // time this process
	if(OutLevel() > 5) {
    cout << "MOERTEL::Interface " << Id() << ": Integration on proc " << gComm().MyPID()
		     << " finished in " << time.ElapsedTime() << " sec\n";
		fflush(stdout);
  }

  //-------------------------------------------------------------------
  return true;
}

/*----------------------------------------------------------------------*
 |  make mortar integration of master/slave side in 3D (2D interface)   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_3D() {

	if(!IsComplete()) {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MOERTEL::Interface::Integrate_3D:\n"
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

	for(scurr = rseg_[sside].begin(); scurr != rseg_[sside].end(); ++scurr) {
    // the segment to be integrated
	Teuchos::RCP<MOERTEL::Segment> actsseg = scurr->second;

#if 0
    cout << "\n\nActive sseg id " << actsseg->Id() << "\n";
#endif

    // check whether I own at least one of the nodes on this slave segment
    const int nnode = actsseg->Nnode();
    MOERTEL::Node** nodes = actsseg->Nodes();
    bool foundone = false;

    for (int i=0; i<nnode; ++i)
			if(NodePID(nodes[i]->Id()) == lComm()->MyPID()) {
        foundone = true;
        break;
      }

    // if none of the nodes belongs to me, do nothing on this segment
    if (!foundone) continue;
    
    // time this process
    //Epetra_Time time(*lComm());
    //time.ResetStartTime();

    // loop over all segments on the master side
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator mcurr;


		for(mcurr = rseg_[mside].begin(); mcurr != rseg_[mside].end(); ++mcurr) {
	  Teuchos::RCP<MOERTEL::Segment> actmseg = mcurr->second;
#if 0
      cout << "Active mseg id " << actmseg->Id() << endl;
#endif

      
      // if there is an overlap, integrate the pair
      // (whether there is an overlap or not will be checked inside)
      Integrate_3D_Section(*actsseg,*actmseg);
      
    } // for (mcurr=rseg_[mside].begin(); mcurr!=rseg_[mside].end(); ++mcurr)  


    //cout << "time for this slave segment: " << time.ElapsedTime() << endl;

  } // for (scurr=rseg_[sside].begin(); scurr!=rseg_[sside].end(); ++scurr)

  return true;
}

/*----------------------------------------------------------------------*
 | integrate the master/slave side's contribution from the overlap      |
 | of 2 segments (3D version) IF there is an overlap                    |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Integrate_3D_Section(MOERTEL::Segment& sseg, 
                                           MOERTEL::Segment& mseg){ 

  if ( (sseg.Type()!=MOERTEL::Segment::seg_BiLinearTri &&
        sseg.Type()!=MOERTEL::Segment::seg_BiLinearQuad    )  || 
       (mseg.Type()!=MOERTEL::Segment::seg_BiLinearTri && 
	         mseg.Type() != MOERTEL::Segment::seg_BiLinearQuad)) {
		std::stringstream oss;
		oss << "***ERR*** MOERTEL::Interface::Integrate_3D_Section:\n"
         << "***ERR*** Integration of other then bilinear triangles/quads not implemented\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
  }
  
  // find whether we want exact values at gaussian points
  bool exactvalues = intparams_->get("exact values at gauss points",true);
  
  // first determine whether there is an overlap between sseg and mseg
  // for this purpose, the 'overlapper' class is used
  // It also builds a triangulation of the overlap polygon if there is any
  MOERTEL::Overlap overlap(sseg,mseg,*this,exactvalues,OutLevel());

  // determine the overlap triangulation if any
  bool ok = overlap.ComputeOverlap();

  if (!ok) return true; // There's no overlap

  // # new segments the overlap polygon was discretized with
  int nseg = overlap.Nseg();
  // view of these segments
  std::vector<Teuchos::RCP<MOERTEL::Segment> > segs;
  overlap.SegmentView(segs);
  
  // integrator object
  int ngp = intparams_->get("number gaussian points 2D",12);
  MOERTEL::Integrator integrator(ngp,IsOneDimensional(),OutLevel());
  
  // loop segments and integrate them
	for(int s = 0; s < nseg; ++s) {
	Teuchos::RCP<MOERTEL::Segment> actseg = segs[s];

    // integrate master and slave part of this segment
    Epetra_SerialDenseMatrix* Ddense = NULL;
    Epetra_SerialDenseMatrix* Mdense = NULL;
    bool ok = integrator.Integrate(actseg,sseg,mseg,&Ddense,&Mdense,overlap,1.0e-04,
                                   exactvalues);
    if (!ok)
      continue;
    
    // assemble temporarily into the nodes
    integrator.Assemble(*this,sseg,*Ddense);    
    integrator.Assemble(*this,sseg,mseg,*Mdense);    
    
    if (Ddense) delete Ddense;

    if (Mdense) delete Mdense;
          
  } // for (int s=0; s<nseg; ++s)

  segs.clear();
  
  return true;
}


/*----------------------------------------------------------------------*
 |  assemble integration of master/slave side in 3D (2D interface)      |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Assemble_3D(Epetra_CrsMatrix& D, Epetra_CrsMatrix& M) {

	if(!IsComplete()) {
    cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
         << "***ERR*** Complete() not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  //-------------------------------------------------------------------
  // loop over all slave nodes
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;

	for(curr = rnode_[sside].begin(); curr != rnode_[sside].end(); ++curr) {
    // loop only my own nodes
    if (NodePID(curr->second->Id()) != lComm()->MyPID())
      continue;
    
    // get std::maps D and M and Mmod from node
	Teuchos::RCP<std::map<int,double> >          Drow = curr->second->GetD();
	Teuchos::RCP<std::map<int,double> >          Mrow = curr->second->GetM();
	Teuchos::RCP<std::vector<std::map<int,double> > > Mmod = curr->second->GetMmod();

    // if there's no D or M there's nothing to do
    if (Drow==Teuchos::null && Mrow==Teuchos::null)
      continue;
    
	Teuchos::RCP<MOERTEL::Node> rowsnode = curr->second;
    int snlmdof = rowsnode->Nlmdof();
    const int* slmdof = rowsnode->LMDof();
    //cout << "Current row snode: " << rowsnode->Id() << endl;

   std::map<int,double>::iterator rowcurr;

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
    // assemble the Drow
		if(Drow != Teuchos::null) {
			for(rowcurr = Drow->begin(); rowcurr != Drow->end(); ++rowcurr) {
        int colnode = rowcurr->first;
        double val  = rowcurr->second;

				if(abs(val) < CONSTRAINT_MATRIX_ZERO)
          continue;
          
        //cout << "Current col snode: " << colnode << endl;
        
        // get the colsnode
		Teuchos::RCP<MOERTEL::Node> colsnode = GetNodeView(colnode);

				if(colsnode == Teuchos::null) {
          cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
               << "***ERR*** interface " << Id_ << ": cannot get view of node " << colnode << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
        // get the primal dofs
        int sndof = colsnode->Ndof();
        const int* sdof = colsnode->Dof();

				if(snlmdof != sndof) {
          cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
               << "***ERR*** interface " << Id_ << ": mismatch in # lagrange multipliers and primal variables\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
				for(int i = 0; i < snlmdof; ++i) {
          int row = slmdof[i];
          int col = sdof[i];
           //cout << "Inserting D row/col:" << row << "/" << col << " val " << val << endl;
          int err = D.SumIntoGlobalValues(row,1,&val,&col);

          if (err)
            err = D.InsertGlobalValues(row,1,&val,&col);

					if(err < 0) {
            cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                 << "***ERR*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
            return false;
          }

					if(err && OutLevel() > 0) {
            cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
                 << "MOERTEL: ***WRN*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                 << "MOERTEL: ***WRN*** indicating that initial guess for memory of D too small\n"
                 << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          } 
        } // for (int i=0; i<snlmdof; ++i)
      } // for (rowcurr=Drow->begin(); rowcurr!=Drow->end(); ++rowcurr)
    }

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
    // assemble the Mrow
		if(Mrow != Teuchos::null) {
			for(rowcurr = Mrow->begin(); rowcurr != Mrow->end(); ++rowcurr) {
        int colnode = rowcurr->first;
        double val  = rowcurr->second;

				if(abs(val) < CONSTRAINT_MATRIX_ZERO)
          continue;
          
        // cout << "Current colmnode: " << colnode << endl;
        
        // get the colmnode
		Teuchos::RCP<MOERTEL::Node> colmnode = GetNodeView(colnode);

				if(colmnode == Teuchos::null) {
          cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
               << "***ERR*** interface " << Id_ << ": cannot get view of node " << colnode << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
        // get the primal dofs
        int mndof = colmnode->Ndof();
        const int* mdof = colmnode->Dof();

				if(snlmdof != mndof) {
          cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
               << "***ERR*** interface " << Id_ << ": mismatch in # lagrange multipliers and primal variables\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
				for(int i = 0; i < snlmdof; ++i) {
          int row = slmdof[i];
          int col = mdof[i];
          // cout << "Inserting M row/col:" << row << "/" << col << " val " << val << endl;
          int err = M.SumIntoGlobalValues(row,1,&val,&col);

          if (err)
            err = M.InsertGlobalValues(row,1,&val,&col);

					if(err < 0) {
            cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                 << "***ERR*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
            return false;
          }

					if(err && OutLevel() > 0) {
            cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
                 << "MOERTEL: ***WRN*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                 << "MOERTEL: ***WRN*** indicating that initial guess for memory of M too small\n"
                 << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          } 
        } // for (int i=0; i<snlmdof; ++i)
      } // for (rowcurr=Mrow->begin(); rowcurr!=Mrow->end(); ++rowcurr)
    }
    
    
    // assemble the Mmod block if there is any
		if(Mmod != Teuchos::null) {
      // loop over he rows of the Mmod block
			for(int lrow = 0; lrow < (int)Mmod->size(); ++lrow) {
//        std::map<int,double>& Mmodrow = (*Mmod)[lrow];
        int row = slmdof[lrow];
        
        // loop over the columns in that row
        // FIXMEL: should this be Mmodrow (above)?
				for(rowcurr = Mrow->begin(); rowcurr != Mrow->end(); ++rowcurr) {
          int col = rowcurr->first;
          double val = rowcurr->second;

					if(abs(val) < CONSTRAINT_MATRIX_ZERO)
            continue;
          
          //cout << "Inserting M row/col:" << row << "/" << col << " val " << val << endl;
          int err = M.SumIntoGlobalValues(row,1,&val,&col);

          if (err)
            err = M.InsertGlobalValues(row,1,&val,&col);

					if(err < 0) {
            cout << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                 << "***ERR*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
            return false;
          }

					if(err && OutLevel() > 0) {
            cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
                 << "MOERTEL: ***WRN*** interface " << Id_ << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                 << "MOERTEL: ***WRN*** indicating that initial guess for memory of M too small\n"
                 << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          } 
          
        
        } // for (rowcurr=Mrow->begin(); rowcurr!=Mrow->end(); ++rowcurr)
        
      } // for (int lrow=0; lrow<(int)Mmod->size(); ++lrow)
    }

  } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  // In case this interface was parallel we might have missed something upto 
  // here. Boundary terms of D and M, Mmod are assembled non-local (that is to
  // close inner-interface nodes). If these inner-interface nodes belong
  // to a different proc values were not assembled.
  // Loop snodes again an check and communicate these entries
	if(lComm()->NumProc() != 1) {
      // note that we miss the communication of Mmod yet
      
      // allocate a sendbuffer for D and M
      int countD=0;
      int countM=0;
	  std::vector<int> colD_s(countD);
	  std::vector<double> valD_s(countD);
	  std::vector<int> colM_s(countM);
	  std::vector<double> valM_s(countM);

		for(curr = rnode_[sside].begin(); curr != rnode_[sside].end(); ++curr) {
        // we've done all my own nodes already
        if (NodePID(curr->second->Id()) == lComm()->MyPID())
          continue;
        
        // check whether we have M or D values here
        // get maps D and M from node
		Teuchos::RCP<std::map<int,double> > Drow = curr->second->GetD();
		Teuchos::RCP<std::map<int,double> > Mrow = curr->second->GetM();
    
        // if there's no D/M there's nothing to do
        if (Drow==Teuchos::null && Mrow==Teuchos::null)
          continue;
        
        //cout << "lProc " << lComm()->MyPID() << " Node " << curr->second->Id() << " unassembled\n"; 
        
        // fill the D sendbuffer
			if(Drow != Teuchos::null) {
          // resize the sendbuffers
          colD_s.resize(colD_s.size()+Drow->size()+2);
          valD_s.resize(valD_s.size()+Drow->size()+2);
          // Add node Id and size
          colD_s[countD]   = curr->second->Id();
          valD_s[countD]   = 0.0;
          colD_s[countD+1] = (int)Drow->size();
          valD_s[countD+1] = 0.0;
          countD += 2;
          // loop D
		  std::map<int,double>::iterator rowcurr;

				for(rowcurr = Drow->begin(); rowcurr != Drow->end(); ++rowcurr) {
            colD_s[countD] = rowcurr->first;
            valD_s[countD] = rowcurr->second;
            ++countD;
          }
        }

        //cout << "lProc " << lComm()->MyPID() << " Node " << curr->second->Id() << " countD " << countD << endl;
        
        // fill the M sendbuffer
			if(Mrow != Teuchos::null) {
          // resize the sendbuffers
          colM_s.resize(colM_s.size()+Mrow->size()+2);
          valM_s.resize(valM_s.size()+Mrow->size()+2);
          // Add node id and size
          colM_s[countM] = curr->second->Id();
          valM_s[countM] = 0.0;
          colM_s[countM+1] = Mrow->size();
          valM_s[countM+1] = 0.0;
          countM += 2;
          // loop M
		  std::map<int,double>::iterator rowcurr;

				for(rowcurr = Mrow->begin(); rowcurr != Mrow->end(); ++rowcurr) {
            colM_s[countM] = rowcurr->first;
            valM_s[countM] = rowcurr->second;
            ++countM;
          }
        }
        
        //cout << "lProc " << lComm()->MyPID() << " Node " << curr->second->Id() << " countM " << countM << endl;
      } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
      
      
      // loop all processes in lComm and communicate and assemble
		for(int proc = 0; proc < lComm()->NumProc(); ++proc) {
        // send sizes
        int countDr = countD;
        int countMr = countM;
        lComm()->Broadcast(&countDr,1,proc);
        lComm()->Broadcast(&countMr,1,proc);
        // allocate receive buffers
		std::vector<int>    colD_r(countDr);
		std::vector<double> valD_r(countDr);
		std::vector<int>    colM_r(countMr);
		std::vector<double> valM_r(countMr);

        // send data
			if(proc == lComm()->MyPID()) {
				for(int i = 0; i < countDr; ++i) {
            colD_r[i] = colD_s[i];
            valD_r[i] = valD_s[i];
          }

				for(int i = 0; i < countMr; ++i) {
            colM_r[i] = colM_s[i];
            valM_r[i] = valM_s[i];
          }
        }
		if(countDr > 0){

			lComm()->Broadcast(&colD_r[0],countDr,proc);
			lComm()->Broadcast(&valD_r[0],countDr,proc);
		}

		if(countMr > 0){

			lComm()->Broadcast(&colM_r[0],countMr,proc);
			lComm()->Broadcast(&valM_r[0],countMr,proc);
		}
        
        // Assemble (remote procs only)
			if(proc != lComm()->MyPID()) {
          // --------------------------------------------------- Assemble D
				for(int i = 0; i < countDr;) {
            int nodeid = colD_r[i];
            int size   = colD_r[i+1];
            i += 2;

            // find whether I am owner of this node
					if(NodePID(nodeid) == lComm()->MyPID()) {
              // get the node
			  Teuchos::RCP<MOERTEL::Node> snode = GetNodeView(nodeid);

						if(snode == Teuchos::null) {
							std::stringstream oss;
							oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                     << "***ERR*** Cannot find view of node " << nodeid << endl
                     << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
							throw ReportError(oss);
              }
              
              // get lagrange multipliers
              int nslmdof = snode->Nlmdof();
              const int* slmdof = snode->LMDof();
              
              // loop colD_r/valD_r and assemble
						for(int j = 0; j < size; ++j) {
                int colsnode = colD_r[i+j];
                double val   = valD_r[i+j];

							if(abs(val) < CONSTRAINT_MATRIX_ZERO)
                  continue;
                
                // get view of column node and primal dofs
				Teuchos::RCP<MOERTEL::Node> colnode = GetNodeView(colsnode);

							if(colnode == Teuchos::null) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                       << "***ERR*** Cannot find view of node " << colsnode << endl
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }

                int       nsdof = colnode->Ndof();
                const int* sdof = colnode->Dof();

							if(nsdof != nslmdof) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                       << "***ERR*** Mismatch in # primal dofs and lagrange mutlipliers\n"
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }
                
							for(int k = 0; k < nslmdof; ++k) {
                  int row = slmdof[k];
                  int col = sdof[k];
                  //cout << "Proc " << lComm()->MyPID() << " inserting D row/col:" << row << "/" << col << " val " << val << endl;
                  int err = D.SumIntoGlobalValues(row,1,&val,&col);

                  if (err)
                    err = D.InsertGlobalValues(row,1,&val,&col);

								if(err < 0) {
									std::stringstream oss;
									oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                         << "***ERR*** Serious error=" << err << " in assembly\n"
                         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
									throw ReportError(oss);
                  }

								if(err && OutLevel() > 0) {
                    cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
                         << "MOERTEL: ***WRN*** interface " << Id() << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                         << "MOERTEL: ***WRN*** indicating that initial guess for memory of D too small\n"
                         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
                  }
                } // for (int k=0; k<nslmdof; ++k)
              } // for (int j=0; j<size; ++j)

              i += size;
            }

            else // I am not owner of this node, skip it
              i += size;
          } // for (int i=0; i<countDr;)
          
          // --------------------------------------------------- Assemble M
				for(int i = 0; i < countMr;) {
            int nodeid = colM_r[i];
            int size   = colM_r[i+1];
            i += 2;

            // find whether I am owner of this node
					if(NodePID(nodeid) == lComm()->MyPID()) {
              // get the node
			  Teuchos::RCP<MOERTEL::Node> snode = GetNodeView(nodeid);

						if(snode == Teuchos::null) {
							std::stringstream oss;
							oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                     << "***ERR*** Cannot find view of node " << nodeid << endl
                     << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
							throw ReportError(oss);
              }
              
              // get the lagrange multipliers
              int nslmdof = snode->Nlmdof();
              const int* slmdof = snode->LMDof();
              
              // loop colM_r/valM_r and assemble
						for(int j = 0; j < size; ++j) {
                int colmnode = colM_r[i+j];
                double val   = valM_r[i+j];

							if(abs(val) < CONSTRAINT_MATRIX_ZERO)
                  continue;
                
                // get view of column node and primal dofs
				Teuchos::RCP<MOERTEL::Node> colnode = GetNodeView(colmnode);

							if(colnode == Teuchos::null) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                       << "***ERR*** Cannot find view of node " << colmnode << endl
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }

                int nmdof = colnode->Ndof();
                const int* mdof = colnode->Dof();

							if(nmdof != nslmdof) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                       << "***ERR*** Mismatch in # primal dofs and lagrange mutlipliers\n"
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }
                
							for(int k = 0; k < nslmdof; ++k) {
                  int row = slmdof[k];
                  int col = mdof[k];
                  //cout << "Proc " << lComm()->MyPID() << " inserting M row/col:" << row << "/" << col << " val " << val << endl;
                  int err = M.SumIntoGlobalValues(row,1,&val,&col);

                  if (err)
                    err = M.InsertGlobalValues(row,1,&val,&col);

								if(err < 0) {
									std::stringstream oss;
									oss << "***ERR*** MOERTEL::Interface::Assemble_3D:\n"
                         << "***ERR*** Serious error=" << err << " in assembly\n"
                         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
									throw ReportError(oss);
                  }

								if(err && OutLevel() > 0) {
                    cout << "MOERTEL: ***WRN*** MOERTEL::Interface::Assemble_3D:\n"
                         << "MOERTEL: ***WRN*** interface " << Id() << ": Epetra_CrsMatrix::InsertGlobalValues returned " << err << "\n"
                         << "MOERTEL: ***WRN*** indicating that initial guess for memory of M too small\n"
                         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
                  }
                } // for (int k=0; k<nslmdof; ++k)
              } // for (int j=0; j<size; ++j)

              i += size;
            }

            else // I am not owner of this node, skip it
              i += size;
          } // for (int i=0; i<countMr;)
        } // if (proc!=lComm()->MyPID())
        
        colD_r.clear();
        valD_r.clear();
        colM_r.clear();
        valM_r.clear();
      } // for (int proc=0; proc<lComm()->NumProc(); ++proc)
  
    colD_s.clear();
    valD_s.clear();
    colM_s.clear();
    valM_s.clear();
  } // if (lComm()->NumProc()!=1)
  

  return true;
}

/*----------------------------------------------------------------------*
 |  assemble integration of master/slave side into the residual vector (JFNK)
 | 
 |
 | Here, each global node (rnode_) is visited on the slave interface side.
 | Each node stores a local D, M, and Mmod.
 *----------------------------------------------------------------------*/
//#define PDANDM
bool MOERTEL::Interface::AssembleJFNKVec(Lmselector *sel) {

	if(!IsComplete()) {
		cout << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
         << "***ERR*** Complete() not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  if (!lComm()) return true;

  // get the sides
  int mside = MortarSide();
  int sside = OtherSide(mside);

  //-------------------------------------------------------------------
  // loop over all slave nodes
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;

#ifdef PDANDM  // Save and print the D and M for debugging
	int size = rnode_[sside].size();
	int cnt = 0;
	std::vector<int> dtable(size);

  for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)

		dtable[cnt++] = curr->second->Id();

	Epetra_Map Dmap(-1, size, &dtable[0], 0, *lComm());
	Epetra_CrsMatrix Dmat(Copy, Dmap, 4);
	Epetra_CrsMatrix Mmat(Copy, Dmap, 4);

#endif

	for(curr = rnode_[sside].begin(); curr != rnode_[sside].end(); ++curr) {

    // loop only my own nodes
    if (NodePID(curr->second->Id()) != lComm()->MyPID())
      continue;

//	cout << curr->second->Id() << endl;
    
    // get maps D and M and Mmod from node
	Teuchos::RCP<std::map<int,double> >          Drow = curr->second->GetD();
	Teuchos::RCP<std::map<int,double> >          Mrow = curr->second->GetM();
	Teuchos::RCP<std::vector<std::map<int,double> > > Mmod = curr->second->GetMmod();

    // if there's no D or M there's nothing to do
    if (Drow==Teuchos::null && Mrow==Teuchos::null)
      continue;
    
	Teuchos::RCP<MOERTEL::Node> rowsnode = curr->second;
    int snlmdof = rowsnode->Nlmdof();
    const int* slmdof = rowsnode->LMDof();
    //cout << "Current row snode: " << rowsnode->Id() << endl;

   std::map<int,double>::iterator rowcurr;

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
    // assemble the Drow
		if(Drow != Teuchos::null) {
			for(rowcurr = Drow->begin(); rowcurr != Drow->end(); ++rowcurr) {
        int colnode = rowcurr->first;
        double val  = rowcurr->second;

				if(abs(val) < CONSTRAINT_MATRIX_ZERO)   // this entry of D is effectively zero
          continue;
#ifdef PDANDM  // Save the row, col, and value
	Dmat.InsertGlobalValues(curr->second->Id(), 1, &val, &colnode);
#endif
          
        //cout << "Current col snode: " << colnode << endl;
        
        // get the colsnode
		Teuchos::RCP<MOERTEL::Node> colsnode = GetNodeView(colnode);

				if(colsnode == Teuchos::null) {
					cout << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
               << "***ERR*** interface " << Id_ << ": cannot get view of node " << colnode << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
        // get the primal dofs
        int sndof = colsnode->Ndof();
        const int* sdof = colsnode->Dof();

				if(snlmdof != sndof) {
					cout << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
               << "***ERR*** interface " << Id_ << ": mismatch in # lagrange multipliers and primal variables\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
				for(int i = 0; i < snlmdof; ++i) {

					if(!sel->EvaluateLM(rowsnode, i)) // Continue if this LM is not active
				  continue;

          int row = slmdof[i];
          int col = sdof[i];
//           cout << "Inserting D row/col:" << row << "/" << col << " val " << val << endl;

		   /*
          int err = D.SumIntoGlobalValues(row,1,&val,&col);
          if (err)
            err = D.InsertGlobalValues(row,1,&val,&col);
			*/

		  // Assemble D times soln
		  // Row of D determines row in rhs
		  // col of D determines row of soln
		  
		  sel->AssembleNodeVal(row, col, val);

        } // for (int i=0; i<snlmdof; ++i)
      } // for (rowcurr=Drow->begin(); rowcurr!=Drow->end(); ++rowcurr)
    }

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
    // assemble the Mrow
		if(Mrow != Teuchos::null) {
			for(rowcurr = Mrow->begin(); rowcurr != Mrow->end(); ++rowcurr) {
        int colnode = rowcurr->first;
        double val  = rowcurr->second;

				if(abs(val) < CONSTRAINT_MATRIX_ZERO)
          continue;
          
        // cout << "Current colmnode: " << colnode << endl;
				//
#ifdef PDANDM  // Save the row, col, and value
	Mmat.InsertGlobalValues(curr->second->Id(), 1, &val, &colnode);
#endif
        
        // get the colmnode
		Teuchos::RCP<MOERTEL::Node> colmnode = GetNodeView(colnode);

				if(colmnode == Teuchos::null) {
					cout << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
               << "***ERR*** interface " << Id_ << ": cannot get view of node " << colnode << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
        // get the primal dofs
        int mndof = colmnode->Ndof();
        const int* mdof = colmnode->Dof();

				if(snlmdof != mndof) {
					cout << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
               << "***ERR*** interface " << Id_ << ": mismatch in # lagrange multipliers and primal variables\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          return false;
        }
        
				for(int i = 0; i < snlmdof; ++i) {

          if(!sel->EvaluateLM(rowsnode, i)) // true if this LM is active
				  continue;

          int row = slmdof[i];
          int col = mdof[i];
          // cout << "Inserting M row/col:" << row << "/" << col << " val " << val << endl;
		  
		  // Assemble M times soln
		  // Row of M determines row in rhs
		  // col of M determines row of soln

		  sel->AssembleNodeVal(row, col, val);

		  /*
          int err = M.SumIntoGlobalValues(row,1,&val,&col);
          if (err)
            err = M.InsertGlobalValues(row,1,&val,&col);
		  */

        } // for (int i=0; i<snlmdof; ++i)
      } // for (rowcurr=Mrow->begin(); rowcurr!=Mrow->end(); ++rowcurr)
    }
    
    
    // assemble the Mmod block if there is any
		if(Mmod != Teuchos::null) {
      // loop over the rows of the Mmod block

			std::stringstream oss;
			oss << "Mmod has entries in it" << endl;
			throw ReportError(oss);

			for(int lrow = 0; lrow < (int)Mmod->size(); ++lrow) {

//        std::map<int,double>& Mmodrow = (*Mmod)[lrow];
//
        int row = slmdof[lrow];
        
        // loop over the columns in that row
        // FIXMEL: should this be Mmodrow (above)?

				for(rowcurr = Mrow->begin(); rowcurr != Mrow->end(); ++rowcurr) {

          if(!sel->EvaluateLM(rowsnode, lrow)) // true if this LM is active
				  continue;

          int col = rowcurr->first;
          double val = rowcurr->second;

					if(abs(val) < CONSTRAINT_MATRIX_ZERO)
            continue;

		  // Assemble D times soln
		  // Row of D determines row in rhs
		  // col of D determines row of soln

		  sel->AssembleNodeVal(row, col, val);

		  /*
          int err = M.SumIntoGlobalValues(row,1,&val,&col);
          if (err)
            err = M.InsertGlobalValues(row,1,&val,&col);
			*/
        
        } // for (rowcurr=Mrow->begin(); rowcurr!=Mrow->end(); ++rowcurr)
        
      } // for (int lrow=0; lrow<(int)Mmod->size(); ++lrow)
    }

		sel->AccumulateRHS(rowsnode);

  } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)

	// Master side - Is this even needed??? GAH

	/*
  for (curr=rnode_[mside].begin(); curr!=rnode_[mside].end(); ++curr){

	Teuchos::RCP<MOERTEL::Node> rowsnode = curr->second;

    sel->AccumulateMRHS(rowsnode);

  }
	*/

  //-------------------------------------------------------------------
  //-------------------------------------------------------------------
  // In case this interface was parallel we might have missed something upto 
  // here. Boundary terms of D and M, Mmod are assembled non-local (that is to
  // close inner-interface nodes). If these inner-interface nodes belong
  // to a different proc values were not assembled.
  // Loop snodes again an check and communicate these entries
	
	if(lComm()->NumProc() != 1) {

      // note that we miss the communication of Mmod yet
      
      // allocate a sendbuffer for D and M
      int countD=0;
      int countM=0;
	  std::vector<int> colD_s(countD);
	  std::vector<double> valD_s(countD);
	  std::vector<int> colM_s(countM);
	  std::vector<double> valM_s(countM);

		for(curr = rnode_[sside].begin(); curr != rnode_[sside].end(); ++curr) {
        // we've done all my own nodes already
        if (NodePID(curr->second->Id()) == lComm()->MyPID())
          continue;
        
        // check whether we have M or D values here
        // get maps D and M from node
		Teuchos::RCP<std::map<int,double> > Drow = curr->second->GetD();
		Teuchos::RCP<std::map<int,double> > Mrow = curr->second->GetM();
    
        // if there's no D/M there's nothing to do
        if (Drow==Teuchos::null && Mrow==Teuchos::null)
          continue;
        
        //cout << "lProc " << lComm()->MyPID() << " Node " << curr->second->Id() << " unassembled\n"; 
        
        // fill the D sendbuffer
			if(Drow != Teuchos::null) {
          // resize the sendbuffers
          colD_s.resize(colD_s.size()+Drow->size()+2);
          valD_s.resize(valD_s.size()+Drow->size()+2);
          // Add node Id and size
          colD_s[countD]   = curr->second->Id();
          valD_s[countD]   = 0.0;
          colD_s[countD+1] = (int)Drow->size();
          valD_s[countD+1] = 0.0;
          countD += 2;
          // loop D
		  std::map<int,double>::iterator rowcurr;

				for(rowcurr = Drow->begin(); rowcurr != Drow->end(); ++rowcurr) {
            colD_s[countD] = rowcurr->first;
            valD_s[countD] = rowcurr->second;
            ++countD;
          }
        }

        //cout << "lProc " << lComm()->MyPID() << " Node " << curr->second->Id() << " countD " << countD << endl;
        
        // fill the M sendbuffer
			if(Mrow != Teuchos::null) {
          // resize the sendbuffers
          colM_s.resize(colM_s.size()+Mrow->size()+2);
          valM_s.resize(valM_s.size()+Mrow->size()+2);
          // Add node id and size
          colM_s[countM] = curr->second->Id();
          valM_s[countM] = 0.0;
          colM_s[countM+1] = Mrow->size();
          valM_s[countM+1] = 0.0;
          countM += 2;
          // loop M
		  std::map<int,double>::iterator rowcurr;

				for(rowcurr = Mrow->begin(); rowcurr != Mrow->end(); ++rowcurr) {
            colM_s[countM] = rowcurr->first;
            valM_s[countM] = rowcurr->second;
            ++countM;
          }
        }
        
        //cout << "lProc " << lComm()->MyPID() << " Node " << curr->second->Id() << " countM " << countM << endl;
      } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
      
      
      // loop all processes in lComm and communicate and assemble
		for(int proc = 0; proc < lComm()->NumProc(); ++proc) {
        // send sizes
        int countDr = countD;
        int countMr = countM;
        lComm()->Broadcast(&countDr,1,proc);
        lComm()->Broadcast(&countMr,1,proc);
        // allocate receive buffers
		std::vector<int>    colD_r(countDr);
		std::vector<double> valD_r(countDr);
		std::vector<int>    colM_r(countMr);
		std::vector<double> valM_r(countMr);

        // send data
			if(proc == lComm()->MyPID()) {
				for(int i = 0; i < countDr; ++i) {
            colD_r[i] = colD_s[i];
            valD_r[i] = valD_s[i];
          }

				for(int i = 0; i < countMr; ++i) {
            colM_r[i] = colM_s[i];
            valM_r[i] = valM_s[i];
          }
        }

        lComm()->Broadcast(&colD_r[0],countDr,proc);
        lComm()->Broadcast(&valD_r[0],countDr,proc);
        lComm()->Broadcast(&colM_r[0],countMr,proc);
        lComm()->Broadcast(&valM_r[0],countMr,proc);
        
        // Assemble (remote procs only)
			if(proc != lComm()->MyPID()) {
          // --------------------------------------------------- Assemble D
				for(int i = 0; i < countDr;) {
            int nodeid = colD_r[i];
            int size   = colD_r[i+1];
            i += 2;

            // find whether I am owner of this node
					if(NodePID(nodeid) == lComm()->MyPID()) {
              // get the node
			  Teuchos::RCP<MOERTEL::Node> snode = GetNodeView(nodeid);

						if(snode == Teuchos::null) {
							std::stringstream oss;
							oss << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
                     << "***ERR*** Cannot find view of node " << nodeid << endl
                     << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
							throw ReportError(oss);
              }
              
              // get lagrange multipliers
              int nslmdof = snode->Nlmdof();
              const int* slmdof = snode->LMDof();
              
              // loop colD_r/valD_r and assemble
						for(int j = 0; j < size; ++j) {
                int colsnode = colD_r[i+j];
                double val   = valD_r[i+j];

							if(abs(val) < CONSTRAINT_MATRIX_ZERO)
                  continue;
                
                // get view of column node and primal dofs
				Teuchos::RCP<MOERTEL::Node> colnode = GetNodeView(colsnode);

							if(colnode == Teuchos::null) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
                       << "***ERR*** Cannot find view of node " << colsnode << endl
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }

                int       nsdof = colnode->Ndof();
                const int* sdof = colnode->Dof();

							if(nsdof != nslmdof) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
                       << "***ERR*** Mismatch in # primal dofs and lagrange mutlipliers\n"
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }
                
							for(int k = 0; k < nslmdof; ++k) {

                  if(!sel->EvaluateLM(snode, k)) // true if this LM is active
						continue;

                  int row = slmdof[k];
                  int col = sdof[k];
                  //cout << "Proc " << lComm()->MyPID() << " inserting D row/col:" << row << "/" << col << " val " << val << endl;
				  
				  // Assemble D times soln
				  // Row of D determines row in rhs
				  // col of D determines row of soln


		  		  sel->AssembleNodeVal(row, col, val);

				  /*
                  int err = D.SumIntoGlobalValues(row,1,&val,&col);
                  if (err)
                    err = D.InsertGlobalValues(row,1,&val,&col);
				  */

                } // for (int k=0; k<nslmdof; ++k)
              } // for (int j=0; j<size; ++j)

              i += size;
            }

            else // I am not owner of this node, skip it
              i += size;
          } // for (int i=0; i<countDr;)
          
          // --------------------------------------------------- Assemble M
				for(int i = 0; i < countMr;) {
            int nodeid = colM_r[i];
            int size   = colM_r[i+1];
            i += 2;

            // find whether I am owner of this node
					if(NodePID(nodeid) == lComm()->MyPID()) {
              // get the node
			  Teuchos::RCP<MOERTEL::Node> snode = GetNodeView(nodeid);

						if(snode == Teuchos::null) {
							std::stringstream oss;
							oss << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
                     << "***ERR*** Cannot find view of node " << nodeid << endl
                     << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
							throw ReportError(oss);
              }
              
              // get the lagrange multipliers
              int nslmdof = snode->Nlmdof();
              const int* slmdof = snode->LMDof();
              
              // loop colM_r/valM_r and assemble
						for(int j = 0; j < size; ++j) {
                int colmnode = colM_r[i+j];
                double val   = valM_r[i+j];

							if(abs(val) < CONSTRAINT_MATRIX_ZERO)
                  continue;
                
                // get view of column node and primal dofs
				Teuchos::RCP<MOERTEL::Node> colnode = GetNodeView(colmnode);

							if(colnode == Teuchos::null) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
                       << "***ERR*** Cannot find view of node " << colmnode << endl
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }

                int nmdof = colnode->Ndof();
                const int* mdof = colnode->Dof();

							if(nmdof != nslmdof) {
								std::stringstream oss;
								oss << "***ERR*** MOERTEL::Interface::AssembleJFNKVec:\n"
                       << "***ERR*** Mismatch in # primal dofs and lagrange mutlipliers\n"
                       << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
								throw ReportError(oss);
                }
                
							for(int k = 0; k < nslmdof; ++k) {
                  int row = slmdof[k];
                  int col = mdof[k];
                  //cout << "Proc " << lComm()->MyPID() << " inserting M row/col:" << row << "/" << col << " val " << val << endl;
				  
		  // Assemble M times soln
		  // Row of M determines row in rhs
		  // col of M determines row of soln

		  sel->AssembleNodeVal(row, col, val);

		  /*
                  int err = M.SumIntoGlobalValues(row,1,&val,&col);
                  if (err)
                    err = M.InsertGlobalValues(row,1,&val,&col);
					*/

                } // for (int k=0; k<nslmdof; ++k)
              } // for (int j=0; j<size; ++j)

              i += size;
            }

            else // I am not owner of this node, skip it
              i += size;
          } // for (int i=0; i<countMr;)
        } // if (proc!=lComm()->MyPID())
        
        colD_r.clear();
        valD_r.clear();
        colM_r.clear();
        valM_r.clear();
      } // for (int proc=0; proc<lComm()->NumProc(); ++proc)
  
    colD_s.clear();
    valD_s.clear();
    colM_s.clear();
    valM_s.clear();
  } // if (lComm()->NumProc()!=1)
  
#ifdef PDANDM
	Dmat.Print(std::cout);
	Mmat.Print(std::cout);
	throw "Done";
#endif

  return true;
}
