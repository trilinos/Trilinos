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
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# 1. Redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the Corporation nor the names of the
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
# EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
# PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "Moertel_ManagerT.hpp"
#include "Teuchos_Time.hpp"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
MoertelT::ManagerT<ST, LO, GO, N>::ManagerT(const Teuchos::RCP<const Teuchos::Comm<LO> >& comm, int outlevel) :
	outlevel_(outlevel),
	comm_(comm),
	dimensiontype_(MoertelT::ManagerT<ST, LO, GO, N>::manager_none),
	problemmap_(Teuchos::null),
	inputmatrix_(Teuchos::null),
	constraintsmap_(Teuchos::null),
	D_(Teuchos::null),
	M_(Teuchos::null),
	saddlemap_(Teuchos::null),
	saddlematrix_(Teuchos::null),
	spdmatrix_(Teuchos::null),
	spdrhs_(Teuchos::null),
	I_(Teuchos::null),
	WT_(Teuchos::null),
	B_(Teuchos::null),
	annmap_(Teuchos::null),
	integrationparams_(Teuchos::null),
	solverparams_(Teuchos::null),
	solver_(Teuchos::null)
{
  // create default integration parameters
  Default_Parameters();
}

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 03/06|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
MoertelT::ManagerT<ST, LO, GO, N>::ManagerT(const Teuchos::RCP<const Teuchos::Comm<LO> >& comm, 
                          MoertelT::ManagerT<ST, LO, GO, N>::DimensionType dimension, 
                          int outlevel) :
outlevel_(outlevel),
comm_(comm),
dimensiontype_(dimension),
problemmap_(Teuchos::null),
inputmatrix_(Teuchos::null),
constraintsmap_(Teuchos::null),
D_(Teuchos::null),
M_(Teuchos::null),
saddlemap_(Teuchos::null),
saddlematrix_(Teuchos::null),
spdmatrix_(Teuchos::null),
spdrhs_(Teuchos::null),
I_(Teuchos::null),
WT_(Teuchos::null),
B_(Teuchos::null),
annmap_(Teuchos::null),
integrationparams_(Teuchos::null),
solverparams_(Teuchos::null),
solver_(Teuchos::null)
{
  // create default integration parameters
  Default_Parameters();
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
MoertelT::ManagerT<ST, LO, GO, N>::~ManagerT()
{
  interface_.clear();
}

/*----------------------------------------------------------------------*
 |  Add an interface (public)                                mwgee 06/05|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::ParameterList& 
MoertelT::ManagerT<ST, LO, GO, N>::Default_Parameters()
{
  if (integrationparams_==Teuchos::null)
  {
    integrationparams_ = rcp(new Teuchos::ParameterList());
    // see comment in mrtr_manager.H
    integrationparams_->set("exact values at gauss points",true);
    // 1D interface possible values are 1,2,3,4,5,6,7,8,10
    integrationparams_->set("number gaussian points 1D",3);
    // 2D interface possible values are 3,6,12,13,16,19,27
    integrationparams_->set("number gaussian points 2D",6);
  }
  return (*integrationparams_.get());
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
std::ostream& operator << (std::ostream& os, const MoertelT::ManagerT<ST, LO, GO, N>& man)
{ 
  man.Print();
  return (os);
}

/*----------------------------------------------------------------------*
 |  print all data                                           mwgee 06/05|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Print() const
{ 
  comm_->barrier();

  if (comm_->getRank()==0)
  std::cout << "\n========================= Mortar Manager =========================\n\n";
  fflush(stdout);
  
  comm_->barrier();

  typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::const_iterator curr;
  for (curr=interface_.begin(); curr!=interface_.end(); ++curr)
  {
	Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > inter = curr->second;
    if (inter==Teuchos::null)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::Print:\n"
           << "***ERR*** found NULL entry in map of interfaces\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    std::cout << *inter;
  }
  comm_->barrier();
  if (problemmap_ != Teuchos::null)
  {
    if (comm_->getRank() == 0)
    std::cout << "\n------------------ Input RowMap of original Problem ---------------\n";
    comm_->barrier();
    std::cout << *problemmap_;
  }
  comm_->barrier();
  if (constraintsmap_ != Teuchos::null)
  {
    if (comm_->getRank() == 0)
    std::cout << "\n------------------ RowMap of Constraints ---------------\n";
    comm_->barrier();
    std::cout << *constraintsmap_;
  }
  comm_->barrier();
  if (D_ != Teuchos::null)
  {
    if (comm_->getRank() == 0)
    std::cout << "\n------------------ Coupling Matrix D ---------------\n";
    comm_->barrier();
    std::cout << *D_;
  }
  comm_->barrier();
  if (M_ != Teuchos::null)
  {
    if (comm_->getRank() == 0)
    std::cout << "\n------------------ Coupling Matrix M ---------------\n";
    comm_->barrier();
    std::cout << *M_;
  }
  comm_->barrier();

  if (comm_->getRank()==0)
  std::cout << "\n========================= End Mortar Manager =========================\n\n";
  fflush(stdout);

  return true;
}

/*----------------------------------------------------------------------*
 |  Add an interface (public)                                mwgee 06/05|
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::AddInterface(const Teuchos::RCP<const MoertelT::InterfaceT<ST, LO, GO, N> >& interface) 
{
  if (!interface->IsComplete())
  {
    std::cout << "***ERR*** MoertelT::ManagerT::AddInterface:\n"
         << "***ERR*** Cannot add segment as Complete() was not called before\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > tmp = 
      Teuchos::rcp(new MoertelT::InterfaceT<ST, LO, GO, N>(*interface));
  interface_.insert(std::pair<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >(tmp->Id(),tmp));
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Choose dofs for lagrange multipliers (private)           mwgee 07/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
Teuchos::RCP<Tpetra::Map<LO, GO, N> > 
MoertelT::ManagerT<ST, LO, GO, N> ::LagrangeMultiplierDofs()
{
  if (problemmap_==Teuchos::null)
  {
    std::cout << "***ERR*** MoertelT::ManagerT::LagrangeMultiplierDofs:\n"
         << "***ERR*** problemmap_==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return Teuchos::null;
  }
  
  // find the largest row number in problemmap_
//  int maxinputGID = problemmap_->MaxAllGID();
  int maxinputGID = problemmap_->getMaxAllGlobalIndex();

  // start numbering lagrange multipliers from maxinputGID+1
  int minLMGID = maxinputGID+1;
  int maxLMGID = 0;
  
  int length = 0;
  
  // 2D problem:
  // loop interfaces and set LM dofs on the slave side for all slave nodes
  // that have a projection
  // 3D problem:
  // loop interfaces and set LM dofs on the slave side for all slave nodes
  // that have an integrated row of D and M projection

  // start with minLMGID and return maxLMGID+1 on a specific interface
  // Note this is collective for ALL procs
  typename std::map<int, Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    length -= minLMGID;
    maxLMGID = curr->second->SetLMDofs(minLMGID);
    if (!maxLMGID && maxLMGID!=minLMGID)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::LagrangeMultiplierDofs:\n"
           << "***ERR*** interface " << curr->second->Id() << " returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return Teuchos::null;
    }
    minLMGID = maxLMGID;
    length += maxLMGID;
  }  

  std::vector<GO> mylmids(length);
  int count=0;
  
  // get a vector of LM dof ids from each proc and each interface
  // and add it to the global one
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
	Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > inter = curr->second;
	std::vector<GO>* lmids = inter->MyLMIds();
    if (count+lmids->size() > mylmids.size())
      mylmids.resize(mylmids.size()+5*lmids->size());
    for (int i=0; i<(int)lmids->size(); ++i)
      mylmids[count++] = (*lmids)[i];
    delete lmids;
  }  
  mylmids.resize(count);
  int lsize = count;
  int gsize = 0;
  Teuchos::reduceAll(*comm_, Teuchos::REDUCE_SUM, 1, &lsize, &gsize);
  
  // create the rowmap for the constraints
  // Note that this map contains the global communicator from the MoertelT::Manager
  // NOT any interface local one

  Teuchos::RCP<Tpetra::Map<LO, GO, N> > map;

  if(lsize == 0) // this PE contributes nothing to the global map

	map = Teuchos::rcp(new Tpetra::Map<LO, GO, N> (gsize,NULL,0,0,comm_));

  else

	map = Teuchos::rcp(new Tpetra::Map<LO, GO, N> (gsize,&(mylmids[0]),lsize,0,comm_));

  // tidy up
  mylmids.clear();
  
  return map;
}

/*----------------------------------------------------------------------*
 |  integrate all mortar interfaces                       (public) 11/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Mortar_Integrate()
{
		bool status;
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() == MoertelT::ManagerT<ST, LO, GO, N>::manager_2D){

    status = Integrate_Interfaces_2D();

	if(!status) return status;

    status = Build_MD_2D();
  }
  else if (Dimension() == MoertelT::ManagerT<ST, LO, GO, N>::manager_3D){

    status = Integrate_Interfaces_3D();

	if(!status) return status;

    status = Build_MD_3D();
  }
  else
  {
    std::cout << "***ERR*** MoertelT::ManagerT::Mortar_Integrate:\n"
         << "***ERR*** problem dimension not set\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    status = false;
  }

  return status;
}

/*----------------------------------------------------------------------*
 |  Integrate all mortar interfaces                        (public) 8/10|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Integrate_Interfaces()
{
		bool status;
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() == MoertelT::ManagerT<ST, LO, GO, N>::manager_2D)
    status = Integrate_Interfaces_2D();
  else if (Dimension() == MoertelT::ManagerT<ST, LO, GO, N>::manager_3D)
    status = Integrate_Interfaces_3D();
  else
  {
    std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces:\n"
         << "***ERR*** problem dimension not set\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    status = false;
  }

  return status;
}

#if 0 // the old version
/*----------------------------------------------------------------------*
 |  integrate all mortar interfaces in 2D                 (public) 07/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
bool MoertelT::Manager::Mortar_Integrate_2D()
{
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() != MoertelT::Manager::manager_2D)
  {
    std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
         << "***ERR*** problem dimension is not 2D?????\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // check whether we have interfaces
  int ninter = Ninterfaces();
  if (!ninter)
    return true;

  //-------------------------------------------------------------------
  // check whether we have an input map
  if (problemmap_==Teuchos::null)
  {
    std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
         << "***ERR*** problemmap_==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // check whether we have a mortar side chosen on each interface or 
  // whether we have to chose it here
  {
    map<int,Teuchos::RCP<MoertelT::Interface> >::iterator curr;
    bool foundit = true;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int mside = curr->second->MortarSide();
      if (mside==-2)
      {
        foundit = false;
        break;
      }
    }
    if (!foundit) // we have to chose mortar sides ourself
      ChooseMortarSide();
  }  

  //-------------------------------------------------------------------
  // check whether functions have been set on interfaces
  // if not, check for functions flag and set them
  {
    bool foundit = true;
	std::map<int,Teuchos::RCP<MoertelT::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int nseg             = curr->second->GlobalNsegment();
      MOERTEL::Segment** segs = curr->second->GetSegmentView();
      for (int i=0; i<nseg; ++i)
        if (segs[i]->Nfunctions() < 2)
        {
          foundit = false;
          break;
        }
      delete [] segs;
      if (!foundit)
        curr->second->SetFunctionsFromFunctionTypes();
    }
  } 
  
  //-------------------------------------------------------------------
  // build projections for all interfaces
  {
	std::map<int,Teuchos::RCP<MoertelT::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok  = curr->second->Project();
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on projection\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }  

  //-------------------------------------------------------------------
  // this is probably the place to put detection of end segments
  // for each end segment, the order of the lagrange multiplier shape
  // function will be reduced by one
#if 0
  {
	std::map<int,Teuchos::RCP<MoertelT::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->DetectEndSegmentsandReduceOrder();
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false from DetectEndSegmentsandReduceOrder\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }
#endif

  //-------------------------------------------------------------------
  // choose dofs for lagrange multipliers and set them to slave nodes
  // build the rowmap for the coupling matrices M and D
  {
    constraintsmap_ = LagrangeMultiplierDofs();
    if (constraintsmap_==Teuchos::null)
    {
      std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
           << "***ERR*** LagrangeMultiplierDofs() returned NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
    
  //-------------------------------------------------------------------
  // build the map for the saddle point problem
  {
    bool ok = BuildSaddleMap();
    if (!ok)
    {
      std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
           << "***ERR*** BuildSaddleMap() returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }


  //-------------------------------------------------------------------
  // build the Tpetra_CrsMatrix D and M
  D_ = rcp(new Tpetra_CrsMatrix(saddlemap_,5));
  M_ = rcp(new Tpetra_CrsMatrix(saddlemap_,40));

  //-------------------------------------------------------------------
  // integrate all interfaces
  {
	std::map<int,Teuchos::RCP<MoertelT::Interface> >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Integrate(*D_,*M_);
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on integration\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }
  
  //-------------------------------------------------------------------
  // call FillComplete() on M_ and D_ 
  D_->fillComplete(saddlemap_, saddlemap_);
//  D_->OptimizeStorage();
  M_->fillComplete(saddlemap_, saddlemap_);
//  M_->OptimizeStorage();
  
  //-------------------------------------------------------------------
  // print this
  if (OutLevel()>9)
    std::cout << *this;

  return true;
}
#endif

/*----------------------------------------------------------------------*
 |  Build M and D matrices in 2D                          (public) 08/10|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Build_MD_2D()
{

  //-------------------------------------------------------------------
  // build the Tpetra_CrsMatrix D and M
  D_ = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(saddlemap_,5));
  M_ = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(saddlemap_,40));

  //-------------------------------------------------------------------
  // now that we have all maps and dofs we can assemble from the nodes
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Assemble(*D_,*M_);
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Build_MD_2D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on assembly\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // call FillComplete() on M_ and D_ 
  D_->fillComplete(saddlemap_, saddlemap_);
  D_ = MoertelT::StripZeros(*D_, 1.e-9);
//  D_->OptimizeStorage();
  M_->fillComplete(saddlemap_, saddlemap_);
//  M_->OptimizeStorage();
  
  //-------------------------------------------------------------------
  // print this
  if (OutLevel()>9)
    std::cout << *this;

  return true;
}

/*----------------------------------------------------------------------*
 |  just integrate all mortar interfaces in 2D                 (public) 07/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Integrate_Interfaces_2D()
{

  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() != MoertelT::ManagerT<ST, LO, GO, N>::manager_2D)
  {
    std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
         << "***ERR*** problem dimension is not 2D?????\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // check whether we have interfaces
  int ninter = Ninterfaces();
  if (!ninter)
    return true;

  //-------------------------------------------------------------------
  // check whether we have an input map
  if (problemmap_==Teuchos::null)
  {
    std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
         << "***ERR*** problemmap_==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // check whether we have a mortar side chosen on each interface or 
  // whether we have to chose it here
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    bool foundit = true;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int mside = curr->second->MortarSide();
      if (mside==-2)
      {
        foundit = false;
        break;
      }
    }
    if (!foundit) // we have to chose mortar sides ourself
      ChooseMortarSide();
  }  

  //-------------------------------------------------------------------
  // check whether functions have been set on interfaces
  // if not, check for functions flag and set them
  {
    bool foundit = true;
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int nseg             = curr->second->GlobalNsegment();
      MOERTEL::Segment** segs = curr->second->GetSegmentView();
      for (int i=0; i<nseg; ++i)
        if (segs[i]->Nfunctions() < 2)
        {
          foundit = false;
          break;
        }
      delete [] segs;
      if (!foundit)
        curr->second->SetFunctionsFromFunctionTypes();
    }
  } 
  
  //-------------------------------------------------------------------
  // build normals for all interfaces
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->BuildNormals();
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on building normals\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }  

  //-------------------------------------------------------------------
  // prepare the boundary modification for 2D interfaces
  // Nodes on the edge of an interface will not carry LMs so they
  // do not conflict with other interfaces
  // The choice of the Mortar side in the case of several interfaces
  // is then arbitrary
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->DetectEndSegmentsandReduceOrder();
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false from DetectEndSegmentsandReduceOrder\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // integrate all interfaces
  // NOTE: 2D and 3D integration differ:
  // the 3D integration works without choosing Lagrange multipliers
  // in advance. All integrated rows of D and M are stored as
  // scalars in the nodes
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Integrate_2D(integrationparams_);
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on integration\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // choose dofs for lagrange multipliers and set them to slave nodes
  // build the rowmap for the coupling matrices M and D
  {
    constraintsmap_ = LagrangeMultiplierDofs();
    if (constraintsmap_==Teuchos::null)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
           << "***ERR*** LagrangeMultiplierDofs() returned NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  
  //-------------------------------------------------------------------
  // build the map for the saddle point problem
  {
    bool ok = BuildSaddleMap();
    if (!ok)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_2D:\n"
           << "***ERR*** BuildSaddleMap() returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |  Build M & D 3D case                                   (public) 11/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Build_MD_3D()
{

  //-------------------------------------------------------------------
  // build the Tpetra_CrsMatrix D and M
  D_ = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(saddlemap_,50));
  M_ = Teuchos::rcp(new Tpetra::CrsMatrix<ST, LO, GO, N>(saddlemap_,100));

  //-------------------------------------------------------------------
  // now that we have all maps and dofs we can assemble from the nodes
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Assemble(*D_,*M_);
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Build_MD_3D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on assembly\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // call FillComplete() on M_ and D_ 
  D_->fillComplete(saddlemap_, saddlemap_);
  // GAH - found by kimliegeois
  // when remeshing, the size of elements can become so small that all the entries, or part 
  // of the meaning entries, are below the threshold and therefore dropped. 
  // D_ = Teuchos::rcp(MOERTEL::StripZeros(*D_,1.e-4));
//  D_->OptimizeStorage();
  M_->fillComplete(saddlemap_, saddlemap_);
//  M_->OptimizeStorage();
  
  //-------------------------------------------------------------------
  // print this
  if (OutLevel()>9)
    std::cout << *this;

  return true;
}

/*----------------------------------------------------------------------*
 |  just integrate all mortar interfaces 3D case               (public) 11/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::Integrate_Interfaces_3D()
{
  //-------------------------------------------------------------------
  // check for problem dimension
  if (Dimension() != MoertelT::ManagerT<ST, LO, GO, N>::manager_3D)
  {
    std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
         << "***ERR*** problem dimension is not 3D?????\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //-------------------------------------------------------------------
  // check whether we have interfaces
  int ninter = Ninterfaces();
  if (!ninter)
    return true;

  //-------------------------------------------------------------------
  // check whether we have an input map
  if (problemmap_==Teuchos::null)
  {
    std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
         << "***ERR*** problemmap_==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // check whether we have a mortar side chosen on each interface or 
  // whether we have to chose it here
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    bool foundit = true;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int mside = curr->second->MortarSide();
      if (mside==-2)
      {
        foundit = false;
        break;
      }
    }
    if (!foundit) // we have to chose mortar sides ourself
      ChooseMortarSide();
  }  

  //-------------------------------------------------------------------
  // check whether functions have been set on interfaces
  // if not, check for functions flags and set functions from them
  {
    bool foundit = true;
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      const int nseg          = curr->second->GlobalNsegment();
      MOERTEL::Segment** segs = curr->second->GetSegmentView();
      for (int i=0; i<nseg; ++i)
        if (segs[i]->Nfunctions() < 2)
        {
          foundit = false;
          break;
        }
      delete [] segs;
      if (!foundit)
        curr->second->SetFunctionsFromFunctionTypes();
    }
  } 
  
  //-------------------------------------------------------------------
  // build normals for all interfaces
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->BuildNormals();
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on building normals\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }  

  //-------------------------------------------------------------------
  // prepare the boundary modification for 3D interfaces
  // Nodes on the edge of an interface will not carry LMs so they
  // do not conflict with other interfaces
  // The choice of the Mortar side in the case of several interfaces
  // is then arbitrary
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->DetectEndSegmentsandReduceOrder();
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false from DetectEndSegmentsandReduceOrder\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // integrate all interfaces
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->Mortar_Integrate(integrationparams_);
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on integration\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  //-------------------------------------------------------------------
  // choose dofs for lagrange multipliers and set them to slave nodes
  // build the rowmap for the coupling matrices M and D
  {
    constraintsmap_ = LagrangeMultiplierDofs();
    if (constraintsmap_==Teuchos::null)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
           << "***ERR*** LagrangeMultiplierDofs() returned NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  
  //-------------------------------------------------------------------
  // build the map for the saddle point problem
  {
    bool ok = BuildSaddleMap();
    if (!ok)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::Integrate_Interfaces_3D:\n"
           << "***ERR*** BuildSaddleMap() returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  return true;
}

/*----------------------------------------------------------------------*
 | Assemble interface nodes times soln into JFNK residual vector
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::AssembleInterfacesIntoResidual(MOERTEL::Lmselector *sel)
{

  //-------------------------------------------------------------------
  // now that we have all maps and dofs we can assemble from the nodes
  {
	typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      bool ok = curr->second->AssembleJFNKVec(sel);
      if (!ok)
      {
        std::cout << "***ERR*** MoertelT::ManagerT::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on assembly\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
  }

  return true;
}


/*----------------------------------------------------------------------*
 |                                                                 09/05|
 |  choose the mortar side                                              |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::ChooseMortarSide()
{
  // find all interfaces
  std::vector<Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > > inter(Ninterfaces());
  typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
  curr=interface_.begin();
  int count = 0;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
      inter[count] = curr->second;
      ++count;
  }
  inter.resize(count);
  
  // call choice of the mortar side for all 1D interfaces
  bool ok = false; 
  if (Dimension() == MoertelT::ManagerT<ST, LO, GO, N>::manager_2D)
    ok = ChooseMortarSide_2D(inter);
  if (Dimension() == MoertelT::ManagerT<ST, LO, GO, N>::manager_3D)
    ok = ChooseMortarSide_3D(inter);

  // tidy up
  inter.clear();  

  return ok;
}


/*----------------------------------------------------------------------*
 |                                                                 12/05|
 |  choose the mortar side                                              |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::ChooseMortarSide_3D(std::vector<Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >& inter)
{
  // loop interfaces and choose the side with less nodes as slave side
  // (only if not already chosen on some interface)
  for (int i=0; i<(int)inter.size(); ++i)
  {
    if (inter[i]->MortarSide() == 0 || inter[i]->MortarSide() == 1) continue;
    if (inter[i]->GlobalNnode(0) <= inter[i]->GlobalNnode(1))
      inter[i]->SetMortarSide(1);
    else
      inter[i]->SetMortarSide(0);
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                                 12/05|
 |  choose the mortar side                                              |
 *----------------------------------------------------------------------*/
template <class ST,
          class LO,
          class GO,
          class N >
bool 
MoertelT::ManagerT<ST, LO, GO, N>::ChooseMortarSide_2D(std::vector<Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >& inter)
{
  // loop interfaces and choose the side with less nodes as slave side
  // (only if not already chosen on some interface)
  for (int i=0; i<(int)inter.size(); ++i)
  {
    if (inter[i]->MortarSide() == 0 || inter[i]->MortarSide() == 1) continue;
    if (inter[i]->GlobalNnode(0) <= inter[i]->GlobalNnode(1))
      inter[i]->SetMortarSide(1);
    else
      inter[i]->SetMortarSide(0);
  }
  return true;
}

#if 0
/*----------------------------------------------------------------------*
 |                                                                 09/05|
 |  choose the mortar side for 1D interfaces                            |
 *----------------------------------------------------------------------*/
bool MoertelT::Manager::ChooseMortarSide_2D(std::vector<Teuchos::RCP<MoertelT::Interface> >& inter)
{
  // number of interfaces
  const int ninter = inter.size();
  if (ninter < 2) 
  {
    if (inter[0]->MortarSide() != 0 && inter[0]->MortarSide() != 1)
      inter[0]->SetMortarSide(0);
    return true;
  }
  
  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
    {
      std::cout << "---MoertelT::Manager: start interface coloring:\n";
      for (int i=0; i<ninter; ++i)
        std::cout << "Interface " << inter[i]->Id() << " : Mortar Side : " << inter[i]->MortarSide() << std::endl;
      std::cout << "---MoertelT::Manager: finding common nodes on interfaces:\n";
    }
    fflush(stdout);
  }
  Comm().Barrier();
  
  // time this process
  Teuchos_Time time(Comm());
  time.ResetStartTime();
  
  
  // get a view of all nodes from all interfaces
  std::vector<MOERTEL::Node**> nodes(ninter);
  std::vector<int>          nnodes(ninter);
  for (int i=0; i<ninter; ++i)
  {
    nnodes[i] = inter[i]->GlobalNnode(); // returns 0 for procs not in lComm()
    nodes[i]  = inter[i]->GetNodeView(); // get vector of ALL nodes on inter i
  }
  
  // loop all interfaces
  for (int i=0; i<ninter; ++i)
  {
    // loop all nodes on that interface 
    // (procs not part of inter[i]->lComm() don't loop because they have nnodes[i]=0)
    for (int j=0; j<nnodes[i]; ++j)
    {
      // do nothing for node that does not belong to me
      if (inter[i]->NodePID(nodes[i][j]->Id()) != inter[i]->lComm()->MyPID())
        continue;
        
      // do nothing for a node that has been flagged cornernode before
      if (nodes[i][j]->IsCorner())
        continue;
        
      // the node we are currently looking for
      int actnodeid       = nodes[i][j]->Id();
      MOERTEL::Node* inode   = nodes[i][j];
      
      // search all other interfaces for this node
      for (int k=0; k<ninter; ++k)
      {
        if (k==i) 
          continue; // don't do the current interface
        
        MOERTEL::Node* knode = NULL;
        for (int l=0; l<nnodes[k]; ++l)
          if (actnodeid==nodes[k][l]->Id())
          {
            knode = nodes[k][l];
            break;
          }
        if (!knode) // node inode is not on interface k 
          continue;
          
        // found node actnodeid on interface i and k
        if (OutLevel()>9)
        {
          std::cout << "Node " << actnodeid << " on interfaces " << inter[i]->Id() << " and " << inter[k]->Id() << std::endl;
          fflush(stdout);
        }
        
        // flag that node on interfaces i and k as cornernode
        inode->SetCorner();
        knode->SetCorner();
      } // for (int k=0; k<ninter; ++k)
    } // for (int j=0; j< nnodes[i]; ++j)
  } // for (int i=0; i<ninter; ++i)

  
  // make the cornernode information redundant
    // loop all interfaces
  for (int i=0; i<ninter; ++i)
  {
    // if I'm not part of inter[i], continue
    if (!inter[i]->lComm())
      continue;
    for (int proc=0; proc<inter[i]->lComm()->NumProc(); ++proc)
    {
      // get the corner nodes I have found
	  std::vector<int> bcast(4);
      int bsize = 0;
      if (proc==inter[i]->lComm()->MyPID())
        for (int j=0; j<nnodes[i]; ++j)
        {
          if (nodes[i][j]->IsCorner())
          {
            if ((int)bcast.size() <= bsize)
              bcast.resize(bsize+10);
            bcast[bsize] = j;
            ++bsize;  
          }
        } 
      inter[i]->lComm()->Broadcast(&bsize,1,proc);
      bcast.resize(bsize);
      inter[i]->lComm()->Broadcast(&(bcast[0]),bsize,proc);
      if (proc!=inter[i]->lComm()->MyPID())
        for (int k=0; k<bsize; ++k)
        {
          int j = bcast[k];
          nodes[i][j]->SetCorner();
        }
      bcast.clear();
    } // for (int proc=0; proc<Comm()->NumProc(); ++proc)
  } // for (i=0; i<ninter; ++i)

  
  // we have list of interfaces inter[i], now we need
  // to build list of cornernodes associated with each interface
  std::vector< std::vector<int> > cornernodes(ninter);
  for (int i=0; i<ninter; ++i) 
  {
    int bprocs = 0; int bprocr = 0;
    int bsize = 0;
    if (inter[i]->lComm())
      if (inter[i]->lComm()->MyPID()==0)
      {
        bprocs = Comm().MyPID();
        for (int j=0; j<nnodes[i]; ++j)
          if (nodes[i][j]->IsCorner())
          {
            if ((int)cornernodes[i].size() <= bsize)
              cornernodes[i].resize(bsize+5);
            cornernodes[i][bsize] = nodes[i][j]->Id();
            ++bsize;
          }
      }
    Comm().MaxAll(&bprocs,&bprocr,1);
    Comm().Broadcast(&bsize,1,bprocr);
    cornernodes[i].resize(bsize);
    Comm().Broadcast(&(cornernodes[i][0]),bsize,bprocr);
  }

  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
      std::cout << "---MoertelT::Manager: finding common nodes: " 
           << time.ElapsedTime() << " sec" << std::endl;
    fflush(stdout);
  }
  time.ResetStartTime();
  
  // we now color the interfaces in such a way that every corner node
  // has lagrange multipliers either never or once but never more then once
  // this exercice is a little tricky as each interface has it's own
  // intra-communicator where a bunch of procs but not necessarily all of them
  // belong to. Therefore some searching has to be done by procs that are 
  // member of a specific intra-comm, while some communication has to be
  // done using the global comm.
  // This loop probably does not scale well in parallel but I'm already happy
  // if it works at all!
  // The coloring allows for some interfaces to have a user-specified mortar side
  // It tries to make use of such a prescribed choice of side. But
  // when specifying a mortar side on one or more interfaces but not on all of them, 
  // there is no guarantee that the coloring works out even if the user chose
  // the mortar sides carefully so the coloring should/could work out!
  // This is due to the fact that there is some heuristics in here (e.g. line 926)
  // and there is more then one valid way to color the problem.
  std::vector<int> flags((2*ninter));
  std::vector<int> flagr((2*ninter));
  for (int i=0; i<ninter; ++i)
  {
     // check whether inter[i] has mortar side assigned
     // if not, go through all cornernodes and check dependency on
     // other interfaces
     // if there is a dependency, choose mortar side of inter[i] approp., 
     // if there is no dependency, choose something

     if (inter[i]->MortarSide() == -2)
     {
       // create a vector of flags as follows:
       // flags[0..ninter-1]        flags for side 0 of inter[i]
       // flags[ninter..2*ninter-1] flags for side 1 of inter[i]
       // flags[k]==  0 : no dependency
       // flags[k]==  1 : dependency but no side chosen yet on inter[k]
       // flags[k]==  2 : side should be chosen master/mortar side
       // flags[k]==  3 : side should be chosen slave/LM side
       for (int k=0; k<2*ninter; ++k) flags[k] = 0;
       
       // actnode[0] holds the active node id
       // actnode[1] holds the side this node is on on inter[i]  
       int actnodes[2]; 
       int actnoder[2]; 
       // loop cornernodes on inter[i]
       for (int j=0; j<(int)cornernodes[i].size(); ++j)
       {
         actnodes[0] = 0; actnodes[1] = 0;    
         if (inter[i]->lComm())
           if (inter[i]->lComm()->MyPID()==0)
           {
             actnodes[0] = cornernodes[i][j];
             actnodes[1] = inter[i]->GetSide(actnodes[0]);
           }
         Comm().MaxAll(actnodes,actnoder,2);
         
         // loop all other interfaces and check for actnoder[0]
         for (int k=0; k<ninter; ++k)
         {
           if (k==i) continue; // don't do inter[i]
           if (inter[k]->lComm())
             if (inter[k]->lComm()->MyPID()==0)
             {
               for (int l=0; l<(int)cornernodes[k].size(); ++l)
               {
                 if (actnoder[0]==cornernodes[k][l])
                 {
                   flags[k]        = 1; // found a dependency
                   flags[ninter+k] = 1; // found a dependancy
                   int mside = inter[k]->MortarSide();
                   if (mside != -2) // side has been chosen before on inter[k]
                   {
                     int nodeside = inter[k]->GetSide(actnoder[0]);
                     // found node on master/mortar side
                     // so set flags such that node is on slave on inter[i]
                     if (nodeside == mside) 
                     {
                       // actnode[1] has to become slave side on inter[i]
                       int sside_i = actnoder[1];
                       int mside_i = inter[k]->OtherSide(sside_i);
                       flags[sside_i*ninter+k] = 3;
                       flags[mside_i*ninter+k] = 2;
                     }
                     // found node on slave/LM side
                     // so set flags that node has to be on master on inter[i]
                     else if (nodeside == inter[k]->OtherSide(mside))
                     {
                       int mside_i = actnoder[1];
                       int sside_i = inter[k]->OtherSide(mside_i);
                       flags[sside_i*ninter+k] = 3;
                       flags[mside_i*ninter+k] = 2;
                     }
                     else
                     {
						 std::stringstream oss;
                       oss << "***ERR*** MoertelT::Manager::ChooseMortarSide:\n"
                            << "***ERR*** weird\n"
                            << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
					   throw ReportError(oss);
                     }
                   }
                   break;
                 }
               } // for (int l=0; l<cornernodes[k].size(); ++l)
             }
         } // for (int k=0; k<ninter; ++k)
       } // for (int j=0; j<cornernodes[i].size(); ++j)

       // sum all flags
       Comm().MaxAll(&(flags[0]),&(flagr[0]),(2*ninter));

#if 0       
       std::cout << "Flags side 0:\n";
       for (int j=0; j<ninter; ++j) std::cout << "inter " << j << " flag " << flagr[j] << std::endl;
       std::cout << "Flags side 1:\n";
       for (int j=ninter; j<2*ninter; ++j) std::cout << "inter " << j << " flag " << flagr[j] << std::endl;
#endif
       
       // loop through flags and make a decision what to chose on inter[i]
       for (int k=0; k<ninter; ++k)
       {
         if (i==k) continue; 
         // I don't care what is chosen
         if (flagr[k]==0 && flagr[ninter+k]==0); 
         // i don't care either
         else if (flagr[k]==1 && flagr[ninter+k]==1);
         // side 0 has to be master, side 1 has to be slave
         else if (flagr[k]==2 && flagr[ninter+k]==3)
         {
           // check whether a side has been set before and this side is different?
           if (inter[i]->MortarSide()==1)
           {
             if (OutLevel()>3)
             std::cout << "***WRN*** MoertelT::Manager::ChooseMortarSide:\n"
                  << "***WRN*** interface " << inter[i]->Id() << ": want to set mortar side to 0 but has been set to 1 before\n"
                  << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             continue;
           }
           inter[i]->SetMortarSide(0);
         }
         // side 0 has to be slave, side 1 has to be master
         else if (flagr[k]==3 && flagr[ninter+k]==2)
         {
           // check whether a side has been set before and this side is differen?
           if (inter[i]->MortarSide()==0)
           {
             if (OutLevel()>3)
             std::cout << "***WRN*** MoertelT::Manager::ChooseMortarSide:\n"
                  << "***WRN*** interface " << inter[i]->Id() << ": want to set mortar side to 1 but has been set to 0 before\n"
                  << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             continue;
           }
           inter[i]->SetMortarSide(1);
         }
         else
         {
			 std::stringstream oss;
           oss << "***ERR*** MoertelT::Manager::ChooseMortarSide:\n"
                << "***ERR*** unknown type of coloring flag: " << flags[k] << "\n"
                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		   throw ReportError(oss);
         }
       }
       // Obviously we don't need to care what side inter[i] is going to have
       if (inter[i]->MortarSide()==-2)
         inter[i]->SetMortarSide(0);
         
       // now we loop through the flags again and color all interfaces,
       // that had a connection to inter[i] but did not care
       // (so they had flag[k]==1 and flag[ninter+k]==1)
       for (int k=0; k<ninter; ++k)
       {
         if (i==k) continue;
         if (flagr[k]==1 && flagr[ninter+k]==1)
         {
           if (inter[k]->MortarSide()!=-2)
           {
			 std::stringstream oss;
             oss << "***ERR*** MoertelT::Manager::ChooseMortarSide:\n"
                  << "***ERR*** weird, this interface is not supposed to already have a mortar side assigned\n"
                  << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
				throw ReportError(oss);
           }
           
           // loop through cornernodes of inter[i]
           for (int j=0; j<(int)cornernodes[i].size(); ++j)
           {
             actnodes[0] = 0; actnodes[1] = 0;    
             if (inter[i]->lComm())
               if (inter[i]->lComm()->MyPID()==0)
               {
                 actnodes[0] = cornernodes[i][j];
                 actnodes[1] = inter[i]->GetSide(actnodes[0]);
               }
             Comm().MaxAll(actnodes,actnoder,2);
             
             // intraproc 0 of inter[k] makes the decision on what to do
             int mortarside_s = -1;
             int mortarside_r;
             if (inter[k]->lComm())
               if (inter[k]->lComm()->MyPID()==0)
               {
                 // loop through inter[k]'s cornernode
                 for (int l=0; l<(int)cornernodes[k].size(); ++l)
                 {
                   if (actnoder[0]==cornernodes[k][l])
                   {
                     // get the mortar and nodal side of inter[i]
                     int mside_i   = inter[i]->MortarSide();
                     int nodeside_i = actnoder[1];
                     
                     // get the nodal side of inter[k]
                     int nodeside_k = inter[k]->GetSide(actnoder[0]);
                     
                     // the node is on the mortar side of i, 
                     // so put it on the slave side on k
                     if (mside_i == nodeside_i)
                     {
                       if (mortarside_s != -1)
                         if (mortarside_s != inter[k]->OtherSide(nodeside_k))
                         {
							std::stringstream oss;
                           oss << "***ERR*** MoertelT::Manager::ChooseMortarSide:\n"
                                << "***ERR*** interface has a conflict that can not be resolved\n"
                                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
							throw ReportError(oss);
                         }
                       mortarside_s = inter[k]->OtherSide(nodeside_k);
                     }
                     // the node is on the slave side of i,
                     // so put it on the mortar side of k
                     else 
                     {
                       if (mortarside_s != -1)
                         if (mortarside_s != nodeside_k)
                         {
							std::stringstream oss;
                           oss << "***ERR*** MoertelT::Manager::ChooseMortarSide:\n"
                                << "***ERR*** interface has a conflict that can not be resolved\n"
                                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
							throw ReportError(oss);
                         }
                       mortarside_s = nodeside_k;
                     }
                   }
                 } // for (l=0; l<cornernodes[k].size(); ++l)
               } // if (inter[k]->lComm()->MyPID()==0)
             Comm().MaxAll(&mortarside_s,&mortarside_r,1);
             if (mortarside_r != -1)
               inter[k]->SetMortarSide(mortarside_r);
           } // for (int j=0; j<cornernodes[i].size(); ++j)
           
         } // if (flagr[k]==1 && flagr[ninter+k]==1)
       } // for (k=0; k<ninter; ++k)

     } // if (inter[i]->MortarSide() == -2)
  } // for (int i=0; i<ninter; ++i)

  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
      std::cout << "---MoertelT::Manager: coloring interfaces : " 
           << time.ElapsedTime() << " sec" << std::endl;
    fflush(stdout);
  }

  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
    {
      std::cout << "---MoertelT::Manager: chosen mortar sides : \n";
      for (int i=0; i<ninter; ++i)
        std::cout << "Interface " << inter[i]->Id() << " : Mortar Side : " << inter[i]->MortarSide() << std::endl;
      std::cout << "---MoertelT::Manager: end interface coloring\n";
    }
    fflush(stdout);
  }
  Comm().Barrier();

  // tidy up
  nnodes.clear();
  for (int i=0; i<ninter; ++i) 
    if (nodes[i]) delete [] nodes[i];
  nodes.clear();
  for (int i=0; i<ninter; ++i) 
    cornernodes[i].clear();
  cornernodes.clear();
  flags.clear();
  flagr.clear();

  return true;
}
#endif

template <class ST,
          class LO,
          class GO,
          class N >
bool
MoertelT::ManagerT<ST, LO, GO, N>::BuildSaddleMap()
{
  // check whether all interfaces are complete and integrated
  typename std::map<int,Teuchos::RCP<MoertelT::InterfaceT<ST, LO, GO, N> > >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      std::cout << "***ERR*** MoertelT::ManagerT::BuildSaddleMap:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  // check whether we have a problemmap_
  if (problemmap_==Teuchos::null)
  {
      std::cout << "***ERR*** MoertelT::ManagerT::BuildSaddleMap:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }

  // check whether we have a constraintsmap_
  if (constraintsmap_==Teuchos::null)
  {
      std::cout << "***ERR*** MoertelT::ManagerT::BuildSaddleMap:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }

  // the saddle point problem rowmap is the problemmap_ + the constraintmap
  Tpetra::global_size_t numglobalelements = problemmap_->getGlobalNumElements() +
                          constraintsmap_->getGlobalNumElements();
  size_t nummyelements     = problemmap_->getNodeNumElements() +
                          constraintsmap_->getNodeNumElements();
  std::vector<GO> myglobalelements(nummyelements);
  size_t count = 0;
  auto inputmyglobalelements = problemmap_->getMyGlobalIndices();
  for (size_t i=0; i<problemmap_->getNodeNumElements(); ++i)
    myglobalelements[count++] = inputmyglobalelements[i];
  auto constraintsmyglobalelements = constraintsmap_->getMyGlobalIndices();
  for (size_t i=0; i<constraintsmap_->getNodeNumElements(); ++i)
    myglobalelements[count++] = constraintsmyglobalelements[i];
  if (count != nummyelements)
  {
    std::cout << "***ERR*** MoertelT::ManagerT::BuildSaddleMap:\n"
         << "***ERR*** Mismatch in dimensions\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  saddlemap_ = Teuchos::rcp(new Tpetra::Map<LO, GO, N>(numglobalelements, &(myglobalelements[0]), 
       nummyelements, 0, comm_));
  myglobalelements.clear();

  return true;
}

