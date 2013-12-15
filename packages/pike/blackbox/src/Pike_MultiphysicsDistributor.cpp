#ifndef PIKE_MULTIPHYSICS_DISTRIBUTOR_IMPL_HPP
#define PIKE_MULTIPHYSICS_DISTRIBUTOR_IMPL_HPP

#include "mpi.h"
#include <iostream>
#include <algorithm>
#include <set>

#include "Pike_MultiphysicsDistributor.hpp"

namespace pike {

  MultiphysicsDistributor::MultiphysicsDistributor(const std::string distributorName) :
    myName_(distributorName),
    applicationAddMethodType_(Unknown),
    setupCalled_(false)
  { }
  
  void 
  MultiphysicsDistributor::addApplication(const ApplicationIndex index, const std::string name, const int numberOfProcesses)
  {
    if (applicationAddMethodType_ == Unknown)
      applicationAddMethodType_ = NumberOfProcesses;

    TEUCHOS_ASSERT(applicationAddMethodType_ == NumberOfProcesses);

    TEUCHOS_TEST_FOR_EXCEPTION( (applications_.find(index) != applications_.end()),std::logic_error,
				"Duplicate Index Error: The index for the application with index " << index << " named \"" << name << "\" has already been used!");

    TEUCHOS_TEST_FOR_EXCEPTION( (applicationNameToIndex_.find(name) != applicationNameToIndex_.end()),std::logic_error,
				"Duplicate Name Error: The name for the application with index " << index << " named \"" << name << "\" has already been used!");

    applicationNameToIndex_[name] = index;
    appRegistrationOrder_.push_back(index);

    ApplicationData app;
    app.name = name;
    app.numProcesses = numberOfProcesses;

    applications_[index] = app;
  }

  void 
  MultiphysicsDistributor::addApplication(const ApplicationIndex index, const std::string name, const std::vector<int> processes)
  {
    if (applicationAddMethodType_ == Unknown)
      applicationAddMethodType_ = VectorOfProcessIndices;

    TEUCHOS_ASSERT(applicationAddMethodType_ != NumberOfProcesses);

    TEUCHOS_TEST_FOR_EXCEPTION( (applications_.find(index) != applications_.end()),std::logic_error,
				"Duplicate Index Error: The index for the application with index " << index << " named \"" << name << "\" has already been used!");

    TEUCHOS_TEST_FOR_EXCEPTION( (applicationNameToIndex_.find(name) != applicationNameToIndex_.end()),std::logic_error,
				"Duplicate Name Error: The name for the application with index " << index << " named \"" << name << "\" has already been used!");

    applicationNameToIndex_[name] = index;
    appRegistrationOrder_.push_back(index);

    ApplicationData app;
    app.name = name;

    // make sure processes are unique
    std::set<int> tmp_processes;
    for (typename std::vector<int>::const_iterator i=processes.begin(); i != processes.end(); ++i)
      tmp_processes.insert(*i);

    app.processes.assign(tmp_processes.begin(),tmp_processes.end());
    app.startProcess = app.processes[0];

    applications_[index] = app;
  }
  
  void
  MultiphysicsDistributor::addApplication(const ApplicationIndex index, const std::string name, const int beginRank, const int endRank)
  {
    if (applicationAddMethodType_ == Unknown)
      applicationAddMethodType_ = VectorOfProcessIndices;

    TEUCHOS_ASSERT(applicationAddMethodType_ != NumberOfProcesses);

    TEUCHOS_TEST_FOR_EXCEPTION( (applications_.find(index) != applications_.end()),std::logic_error,
				"Duplicate Index Error: The index for the APPLICATION with index " << index << " named \"" << name << "\" has already been used!");

    TEUCHOS_TEST_FOR_EXCEPTION( (applicationNameToIndex_.find(name) != applicationNameToIndex_.end()),std::logic_error,
				"Duplicate Name Error: The name for the APPLICATION with index " << index << " named \"" << name << "\" has already been used!");

    applicationNameToIndex_[name] = index;
    appRegistrationOrder_.push_back(index);

    ApplicationData app;
    app.name = name;
    TEUCHOS_ASSERT(beginRank <= endRank);
    for (int p = beginRank; p < endRank+1; ++p) 
      app.processes.push_back(p);
    app.startProcess = beginRank;
    applications_[index] = app;
  }
  
  void
  MultiphysicsDistributor::addTransfer(const TransferIndex index, const ApplicationIndex a, const ApplicationIndex b, const std::string name)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (transfers_.find(index) != transfers_.end()), std::logic_error,
				"Duplicate Index Error: The index for the TRANSFER with index " << index << " named \"" << name << "\" has already been used!");

    TEUCHOS_TEST_FOR_EXCEPTION( (transferNameToIndex_.find(name) != transferNameToIndex_.end()),std::logic_error,
				"Duplicate Name Error: The name for the TRANSFER with index " << index << " named \"" << name << "\" has already been used!");

    transfers_[index].push_back(a);
    transfers_[index].push_back(b);

    transferNames_[index] = name;
    transferNameToIndex_[name] = index;
  }

  void
  MultiphysicsDistributor::addTransfer(const TransferIndex index, const std::vector<ApplicationIndex>& appIndices, const std::string name)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (transfers_.find(index) != transfers_.end()), std::logic_error,
				"Duplicate Index Error: The index for the TRANSFER with index " << index << " named \"" << name << "\" has already been used!");

    TEUCHOS_TEST_FOR_EXCEPTION( (transferNameToIndex_.find(name) != transferNameToIndex_.end()),std::logic_error,
				"Duplicate Name Error: The name for the TRANSFER with index " << index << " named \"" << name << "\" has already been used!");

    transfers_[index] = appIndices;

    transferNames_[index] = name;
    transferNameToIndex_[name] = index;
  }

  void
  MultiphysicsDistributor::setup(const Teuchos::RCP<const Teuchos::Comm<int> >& globalComm, bool printCommDistribution)
  {
    TEUCHOS_ASSERT(!setupCalled_);
    TEUCHOS_ASSERT(nonnull(globalComm));

    globalComm_ = globalComm;

    TEUCHOS_TEST_FOR_EXCEPTION((applications_.size() < 1),std::logic_error,"Error: No apps were registered with the distributor!");
    
    this->buildComms();

    this->buildOStreams();

    //*************************************
    // Print the comm details
    //*************************************
    if (printCommDistribution) {

      for (typename std::vector<ApplicationIndex>::const_iterator app = appRegistrationOrder_.begin(); app != appRegistrationOrder_.end(); ++app) {
	const std::vector<int>& appRanks = applications_[*app].processes; 
	*pout_ << "Application Ranks(" << applications_[*app].name << ") = ";
	for (std::vector<int>::const_iterator r=appRanks.begin(); r != appRanks.end(); ++r) {
	  if (r != appRanks.begin())
	    *pout_ << ",";
	  *pout_ << *r;
	}
	*pout_ << std::endl;
      }
      
      for (typename std::map<TransferIndex,std::vector<int> >::const_iterator t = transferRanks_.begin(); t != transferRanks_.end(); ++t) {
	const std::vector<int>& transRanks = t->second;
	*pout_ << "Transfer Ranks(" << transferNames_[t->first] << ") = ";
	for (std::vector<int>::const_iterator r=transRanks.begin(); r != transRanks.end(); ++r) {
	  if (r != transRanks.begin())
	    *pout_ << ",";
	  *pout_ << *r;
	}
	*pout_ << std::endl;
      }
      
    }

    setupCalled_ = true;
  }

  void
  MultiphysicsDistributor::buildComms()
  {
    TEUCHOS_ASSERT(applications_.size() == appRegistrationOrder_.size());

    //std::map<ApplicationIndex,std::pair<int,int> > appStartEndProcesses;

    if (applicationAddMethodType_ == NumberOfProcesses) {

      // Make sure the total number of processes equals global comm size
      int totalNumProcesses = 0;
      for (typename std::vector<ApplicationIndex>::const_iterator app = appRegistrationOrder_.begin();
	   app != appRegistrationOrder_.end(); ++app) {
	TEUCHOS_TEST_FOR_EXCEPTION((applications_[*app].numProcesses < 1),
				   std::logic_error,
				   "Error: the application \""<< applications_[*app].name 
				   << "\" requires at least one process allocated.  It was registered with " << applications_[*app].numProcesses << ".");
	totalNumProcesses += applications_[*app].numProcesses;
      }
      
      TEUCHOS_ASSERT(totalNumProcesses == globalComm_->getSize());
      
      // Add application process ranks to processes array
      int globalProcessCount = 0;
      for (typename std::vector<ApplicationIndex>::const_iterator app = appRegistrationOrder_.begin();
	   app != appRegistrationOrder_.end(); ++app) {
	for (int i = 0; i < applications_[*app].numProcesses; ++i,++globalProcessCount) {
	  if (i == 0)
	    applications_[*app].startProcess = globalProcessCount;
	  applications_[*app].processes.push_back(globalProcessCount);
	}
      }
    }
    
    // Build the application comms
    for (typename std::vector<ApplicationIndex>::const_iterator app = appRegistrationOrder_.begin();
	 app != appRegistrationOrder_.end(); ++app) {
      applicationComms_[*app] = globalComm_->createSubcommunicator(Teuchos::arrayViewFromVector(applications_[*app].processes));
    }

    // Build transfer comms
    for (typename std::map<TransferIndex,std::vector<ApplicationIndex> >::const_iterator t = transfers_.begin(); t != transfers_.end(); ++t) {
      for (typename std::vector<ApplicationIndex>::const_iterator a = t->second.begin(); a != t->second.end(); ++a) {
	transferRanks_[t->first].insert(transferRanks_[t->first].end(), applications_[*a].processes.begin(), applications_[*a].processes.end());
      }
      std::sort(transferRanks_[t->first].begin(), transferRanks_[t->first].end());
      transferRanks_[t->first].erase( std::unique( transferRanks_[t->first].begin(), transferRanks_[t->first].end() ), transferRanks_[t->first].end() );
    }

    for (typename std::map<TransferIndex,std::vector<int> >::const_iterator t = transferRanks_.begin(); t != transferRanks_.end(); ++t) {
	transferComms_[t->first] = globalComm_->createSubcommunicator(Teuchos::arrayViewFromVector(t->second));
    }


  }

  void
  MultiphysicsDistributor::buildOStreams()
  {

    {
      out_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      out_->setOutputToRootOnly(this->mapRankToCommWorldRank(0,*globalComm_));
      out_->pushLinePrefix(myName_);
      out_->setShowLinePrefix(true);
    }

    for (typename std::map<ApplicationIndex,ApplicationData>::const_iterator a = applications_.begin(); a != applications_.end(); ++a) {
      aout_[a->first] = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      aout_[a->first]->pushLinePrefix(a->second.name);
      aout_[a->first]->setOutputToRootOnly(this->mapRankToCommWorldRank(a->second.startProcess,*globalComm_));
      aout_[a->first]->setShowLinePrefix(true);
      aout_[a->first]->setShowProcRank(true);
    }

    {
      pout_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      for (typename std::map<ApplicationIndex,ApplicationData>::const_iterator a = applications_.begin(); a != applications_.end(); ++a) {
	if (this->appExistsOnProcess(a->first)) {
	  pout_->pushLinePrefix(a->second.name);
	}
      }
      pout_->setShowLinePrefix(true);
      pout_->setShowProcRank(true);
    }

  }

  int MultiphysicsDistributor::mapRankToCommWorldRank(int myRank, const Teuchos::Comm<int>& myTeuchosComm)
  {
    int worldRank = -1;
    
    const Teuchos::MpiComm<int>* myTeuchosMpiComm = dynamic_cast<const Teuchos::MpiComm<int>* >(&myTeuchosComm);
    TEUCHOS_ASSERT(myTeuchosMpiComm != 0);

    MPI_Comm myRawMpiComm = (*myTeuchosMpiComm->getRawMpiComm())();
    MPI_Group myGroup;
    MPI_Comm_group(myRawMpiComm,&myGroup);

    MPI_Group worldGroup;
    MPI_Comm_group(MPI_COMM_WORLD,&worldGroup);

    MPI_Group_translate_ranks(myGroup,1,&myRank,worldGroup,&worldRank);
   
    return worldRank;
  }

  std::string MultiphysicsDistributor::getApplicationName(const ApplicationIndex appIndex) const
  {
    typename std::map<ApplicationIndex,ApplicationData>::const_iterator n = applications_.find(appIndex);
    TEUCHOS_ASSERT(n != applications_.end());
    return n->second.name;
  }
  
  std::string MultiphysicsDistributor::getTransferName(const TransferIndex transferIndex) const
  {
    typename std::map<TransferIndex,std::string>::const_iterator n = transferNames_.find(transferIndex);
    TEUCHOS_ASSERT(n != transferNames_.end());
    return n->second;
  }

  MultiphysicsDistributor::ApplicationIndex MultiphysicsDistributor::getApplicationIndex(const std::string appName) const
  {
    typename std::map<std::string,ApplicationIndex>::const_iterator n = applicationNameToIndex_.find(appName);
    TEUCHOS_ASSERT(n != applicationNameToIndex_.end());
    return n->second;
  }
  
  MultiphysicsDistributor::TransferIndex MultiphysicsDistributor::getTransferIndex(const std::string transferName) const
  {
    typename std::map<std::string,MultiphysicsDistributor::TransferIndex>::const_iterator n = transferNameToIndex_.find(transferName);
    TEUCHOS_ASSERT(n != transferNameToIndex_.end());
    return n->second;
  }

  bool MultiphysicsDistributor::appExistsOnProcess(const MultiphysicsDistributor::ApplicationIndex index) const
  {
    typename std::map<ApplicationIndex,Teuchos::RCP<const Teuchos::Comm<int> > >::const_iterator comm = applicationComms_.find(index);
    TEUCHOS_ASSERT(comm != applicationComms_.end());
    return nonnull(comm->second);
  }

  bool MultiphysicsDistributor::transferExistsOnProcess(const MultiphysicsDistributor::TransferIndex index) const
  {
    typename std::map<TransferIndex,Teuchos::RCP<const Teuchos::Comm<int> > >::const_iterator comm = transferComms_.find(index);
    TEUCHOS_ASSERT(comm != transferComms_.end());
    return nonnull(comm->second);
  }

  Teuchos::RCP<const Teuchos::Comm<int> > MultiphysicsDistributor::getGlobalComm() const
  {
    return globalComm_;
  }

  Teuchos::RCP<const Teuchos::Comm<int> >
  MultiphysicsDistributor::getAppComm(const ApplicationIndex applicationIndex) const
  {
    typename std::map<ApplicationIndex,Teuchos::RCP<const Teuchos::Comm<int> > >::const_iterator ac =
      applicationComms_.find(applicationIndex);
    TEUCHOS_ASSERT(ac != applicationComms_.end())
    return ac->second;
  }

  Teuchos::RCP<const Teuchos::Comm<int> >
  MultiphysicsDistributor::getTransferComm(const TransferIndex transferIndex) const
  {
    typename std::map<TransferIndex,Teuchos::RCP<const Teuchos::Comm<int> > >::const_iterator tc = 
      transferComms_.find(transferIndex);
    return tc->second;
  }

  int MultiphysicsDistributor::getPrintRank(const ApplicationIndex index) const
  {
    typename std::map<ApplicationIndex,ApplicationData>::const_iterator s = applications_.find(index);
    TEUCHOS_ASSERT(s != applications_.end());
    return s->second.startProcess;
  }

  Teuchos::RCP<Teuchos::FancyOStream> MultiphysicsDistributor::getSerialOStream() const
  { return out_; }

  Teuchos::RCP<Teuchos::FancyOStream> MultiphysicsDistributor::getApplicationOStream(const ApplicationIndex applicationIndex) const
  { return aout_.find(applicationIndex)->second; }

  Teuchos::RCP<Teuchos::FancyOStream> MultiphysicsDistributor::getParallelOStream() const
  { return pout_; }

}

#endif
