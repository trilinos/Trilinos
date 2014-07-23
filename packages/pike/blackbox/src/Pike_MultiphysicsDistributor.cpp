#ifndef PIKE_MULTIPHYSICS_DISTRIBUTOR_IMPL_HPP
#define PIKE_MULTIPHYSICS_DISTRIBUTOR_IMPL_HPP

#include "mpi.h"
#include <iostream>
#include <algorithm>
#include <set>

#include "Pike_BlackBox_config.hpp"  // for debug define
#include "Pike_MultiphysicsDistributor.hpp"

namespace pike {
  int translateMpiRank(const int& rankA, 
		       const Teuchos::Comm<int>& commA,
		       const Teuchos::Comm<int>& commB)
  {
    MPI_Group groupA;
    {
      const Teuchos::MpiComm<int>* teuchosMpiCommA = dynamic_cast<const Teuchos::MpiComm<int>* >(&commA);
      TEUCHOS_ASSERT(teuchosMpiCommA != 0);
      MPI_Comm rawMpiCommA = (*teuchosMpiCommA->getRawMpiComm())();
      MPI_Comm_group(rawMpiCommA,&groupA);
    }
    
    MPI_Group groupB;
    {
      const Teuchos::MpiComm<int>* teuchosMpiCommB = dynamic_cast<const Teuchos::MpiComm<int>* >(&commB);
      TEUCHOS_ASSERT(teuchosMpiCommB != 0);
      MPI_Comm rawMpiCommB = (*teuchosMpiCommB->getRawMpiComm())();
      MPI_Comm_group(rawMpiCommB,&groupB);
    }
    
    int rankB = -1;
    MPI_Group_translate_ranks(groupA,1,const_cast<int*>(&rankA),groupB,&rankB);
    
    return rankB;
  }

  MultiphysicsDistributor::MultiphysicsDistributor(const std::string& distributorName) :
    myName_(distributorName),
    setupCalled_(false)
  { }

  MultiphysicsDistributor::ApplicationIndex
  MultiphysicsDistributor::addApplication(const std::string& name,
					  const std::vector<int>& processes)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (applicationNameToIndex_.find(name) != applicationNameToIndex_.end()),
				std::logic_error,
				"Duplicate Name Error: The application name \"" << name 
				<< "\" has already been used! Each application must have a unique name.");

    const ApplicationIndex index = applications_.size(); 
    applicationNameToIndex_[name] = index;

    ApplicationData app;
    app.name = name;

    // make sure processes are unique by using set
    std::set<int> tmp_processes;
    for (std::vector<int>::const_iterator i=processes.begin(); i != processes.end(); ++i)
      tmp_processes.insert(*i);

    app.processes.assign(tmp_processes.begin(),tmp_processes.end());
    app.startProcess = app.processes[0];

    applications_.push_back(app);

    return index;
  }

  MultiphysicsDistributor::ApplicationIndex
  MultiphysicsDistributor::addApplication(const std::string& name,
					  const int beginRank,
					  const int endRank)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (applicationNameToIndex_.find(name) != applicationNameToIndex_.end()),
				std::logic_error,
				"Duplicate Name Error: The application name \"" << name 
				<< "\" has already been used! Each application must have a unique name.");

    const ApplicationIndex index = applications_.size(); 
    applicationNameToIndex_[name] = index;

    ApplicationData app;
    app.name = name;
    TEUCHOS_ASSERT(beginRank <= endRank);
    for (int p = beginRank; p < endRank+1; ++p) 
      app.processes.push_back(p);
    app.startProcess = beginRank;
    applications_.push_back(app);
    return index;
  }

  MultiphysicsDistributor::TransferIndex
  MultiphysicsDistributor::addTransfer(const std::string& name,
				       const ApplicationIndex a,
				       const ApplicationIndex b)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (transferNameToIndex_.find(name) != transferNameToIndex_.end()),
				std::logic_error,
				"Duplicate Name Error: The transfer named \"" << name 
				<< "\" has already been used! Each transfer requires a unique name.");

    const TransferIndex index = transfers_.size();
    std::vector<ApplicationIndex> apps;
    apps.push_back(a);
    apps.push_back(b);
    transfers_.push_back(apps);
    transferNames_.push_back(name);
    transferNameToIndex_[name] = index;
    std::vector<int> dummy;
    transferRanks_.push_back(dummy); // Actual ranks determined during setup
    return index;
  }

  MultiphysicsDistributor::TransferIndex
  MultiphysicsDistributor::addTransfer(const std::string& name,
				       const std::vector<ApplicationIndex>& appIndices)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (transferNameToIndex_.find(name) != transferNameToIndex_.end()),
				std::logic_error,
				"Duplicate Name Error: The application name \"" << name 
				<< "\" has already been used! Each application must have a unique name.");

    const TransferIndex index = transfers_.size();
    transfers_.push_back(appIndices);
    transferNames_.push_back(name);
    transferNameToIndex_[name] = index;
    std::vector<int> dummy;
    transferRanks_.push_back(dummy); // Actual ranks determined during setup
    return index;
  }
  
  MultiphysicsDistributor::TransferIndex 
  MultiphysicsDistributor::addTransferByRanks(const std::string& name, const std::vector<int>& mpiRanks)
  {
    TEUCHOS_TEST_FOR_EXCEPTION( (transferNameToIndex_.find(name) != transferNameToIndex_.end()),
				std::logic_error,
				"Duplicate Name Error: The application name \"" << name 
				<< "\" has already been used! Each application must have a unique name.");
    

    const TransferIndex index = transfers_.size();
    std::vector<ApplicationIndex> dummy;
    transfers_.push_back(dummy);  // don't need application indices - ranks explicitly set
    transferNames_.push_back(name);
    transferNameToIndex_[name] = index;
    transferRanks_.push_back(mpiRanks);
    return index;
  }

  void
  MultiphysicsDistributor::setup(const Teuchos::RCP<const Teuchos::Comm<int> >& globalComm,
				 bool printCommDistribution)
  {
    TEUCHOS_ASSERT(!setupCalled_);
    TEUCHOS_ASSERT(nonnull(globalComm));

    // Duplicate the global comm to get a separate context for this
    // objects's communication
    globalComm_ = globalComm->duplicate();

    TEUCHOS_TEST_FOR_EXCEPTION((applications_.size() < 1),
			       std::logic_error,
			       "Error: No apps were registered with the distributor!");
    
    this->buildComms();

    this->buildOStreams();

    //*************************************
    // Print the comm details
    //*************************************
    if (printCommDistribution) {

      for (std::size_t app = 0; app < applications_.size(); ++app) {
	const std::vector<int>& appRanks = applications_[app].processes; 
	*pout_ << "Application Ranks(" << applications_[app].name << ") = ";
	for (std::vector<int>::const_iterator r=appRanks.begin(); r != appRanks.end(); ++r) {
	  if (r != appRanks.begin())
	    *pout_ << ",";
	  *pout_ << *r;
	}
	*pout_ << std::endl;
      }
      
      for (std::size_t t = 0; t < transferRanks_.size(); ++t) {
	*pout_ << "Transfer Ranks(" << transferNames_[t] << ") = ";
	for (std::vector<int>::const_iterator r=transferRanks_[t].begin(); r != transferRanks_[t].end(); ++r) {
	  if (r != transferRanks_[t].begin())
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
    // Application comms
    applicationComms_.resize(applications_.size());
    for (std::size_t app = 0; app != applications_.size(); ++app) {
      applicationComms_[app] = 
	globalComm_->createSubcommunicator(Teuchos::arrayViewFromVector(applications_[app].processes));
    }

    // Build transfer comms
    TEUCHOS_ASSERT(transfers_.size() == transferRanks_.size());
    TEUCHOS_ASSERT(transfers_.size() == transferNames_.size());
    TEUCHOS_ASSERT(transfers_.size() == transferNameToIndex_.size());
    for (std::size_t t = 0; t < transfers_.size(); ++t) {
      for (std::vector<ApplicationIndex>::const_iterator a = transfers_[t].begin(); a != transfers_[t].end(); ++a) {
	transferRanks_[t].insert(transferRanks_[t].end(), applications_[*a].processes.begin(), applications_[*a].processes.end());
      }
      std::sort(transferRanks_[t].begin(), transferRanks_[t].end());
      transferRanks_[t].erase( std::unique( transferRanks_[t].begin(), transferRanks_[t].end() ), transferRanks_[t].end() );
    }

    transferComms_.resize(transfers_.size());
    for (std::size_t t = 0; t < transferRanks_.size(); ++t)
      transferComms_[t] = globalComm_->createSubcommunicator(Teuchos::arrayViewFromVector(transferRanks_[t]));

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

    aout_.resize(applications_.size());
    for (std::size_t a = 0; a < applications_.size(); ++a) {
      aout_[a] = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      aout_[a]->pushLinePrefix(applications_[a].name);
      aout_[a]->setOutputToRootOnly(this->mapRankToCommWorldRank(applications_[a].startProcess,*globalComm_));
      aout_[a]->setShowLinePrefix(true);
      aout_[a]->setShowProcRank(true);
    }

    {
      pout_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
      for (std::size_t a = 0; a < applications_.size(); ++a) {
	if (this->appExistsOnProcess(a)) {
	  pout_->pushLinePrefix(applications_[a].name);
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
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (appIndex >= 0) && (appIndex < applications_.size()) );
#endif
    return applications_[appIndex].name;
  }
  
  std::string MultiphysicsDistributor::getTransferName(const TransferIndex transferIndex) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (transferIndex >= 0) && (transferIndex < transfers_.size()) );
#endif
    return transferNames_[transferIndex];
  }

  MultiphysicsDistributor::ApplicationIndex MultiphysicsDistributor::getApplicationIndex(const std::string& appName) const
  {
    std::map<std::string,ApplicationIndex>::const_iterator n = applicationNameToIndex_.find(appName);
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT(n != applicationNameToIndex_.end());
#endif
    return n->second;
  }
  
  MultiphysicsDistributor::TransferIndex MultiphysicsDistributor::getTransferIndex(const std::string& transferName) const
  {
    std::map<std::string,TransferIndex>::const_iterator n = transferNameToIndex_.find(transferName);
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT(n != transferNameToIndex_.end());
#endif
    return n->second;
  }

  bool MultiphysicsDistributor::appExistsOnProcess(const MultiphysicsDistributor::ApplicationIndex index) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (index >= 0) && (index < applications_.size()) );    
#endif
    return nonnull(applicationComms_[index]);
  }

  bool MultiphysicsDistributor::transferExistsOnProcess(const MultiphysicsDistributor::TransferIndex index) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (index >= 0) && (index < transfers_.size()) );
#endif
    return nonnull(transferComms_[index]);
  }

  Teuchos::RCP<const Teuchos::Comm<int> > MultiphysicsDistributor::getGlobalComm() const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT(nonnull(globalComm_));
#endif
    return globalComm_;
  }

  Teuchos::RCP<const Teuchos::Comm<int> >
  MultiphysicsDistributor::getAppComm(const ApplicationIndex index) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (index >= 0) && (index < applications_.size()) );    
#endif
    return applicationComms_[index];
  }

  Teuchos::RCP<const Teuchos::Comm<int> >
  MultiphysicsDistributor::getTransferComm(const TransferIndex index) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (index >= 0) && (index < transfers_.size()) );
#endif
    return transferComms_[index];
  }

  int MultiphysicsDistributor::getPrintRank(const ApplicationIndex index) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (index >= 0) && (index < applications_.size()) );    
#endif
    return applications_[index].startProcess;
  }

  Teuchos::RCP<Teuchos::FancyOStream> MultiphysicsDistributor::getSerialOStream() const
  { return out_; }

  Teuchos::RCP<Teuchos::FancyOStream> MultiphysicsDistributor::getApplicationOStream(const ApplicationIndex index) const
  {
#ifdef HAVE_PIKE_DEBUG
    TEUCHOS_ASSERT( (index >= 0) && (index < applications_.size()) );    
#endif
    return aout_[index];
  }

  Teuchos::RCP<Teuchos::FancyOStream> MultiphysicsDistributor::getParallelOStream() const
  { return pout_; }

}

#endif
