#ifndef PIKE_MULTIPHYSICS_DISTRIBUTOR_HPP
#define PIKE_MULTIPHYSICS_DISTRIBUTOR_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_DefaultMpiComm.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <map>
#include <vector>
#include <string>

namespace pike {

  /** \brief Translates rank on comm A to its corresponding rank on comm B.
      
      @param rankA Rank of a process on Comm A
      @param commA MPI communicator A
      @param commB MPI communicator B
      @returns Rerturn the rank of rankA in the Comm B space

      \relates pike::MultiphysicsDistributor
   */
  int translateMpiRank(const int& rankA, 
		       const Teuchos::Comm<int>& commA,
		       const Teuchos::Comm<int>& commB);

  /** \brief Multiphysics driver utility that builds MPI
      sub-communicators and specialized ostreams for applications and
      data transfers.

      NOTE: All applications and transfers must be registered with the
      distributor on all processes of the global comm but do not have
      to actually exist on all processes of the global comm.
   */
  class MultiphysicsDistributor : public Teuchos::Describable,
				  public Teuchos::VerboseObject<pike::MultiphysicsDistributor> {

  public:

    typedef std::size_t ApplicationIndex;
    typedef std::size_t TransferIndex;

    MultiphysicsDistributor(const std::string& distributorName = "Pike");

    /** \brief Register a new application with this driver.
	
	This method allows applications to exist on the same set or a subset of MPI processes.

	\param[in] name Name of the application
	\param[in] processes A vector containing mpi processes that this application will be run on.  The process ranks are associated with the global comm used in the setup() method.
	\returns Index of the application
     */
    ApplicationIndex addApplication(const std::string& name, const std::vector<int>& processes);

    /** \brief Register a new application  with this driver given a range of ranks to exist on.
	
	This method allows applications to exist on the same set or a subset of MPI processes.

	\param[in] name Name of the application.
	\param[in] begin_rank The beginning of a range of processes that this application will exist on.  The range is inclusive of the end points, [begin_rank,end_rank].  The process ranks are associated with the global comm used in the setup() method.
	\param[in] end_rank The end of a range of processes that this application will exist on.  The range is inclusive of the end points, [begin_rank,end_rank].  The process ranks are associated with the global comm used in the setup() method.
	\returns Index of the application.
     */
    ApplicationIndex addApplication(const std::string& name, const int beginRank, const int endRank);

    /** \brief Tells this object that an active coupling between two physics exisits and that a union of the two applicaiton subcommunicators should be built for coupled data transfer.

       \param[in] a Index of the fist application involved in the transfer.
       \param[in] b Index of the second application involved in the transfer.
       \param[in] name Name of the transfer. Must be unique.
       \returns The transfer index.

       This is a simplification of the general addTranfer that takes a std::vector as its argument.  Most couplings are between two codes, and this case comes up so often that we have a specialized ctor for it.
     */
    TransferIndex addTransfer(const std::string& name, const ApplicationIndex a, const ApplicationIndex b);

    /** \brief Tells this object that an active coupling between multiple physics exisits and that a union of the applicaiton subcommunicators should be built for coupled data transfer.
 
       \param[in] appIndices Indices of the applications involved in the transfer.
       \param[in] name Name of the transfer.  Must be unique.
       \returns The transfer index.
    */
    TransferIndex addTransfer(const std::string& name, const std::vector<ApplicationIndex>& appIndices);

    /** \brief Tells this object that an active coupling between multiple physics exisits and that a subcommunicator should be built using the specified ranks.
 
       \param[in] name Name of the transfer.  Must be unique.
       \param[in] mpiRanks MPI Ranks corresponding to the global comm that are involved in the transfer.
       \returns The transfer index.
    */
    TransferIndex addTransferByRanks(const std::string& name, const std::vector<int>& mpiRanks);

    /** \brief Builds the application subcommunicators and any coupling subcommunicators. 
  
       \param[in] globalComm The global comm across which all applications are contained.  Apps do not have to exist on every process of the global comm, but all application subcomms must be contained in this comm.

       NOTE: The globalComm supplied by the method call is duplicated() in the setup to get a unique context space for MultiphysicsDistributor opperations.  the call to getGlobalComm returns the internal duplicated comm that is separate form the one supplied by the user.   
     */
    void setup(const Teuchos::RCP<const Teuchos::Comm<int> >& globalComm, bool printCommDistribution = false);

    /** \brief Returns the application name given the application index. */
    std::string getApplicationName(const ApplicationIndex appIndex) const;

    /** \brief Returns the transfer name given the transfer index. */
    std::string getTransferName(const TransferIndex transferIndex) const;

    /** \brief Returns the application index given the application string identifier. */
    ApplicationIndex getApplicationIndex(const std::string& appName) const;

    /** \brief Returns the transfer index given the transfer string identifier. */
    TransferIndex getTransferIndex(const std::string& transferName) const;

    /** \brief Returns true if the application is active on this MPI process. */
    bool appExistsOnProcess(const ApplicationIndex index) const;

    /** \brief Returns true if the application is active on this MPI process. */
    bool transferExistsOnProcess(const TransferIndex index) const;

    /** \brief Returns the global comm. */
    Teuchos::RCP<const Teuchos::Comm<int> > getGlobalComm() const;

    /** \brief Returns the subcommunicator for the application index.  Returns a null RCP if app is not active on this process. */ 
    Teuchos::RCP<const Teuchos::Comm<int> > getAppComm(const ApplicationIndex applicationIndex) const;

    /** \brief Returns the subcommunicator for the transfer index.  Returns a null RCP if neither of the coupled applications is active on this process.*/ 
    Teuchos::RCP<const Teuchos::Comm<int> > getTransferComm(const TransferIndex transferIndex) const;

    /** \brief Returns the rank of the print process for the given application index.  This rank corresponds to the global communicator. */
    int getPrintRank(const ApplicationIndex index) const;

    //! Prints to only the print process of the global communicator.
    Teuchos::RCP<Teuchos::FancyOStream> getSerialOStream() const;
    
    //! Prints to a single process of the application communicator.
    Teuchos::RCP<Teuchos::FancyOStream> getApplicationOStream(const ApplicationIndex applicationIndex) const;

    //! Prints to all processes
    Teuchos::RCP<Teuchos::FancyOStream> getParallelOStream() const;

  private:

    struct ApplicationData {

      std::string name;
      std::vector<int> processes;
      int numProcesses;
      int startProcess;

      ApplicationData() : name(""),processes(0),numProcesses(-1),startProcess(-1)
      { }

    };

    //! Builds all communicators.
    void buildComms();

    //! Builds all ostreams.
    void buildOStreams();

    /** Teuchos::FancyOStream uses mpi comm world to determine if the
	current processt is the print process.  So when settting print
	processes you need to map it to the world rank.  This function
	does that mapping.
     */ 
    int mapRankToCommWorldRank(int rank, const Teuchos::Comm<int>& comm);

  private:
    std::string myName_;

    Teuchos::RCP<const Teuchos::Comm<int> > globalComm_;

    bool setupCalled_;

    std::vector<ApplicationData> applications_;

    std::map<std::string,ApplicationIndex> applicationNameToIndex_;

    //! A nonnull RCP means that this application is instantiatied on this process.
    std::vector<Teuchos::RCP<const Teuchos::Comm<int> > > applicationComms_;

    //! Each transfer contains a vector of application indices associated with the transfer.
    std::vector<std::vector<ApplicationIndex> > transfers_;

    //! A nonnull RCP means that this transfer is instantiatied on this process.  A transfer comm must represent at least the union of all application comms.
    std::vector<Teuchos::RCP<const Teuchos::Comm<int> > > transferComms_;

    std::vector<std::vector<int> > transferRanks_;

    std::vector<std::string> transferNames_;

    std::map<std::string,TransferIndex> transferNameToIndex_;

    //! Serial ostream that prints to global process 0.
    Teuchos::RCP<Teuchos::FancyOStream> out_;
    //! Application ostream that prints on the first process of each application comm.
    std::vector<Teuchos::RCP<Teuchos::FancyOStream> > aout_;
    //! Transfer opstream that prints on the first process of each tranfer comm.
    std::vector<Teuchos::RCP<Teuchos::FancyOStream> > tout_;
    //! Parallel ostream that prints on each process.
    Teuchos::RCP<Teuchos::FancyOStream> pout_;
  };

}

#endif
