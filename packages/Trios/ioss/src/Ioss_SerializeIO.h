#ifndef SIERRA_Ioss_SerializeIO_h
#define SIERRA_Ioss_SerializeIO_h

#include <Ioss_CodeTypes.h>
#include <string>

namespace Ioss {

class DatabaseIO;

/**
 * @brief Class <code>SerializeIO</code> is a sentry class which performs serialization
 * for mesh database I/O.
 *
 * This sentry guards serialization of parallel I/O routines.  At construction, it
 * blocks the processes via an MPI barrier, releasing them to execute in groups specified
 * by <code>s_groupSize</code>.  At destruction, it continues to block via MPI barriers
 * until all the processor have been released by the constructor.
 *
 * In the case where the constructor is called, and the sentry is already active and owned
 * by the processes group, the constructor and destrutor simply fall through since the
 * serialization is already in place at a higher level.
 *
 * The MPI
 *
 */
class SerializeIO
{
public:
  /**
   * Creates a new <code>SerializeIO</code> instance.
   *
   * @param database_io	a <code>DatabaseIO</code> variable ...
   */
  explicit SerializeIO(const DatabaseIO *database_io, int manual_owner_processor = -1);
  ~SerializeIO();

  inline static int getOwner() {
    return s_owner;
  }

  inline static int getRank() {
    return s_rank;
  }

  inline static int getSize() {
    return s_size;
  }

  inline static int getGroupRank() {
    return s_groupRank;
  }

  inline static int getGroupSize() {
    return s_groupSize;
  }

  static void setGroupFactor(int factor);

  inline static bool isEnabled() {
    return s_groupFactor != 0;
  }

  inline static bool inBarrier() {
    return s_owner != -1;
  }

  inline static bool inMyGroup() {
    return s_owner == s_groupRank;
  }

private:
  SerializeIO(const SerializeIO& from); // do not implement
  SerializeIO& operator=(const SerializeIO& from); // do not implement

  const DatabaseIO *	m_databaseIO;		///< Database I/O pointer
  const bool		m_activeFallThru;	///< No barries since my group is running
  int			m_manualOwner;		///< Manually specified owner

  static int		s_groupFactor;		///< Grouping factor
  static int		s_size;			///< Number of processors
  static int		s_rank;			///< My processor rank
  static int		s_groupSize;		///< Number of groups
  static int		s_groupRank;		///< My group rank
  static int		s_owner;		///< Group currently running
};

} // namespace Ioss

#endif // SIERRA_Ioss_SerializeIO_h
