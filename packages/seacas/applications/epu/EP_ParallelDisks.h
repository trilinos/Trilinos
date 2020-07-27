/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#ifndef ParallelDisksH
#define ParallelDisksH

#include <string>

/*****************************************************************************/
namespace Excn {
  class ParallelDisks
  {
    //: This class declares the unique, global object that represents the
    //: disk farm on a parallel machine.  Its function is to map
    //: processors onto disk files.

  public:
    ParallelDisks();
    ~ParallelDisks();

    void Number_of_Raids(int /*i*/);
    void Raid_Offset(int /*i*/);

    int Number_of_Raids() const;
    int Raid_Offset() const;

    static void Create_IO_Filename(std::string & /*name*/, int processor, int num_processors);

    void rename_file_for_mp(const std::string &rootdir, const std::string &subdir,
                            std::string &name, int node, int numproc) const;

  private:
    std::vector<std::string> disk_names{};

    int number_of_raids{0};
    int raid_offset{0};

    void create_disk_names();

    // not defined (the parallel disks object is unique!)
    ParallelDisks(const ParallelDisks &);
    ParallelDisks &operator=(const ParallelDisks &);
  };
} // namespace Excn
#endif
