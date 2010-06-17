/*--------------------------------------------------------------------*/
/*    Copyright 2000 - 2009 Sandia Corporation.                       */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef IOSS_Ioss_IOFactory_h
#define IOSS_Ioss_IOFactory_h

#include <Ioss_CodeTypes.h>
#include <string>

#include <Ioss_DBUsage.h>

#include <map>
#include <vector>

namespace Ioss {

class IOFactory;

typedef std::vector<std::string> NameList;
typedef std::map<std::string, IOFactory*, std::less<std::string> > IOFactoryMap;
typedef IOFactoryMap::value_type IOFactoryValuePair;

class DatabaseIO;
class IOFactory {
public:
  virtual ~IOFactory() {};
  static DatabaseIO* create(const std::string& type,
                            const std::string& filename,
                            DatabaseUsage db_usage,
                            MPI_Comm communicator = MPI_COMM_WORLD);

  static int describe(NameList *names);
  static void clean();
    
protected:
  explicit IOFactory(const std::string& type);

  virtual DatabaseIO* make_IO(const std::string& filename,
                              DatabaseUsage db_usage,
                              MPI_Comm communicator) const = 0;

  static void alias(const std::string& base, const std::string& syn);

private:
  //    static IOFactoryMap* registry_;
  static IOFactoryMap* registry();
};

}
#endif
