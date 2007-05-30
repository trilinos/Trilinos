

#include "Teuchos_FileInputSource.hpp"

// Boost initialization
#include <boost/python.hpp>
using namespace boost::python;
#include "teuchos_call_policies.hpp"


void expose_fileinputsource()
{
      // FileInputSource(const string& filename);
    class_<Teuchos::FileInputSource>("FileInputSource",
            "Definition of XMLInputSource derived class for reading XML from a file",
            init<string>() )
        // virtual RefCountPtr<XMLInputStream> stream() const;
        .def( "stream", &Teuchos::FileInputSource::stream ,
                        del_dangleing_rcp<>(),
                        "Create a FileInputStream")
    ;
}
