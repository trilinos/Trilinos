#include "Teuchos_FileInputSource.hpp"
#include "Teuchos_FileInputStream.hpp"

using namespace Teuchos;


FileInputSource::FileInputSource(const string& filename)
	: XMLInputSource(), filename_(filename)
{;}

RefCountPtr<XMLInputStream> FileInputSource::stream() const 
{
	return rcp(new FileInputStream(filename_));
}

