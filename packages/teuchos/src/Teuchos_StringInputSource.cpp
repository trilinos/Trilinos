#include "Teuchos_StringInputSource.hpp"
#include "Teuchos_StringInputStream.hpp"

using namespace Teuchos;


StringInputSource::StringInputSource(const string& text)
	: XMLInputSource(), text_(text)
{;}

RefCountPtr<XMLInputStream> StringInputSource::stream() const 
{
	return rcp(new StringInputStream(text_));
}

