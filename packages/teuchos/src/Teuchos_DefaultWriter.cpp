#include "Teuchos_DefaultWriter.hpp"
#include "Teuchos_Utils.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_StrUtils.hpp"

using namespace Teuchos;

string& DefaultWriter::header()
{
	static string rtn = "<p=";
	return rtn;
}

DefaultWriter::DefaultWriter()
	: os_(std::cerr)
{}

DefaultWriter::DefaultWriter(std::ostream& os)
	: os_(os)
{}

DefaultWriter::DefaultWriter(const RefCountPtr<ostream>& os_ptr)
	: os_ptr_(os_ptr), os_(*os_ptr_)
{}

void DefaultWriter::print(const std::string& msg)
{
	TeuchosOStringStream ss;
	string head = header() + Utils::toString(MPISession::getRank()) + "> ";
	unsigned int maxLineSize = 78 - head.length();

	Array<string> tokens = StrUtils::getTokensPlusWhitespace(msg);

	unsigned int lineSize = 0;
	ss << head;
	for (int i=0; i<tokens.length(); i++)
		{
			if (lineSize+tokens[i].length() > maxLineSize)
				{
					if (StrUtils::isWhite(tokens[i])) continue;
					ss << std::endl << head << tokens[i];
					lineSize = 0;
				}
			else
				{
					ss << tokens[i];
				}
			lineSize += tokens[i].length();
		}
#if !HAVE_SSTREAM
	ss << std::ends;
#endif
	os_ << ss.str();
}

void DefaultWriter::println(const std::string& msg)
{
	print(msg);
	os_ << std::endl;
}
