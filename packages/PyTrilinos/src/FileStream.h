#include <streambuf>

class FileStream : public std::streambuf {

public:
  FileStream(std::FILE* file);

protected:
  std::streambuf::int_type overflow(std::streambuf::int_type c);
  FILE* self_file;

};
