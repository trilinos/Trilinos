#include <streambuf>

class FileStream : public std::streambuf {

  public:

  FileStream(std::FILE* file): self_file(file) {}

  protected:

    std::streambuf::int_type overflow(std::streambuf::int_type c) {
      return std::fputc(c, self_file) == EOF?
	std::streambuf::traits_type::eof(): c;
    }

    FILE* self_file;
};
