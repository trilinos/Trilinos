/*--------------------------------------------------------------------*/
/*    Copyright 2002 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef SIERRA_Ioss_FileInfo_h
#define SIERRA_Ioss_FileInfo_h

#include <string>
#include <iosfwd>

#include <string>
#include <sys/types.h>

namespace Ioss {

  /*! \class FileInfo 
   *  \author Greg Sjaardema
   *  \brief  Return information about the specified file.
   *
   *  A very minimal class (at least it used to be) for providing
   *  information about a file.  FileInfo provides information about a
   *  file's name, path, and type (directory, symbolic link, file).  Other
   *  information could be added as needed.  It currently does not cache
   *  any information, so if it is heavily used, a caching capability
   *  should be added.  See the Qt Toolkit QFileInfo class for a richer
   *  file class.
   */

  class FileInfo
  {
  public:
    //! Empty class referring to no file.
    FileInfo();
  
    //! Create object referring to file with name \a filename
    //! \param filename name of file
    explicit FileInfo(const std::string &filename);

    //! Create object referring to file with name \a filename
    //! \param filename name of file
    explicit FileInfo(const char   *filename);

    //! Copy constructor
    FileInfo(const FileInfo&);

    //! Constructor
    //! \param dirpath Directory Path
    //! \param filename base filename
    FileInfo(const std::string &dirpath, const std::string &filename);
  
    ~FileInfo();

    bool exists()        const; //!< returns True if file exists, false if nonexistant
    bool is_readable()   const; //!< Exists and is readable
    bool is_writable()   const; //!< Exists and is writable
    bool is_executable() const; //!< Exists and is executable

    bool is_file()       const; //!< Is a plain file
    bool is_dir()        const; //!< Is a directory
    bool is_symlink()    const; //!< Is a symbolic link to a file or directory
  
    time_t modified()    const; //!< Time of last data modification. See 'man stat(2)'
    time_t accessed()    const; //!< Time of last access
    time_t created()     const; //!< Time of last status change. (creation, chmod, ...)
  
    off_t  size()        const; //!< File size in bytes. Only if is_file() == true
  
    const std::string filename()  const; //!< Complete filename including path
    const std::string basename()  const; //!< strip path and extension
    const std::string tailname()  const; //!< basename() + extension()
    const std::string extension() const; //!< file extension.
    const std::string pathname()  const; //!< directory path, no filename

    void  set_filename(const std::string &name);
    void  set_filename(const char *name);
  
    bool  operator==(const FileInfo &other) const
    { return filename_ == other.filename_; }
  
    bool  operator!=(const FileInfo &other) const
    { return filename_ != other.filename_; }
  
    bool remove_file();
    
  private:
    std::string filename_;
    bool exists_;   //<! this is used frequently, check on creation
    bool readable_; //<! this is used frequently, check on creation
  };
}
#endif // SIERRA_Ioss_FileInfo_h

