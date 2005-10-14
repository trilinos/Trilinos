// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ***********************************************************************
// @HEADER

#ifndef TEUCHOS_FANCY_O_STREAM_HPP
#define TEUCHOS_FANCY_O_STREAM_HPP

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_GlobalMPISession.hpp"

namespace Teuchos {


/** \brief Stream buffering class that performs the magic of indenting
 * data sent to an std::ostream object.
 *
 * Note, this is not a user-level class.  Users should use
 * <tt>basic_FancyOStream</tt>.
 */
template<typename CharT, typename Traits>
class basic_FancyOStream_buf : public std::basic_streambuf<CharT,Traits>
{
public:
  
  /** \brief . */
  typedef CharT                             char_type;
  /** \brief . */
  typedef Traits  					                traits_type;
  /** \brief . */
  typedef typename traits_type::int_type 		int_type;
  /** \brief . */
  typedef typename traits_type::pos_type 		pos_type;
  /** \brief . */
  typedef typename traits_type::off_type 		off_type;

  /** \brief . */
  basic_FancyOStream_buf(
    const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr
    ,const int                                                     startingTab
    ,const bool                                                    showLinePrefix
    ,const int                                                     maxLenLinePrefix
    ,const bool                                                    showTabCount
    ,const bool                                                    showProcRank
    );

  /** \brief . */
  void initialize(
    const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr
    ,const int                                                     startingTab
    ,const bool                                                    showLinePrefix
    ,const int                                                     maxLenLinePrefix
    ,const bool                                                    showTabCount
    ,const bool                                                    showProcRank
    );

  /** \brief . */
  RefCountPtr<std::basic_ostream<char_type,traits_type> > getOStream();

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTabIndentStr();

  /** \brief .*/
  void setShowLinePrefix(const bool showLinePrefix);

  /** \brief .*/
  void setMaxLenLinePrefix(const bool maxLenLinePrefix);

  /** \brief . */
  void setShowTabCount(const bool showTabCount);

  /** \brief . */
  void setShowProcRank(const bool showProcRank);

  /** \brief .*/
  void setProcRankAndSize( const int procRank, const int numProcs );

  /** \brief . */
  void pushTab(const int tabs);

  /** \brief . */
  void popTab();

  /** \brief . */
  void pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix);

  /** \brief . */
  void popLinePrefix();

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTopLinePrefix() const;

  /** \brief . */
  void pushDisableTabbing();

  /** \brief . */
  void popDisableTabbing();
  
protected:
  
  /** \name Protected overridden functions from std::basic_streambuf */
  //@{

  /** \brief . */
  std::streamsize xsputn(const char_type* s, std::streamsize n);

  /** \brief . */
  int_type overflow(int_type c);

#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS

  void imbue(const locale& l) 
    {
      std::cerr << "\ncalled imbue()\n";
      std::basic_streambuf<CharT,Traits>::imbue(l);
    }
  
  basic_streambuf<char_type,Traits>* 
  setbuf(char_type* s, streamsize n)
    {
      std::cerr << "\ncalled setbuf()\n";
      return std::basic_streambuf<CharT,Traits>::setbuf(s,n);
    }
      
  pos_type 
  seekoff(off_type a, ios_base::seekdir b,ios_base::openmode c)
    {
      std::cerr << "\ncalled seekoff()\n";
      return std::basic_streambuf<CharT,Traits>::seekoff(a,b,c);
    }

  pos_type 
  seekpos(pos_type a, ios_base::openmode b)
    {
      std::cerr << "\ncalled seekpos()\n";
      return std::basic_streambuf<CharT,Traits>::seekpos(a,b);
    }
  
  int 
  sync()
    {
      std::cerr << "\ncalled sync()\n";
      return std::basic_streambuf<CharT,Traits>::sync();
    }
  
  streamsize 
  showmanyc()
    {
      std::cerr << "\ncalled showmanyc()\n";
      return std::basic_streambuf<CharT,Traits>::showmanyc();
    }
  
  streamsize 
  xsgetn(char_type* s, streamsize n)
    {
      std::cerr << "\ncalled xsgetn()\n";
      return std::basic_streambuf<CharT,Traits>::xsgetn(s,n);
    }
  
  int_type 
  underflow()
    {
      std::cerr << "\ncalled underflow()\n";
      return std::basic_streambuf<CharT,Traits>::underflow();
    }

  int_type 
  uflow() 
    {
      std::cerr << "\ncalled uflow()\n";
      return std::basic_streambuf<CharT,Traits>::uflow();
    }

  int_type 
  pbackfail(int_type c = traits_type::eof())
    {
      std::cerr << "\ncalled pbackfail()\n";
      return std::basic_streambuf<CharT,Traits>::pbackfail(c);
    }

#endif // TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS

  //@}

private:

  // ////////////////////////
  // Private types

  typedef std::basic_string<char_type,traits_type>   string_t;
  typedef std::deque<int>                            tabIndentStack_t;
  typedef std::deque<string_t>                       linePrefixStack_t;

  // ////////////////////////
  // Private data members

  RefCountPtr<std::basic_ostream<char_type,traits_type> >  oStream_;
  std::basic_string<char_type,traits_type>                 tabIndentStr_;
  bool                                                     showLinePrefix_;
  int                                                      maxLenLinePrefix_;
  bool                                                     showTabCount_;
  bool                                                     showProcRank_;
  int                                                      procRank_;
  int                                                      numProcs_;
  int                                                      rankPrintWidth_;
  
  int                      tabIndent_;
  tabIndentStack_t         tabIndentStack_;
  linePrefixStack_t        linePrefixStack_;
  int                      enableTabbingStack_;
  
  bool                     wroteNewline_;

  // ////////////////////////
  // Private member functions

  void writeChars( const char_type s[], std::streamsize n );

  void writeFrontMatter();
  
};

/** \brief std::ostream subclass that performs the magic of indenting data
 * sent to an std::ostream object among other things.
 *
 * Use the typedef <tt>FancyOStream</tt> for support for the <tt>char</tt>
 * character type.
 *
 * Indentation of the stream is accomplished through creating
 * <tt>basic_OSTab</tt> objects.
 *
 * In addition to indenting output, this stream object can also print various
 * types of information at the beginning of each line.  The type of information
 * supported is:
 * <ul>
 * <li> Processor rank: Set using <tt>setShowProcRank()</tt>.
 * <li> Line prefix name: Set using <tt>showLinePrefix()</tt> and <tt>OSTab::OSTab()</tt>.
 * <li> Tab counts (useful for debugging): Set using <tt>setShowTabCount()</tt>.
 * </ul>
 *
 * See <tt>FancyOutputting_test.cpp</tt> for examples of how this class is
 * used and the output it generates.
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_FancyOStream : public basic_ostream<CharT, Traits>
{
public:

  /** \name Public types */
  //@{

  /** \brief . */
  typedef CharT 					                  char_type;
  /** \brief . */
  typedef Traits 					                  traits_type;
  /** \brief . */
  typedef typename traits_type::int_type 		int_type;
  /** \brief . */
  typedef typename traits_type::pos_type 		pos_type;
  /** \brief . */
  typedef typename traits_type::off_type 		off_type;
  /** \brief . */

  /** \brief . */
  typedef basic_FancyOStream_buf<CharT,Traits> 	streambuf_t;
  /** \brief . */
  typedef basic_ostream<char_type, traits_type>	ostream_t;

  //@}

  /** \name Public client functions */
  //@{

  /** \brief . */
  explicit
  basic_FancyOStream(
    const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr     = "\t"
    ,const int                                                     startingTab       = 0
    ,const bool                                                    showLinePrefix    = false
    ,const int                                                     maxLenLinePrefix  = 10
    ,const bool                                                    showTabCount      = false
    ,const bool                                                    showProcRank      = false
    );

  /** \brief . */
  void initialize(
    const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr     = "\t"
    ,const int                                                     startingTab       = 0
    ,const bool                                                    showLinePrefix    = false
    ,const int                                                     maxLenLinePrefix  = 10
    ,const bool                                                    showTabCount      = false
    ,const bool                                                    showProcRank      = false
    );

  /** \brief. */
  RefCountPtr<std::basic_ostream<char_type,traits_type> > getOStream();

  /** \brief. */
  const std::basic_string<char_type,traits_type>& getTabIndentStr();

  /** \brief Set if processor rank, line prefixes, and tab counts are shown or not .*/
  basic_FancyOStream& setShowAllFrontMatter(const bool showAllFrontMatter);

  /** \brief .*/
  basic_FancyOStream& setShowLinePrefix(const bool showLinePrefix);

  /** \brief .*/
  basic_FancyOStream& setMaxLenLinePrefix(const bool maxLenLinePrefix);

  /** \brief . */
  basic_FancyOStream& setShowTabCount(const bool showTabCount);

  /** \brief . */
  basic_FancyOStream& setShowProcRank(const bool showProcRank);

  /** \brief . */
  basic_FancyOStream& setProcRankAndSize( const int procRank, const int numProcs );

  //@}

  /** \name Functions designed to be used by basic_OSTab */
  //@{

  /** \brief. */
  void pushTab(const int tabs = 1);

  /** \brief. */
  void popTab();

  /** \brief . */
  void pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix);

  /** \brief . */
  void popLinePrefix();

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTopLinePrefix() const;

  /** \brief . */
  void pushDisableTabbing();

  /** \brief . */
  void popDisableTabbing();

  //@}
  
private:

  streambuf_t	streambuf_;

};

/** \brief Tabbing class for helping to create formated, indented output for a
 * <tt>basic_FancyOStream</tt> object.
 *
 * Use the typedef <tt>OSStream</tt> for support for the <tt>char</tt>
 * character type.
 *
 * This class is used to create tab indents and set line prefix names for
 * output that is generated by a <tt>basic_FancyOStream</tt> object.
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_OSTab
{
public:

  /** \brief . */
  static const int DISABLE_TABBING = -99999; // This magic number should be just fine!
  /** \brief. */
  basic_OSTab(
    const RefCountPtr<basic_FancyOStream<CharT,Traits> >   &fancyOStream
    ,const int                                             tabs = 1
    ,const std::basic_string<CharT,Traits>                 linePrefix = ""
    )
    :fancyOStream_(fancyOStream)
    ,tabs_(tabs)
    ,linePrefix_(linePrefix)
    {
      updateState();
    }
  /** \brief. */
  basic_OSTab( const basic_OSTab &osTab )
    :fancyOStream_(osTab.fancyOStream_)
    ,tabs_(osTab.tabs_)
    {
      updateState();
    }
  /** \brief. */
  ~basic_OSTab()
    {
      if(tabs_==DISABLE_TABBING)
        fancyOStream_->popDisableTabbing();
      else
        fancyOStream_->popTab();
      if(linePrefix_.length()) fancyOStream_->popLinePrefix();
    }
  /** \brief. */
  basic_OSTab<CharT,Traits>& operator=( const basic_OSTab &osTab )
    {
      fancyOStream_ = osTab.fancyOStream_;
      tabs_ = osTab.tabs_;
      updateState();
      return *this;
    }
  
private:
  
  RefCountPtr<basic_FancyOStream<CharT,Traits> >  fancyOStream_;
  int                                             tabs_;
  std::basic_string<CharT,Traits>                 linePrefix_;

  void updateState()
    {
      if(tabs_==DISABLE_TABBING)
        fancyOStream_->pushDisableTabbing();
      else
        fancyOStream_->pushTab(tabs_);
      if(linePrefix_.length()) fancyOStream_->pushLinePrefix(linePrefix_);
    }
  
};

// ///////////////////////////////
// Typedefs

/** \brief . */
typedef basic_FancyOStream<char> FancyOStream;

/** \brief . */
typedef basic_OSTab<char> OSTab;

/** \brief .*/
#define TEUCHOS_OSTAB ::Teuchos::OSTab __localThisTab = this->getOSTab()

// ////////////////////////////////
// Defintions

//
// basic_FancyOStream_buf
//

template<typename CharT, typename Traits>
basic_FancyOStream_buf<CharT,Traits>::basic_FancyOStream_buf(
  const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
  )
{
  this->initialize(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::initialize(
  const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
  )
{
  oStream_ = oStream;
  tabIndentStr_ = tabIndentStr;
  showLinePrefix_ = showLinePrefix;
  maxLenLinePrefix_ = maxLenLinePrefix;
  showTabCount_ = showTabCount;
  showProcRank_ = showProcRank;
  procRank_ = GlobalMPISession::getRank();
  numProcs_ = GlobalMPISession::getNProc();
  rankPrintWidth_ = int(std::log10(float(numProcs_)))+1;
  tabIndent_ = startingTab;
  tabIndentStack_.resize(0);
  linePrefixStack_.resize(0);
  wroteNewline_ = true;
  enableTabbingStack_ = 0;
}

template<typename CharT, typename Traits>
RefCountPtr<std::basic_ostream<CharT,Traits> >
basic_FancyOStream_buf<CharT,Traits>::getOStream()
{
  return oStream_;
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream_buf<CharT,Traits>::getTabIndentStr()
{
  return tabIndentStr_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowLinePrefix(const bool showLinePrefix)
{
  showLinePrefix_ = showLinePrefix;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setMaxLenLinePrefix(const bool maxLenLinePrefix)
{
  TEST_FOR_EXCEPT( maxLenLinePrefix >= 5 );
  maxLenLinePrefix_ = maxLenLinePrefix;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowTabCount(const bool showTabCount)
{
  showTabCount_ = showTabCount;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowProcRank(const bool showProcRank)
{
  showProcRank_ = showProcRank;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setProcRankAndSize( const int procRank, const int numProcs )
{
  procRank_ = procRank;
  numProcs_ = numProcs;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushTab(const int tabs)
{
  tabIndentStack_.push_back(tabs);
  tabIndent_ += tabs;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popTab()
{
  tabIndent_ -= tabIndentStack_.back();
  tabIndentStack_.pop_back();
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix)
{
  TEST_FOR_EXCEPT(static_cast<int>(linePrefix.length()) > maxLenLinePrefix_);
  linePrefixStack_.push_back(linePrefix);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popLinePrefix()
{
  linePrefixStack_.pop_back();
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream_buf<CharT,Traits>::getTopLinePrefix() const
{
  return linePrefixStack_.back();
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushDisableTabbing()
{
  ++enableTabbingStack_;
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popDisableTabbing()
{
  --enableTabbingStack_;
}

// protected

template<typename CharT, typename Traits>
std::streamsize basic_FancyOStream_buf<CharT,Traits>::xsputn(const char_type* s, std::streamsize n)
{
#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS
  std::cerr << "\ncalled xsputn()\n";
#endif
  writeChars(s,n);
  return n;
}

template<typename CharT, typename Traits>
typename basic_FancyOStream_buf<CharT,Traits>::int_type 
basic_FancyOStream_buf<CharT,Traits>::overflow(int_type c)
{
#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS
  std::cerr << "\ncalled overflow()\n";
#endif
  if(c != traits_type::eof()) {
    const char_type cc[] = { traits_type::to_char_type(c) };
    this->writeChars(cc,1);
  }
  return traits_type::not_eof(c);
  //return std::basic_streambuf<CharT,Traits>::overflow(c);
}

// private

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::writeChars( const char_type s[], std::streamsize n )
{
  if(n == 0) return;
  std::streamsize p = 0, first_p = 0;
  bool done_outputting = false;
  const char_type newline = '\n';
  while( !done_outputting ) {
    // Find the next newline
    for( p = first_p; p < n; ++p ) {
      if(s[p] == newline) {
        break;
      }
    }
    if(p == n) {
      // We did not find a newline at the end!
      --p;
      done_outputting = true;
    }
    else if( p == n-1 && s[p] == newline ) {
      // The last character in the string is a newline
      done_outputting = true;
    }
    // Write the beginning of the line if we need to
    if(wroteNewline_) {
      writeFrontMatter();
      wroteNewline_ = false;
    }
    // Write up to the newline or the end of the string
    oStream_->write(s+first_p,p-first_p+1);
    if(s[p] == newline) {
      wroteNewline_ = true;
    }
    // Update for next search
    if(!done_outputting)
      first_p = p+1; 
  }
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::writeFrontMatter()
{
  bool didOutput = false;
  if(showProcRank_) {
    *oStream_ << "p=" << std::right << std::setw(rankPrintWidth_) << procRank_;
    didOutput = true;
  }
  if(showLinePrefix_) {
    if(didOutput)
      *oStream_ << ", ";
    *oStream_ << std::left << std::setw(maxLenLinePrefix_);
    if(linePrefixStack_.size()) {
      *oStream_ << this->getTopLinePrefix();
    }
    else {
      *oStream_ << "";
    }
    didOutput = true;
  }
  if(showTabCount_) {
    if(didOutput)
      *oStream_ << ", ";
    *oStream_ << "tabs=" << std::right << std::setw(2) << tabIndent_;
    didOutput = true;
  }
  // ToDo: Add the Prefix name if asked
  // ToDo: Add the processor number if asked
  // ToDo: Add the number of indents if asked
  if(didOutput) {
    *oStream_ << " |" << tabIndentStr_;
  }
  if(enableTabbingStack_==0) {
    for( int i = 0; i < tabIndent_; ++i )
      *oStream_ << tabIndentStr_;
  }
}

//
// basic_FancyOStream
//

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>::basic_FancyOStream(
  const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
 )
  : ostream_t(NULL), streambuf_(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank)
{
  this->init(&streambuf_);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::initialize(
  const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  ,const bool                                                    showProcRank
  )
{
  streambuf_.initialize(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank);
  this->init(&streambuf_);
}

template<typename CharT, typename Traits>
RefCountPtr<std::basic_ostream<CharT,Traits> >
basic_FancyOStream<CharT,Traits>::getOStream()
{
  return streambuf_.getOStream();
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::getTabIndentStr()
{
  return streambuf_.getTabIndentStr();
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowAllFrontMatter(const bool showAllFrontMatter)
{
  streambuf_.setShowLinePrefix(showAllFrontMatter);
  streambuf_.setShowTabCount(showAllFrontMatter);
  streambuf_.setShowProcRank(showAllFrontMatter);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowLinePrefix(const bool showLinePrefix)
{
  streambuf_.setShowLinePrefix(showLinePrefix);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
 basic_FancyOStream<CharT,Traits>::setMaxLenLinePrefix(const bool maxLenLinePrefix)
{
  streambuf_.setMaxLenLinePrefix(maxLenLinePrefix);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowTabCount(const bool showTabCount)
{
  streambuf_.setShowTabCount(showTabCount);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowProcRank(const bool showProcRank)
{
  streambuf_.setShowProcRank(showProcRank);
  return *this;
}

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setProcRankAndSize( const int procRank, const int numProcs )
{
  streambuf_.setProcRankAndSize(procRank,numProcs);
  return *this;
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushTab(const int tabs)
{
  streambuf_.pushTab(tabs);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popTab()
{
  streambuf_.popTab();
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix)
{
  streambuf_.pushLinePrefix(linePrefix);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popLinePrefix()
{
  streambuf_.popLinePrefix();
}

template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::getTopLinePrefix() const
{
  return streambuf_.getTopLinePrefix();
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushDisableTabbing()
{
  streambuf_.pushDisableTabbing();
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popDisableTabbing()
{
  return streambuf_.popDisableTabbing();
}

//
// OSTab
//

} // namespace Teuchos

#endif // TEUCHOS_FANCY_O_STREAM_HPP
