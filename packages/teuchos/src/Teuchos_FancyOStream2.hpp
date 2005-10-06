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
    );

  /** \brief . */
  void initialize(
    const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr
    ,const int                                                     startingTab
    ,const bool                                                    showLinePrefix
    ,const int                                                     maxLenLinePrefix
    ,const bool                                                    showTabCount
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
  void pushTab(const int tabs);

  /** \brief . */
  void popTab();

  /** \brief . */
  void pushLinePrefix(const std::basic_string<char_type,traits_type> &linePrefix);

  /** \brief . */
  void popLinePrefix();

  /** \brief . */
  const std::basic_string<char_type,traits_type>& getTopLinePrefix() const;
  
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

  typedef std::basic_string<char_type,traits_type> string_t;

  typedef std::deque<int>       tabIndentStack_t;
  typedef std::deque<string_t>  linePrefixStack_t;

  // ////////////////////////
  // Private data members

  RefCountPtr<std::basic_ostream<char_type,traits_type> >  oStream_;
  std::basic_string<char_type,traits_type>                 tabIndentStr_;
  bool                                                     showLinePrefix_;
  int                                                      maxLenLinePrefix_;
  bool                                                     showTabCount_;
  
  int                      tabIndent_;
  tabIndentStack_t         tabIndentStack_;
  linePrefixStack_t        linePrefixStack_;
  
  bool                     wroteNewline_;

  // ////////////////////////
  // Private member functions

  void writeChars( const char_type s[], std::streamsize n );

  void writeBeginningOfLine();
  
};

/** \brief std::ostream subclass that performs the magic of indenting data
 * sent to an std::ostream object.
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_FancyOStream : public basic_ostream<CharT, Traits>
{
public:

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

  /** \brief . */
  explicit
  basic_FancyOStream(
    const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr     = "\t"
    ,const int                                                     startingTab       = 0
    ,const bool                                                    showLinePrefix    = false
    ,const int                                                     maxLenLinePrefix  = 10
    ,const bool                                                    showTabCount      = false
    );

  /** \brief . */
  void initialize(
    const RefCountPtr< std::basic_ostream<char_type,traits_type> > &oStream
    ,const std::basic_string<char_type,traits_type>                &tabIndentStr     = "\t"
    ,const int                                                     startingTab       = 0
    ,const bool                                                    showLinePrefix    = false
    ,const int                                                     maxLenLinePrefix  = 10
    ,const bool                                                    showTabCount      = false
    );

  /** \brief. */
  RefCountPtr<std::basic_ostream<char_type,traits_type> > getOStream();

  /** \brief. */
  const std::basic_string<char_type,traits_type>& getTabIndentStr();

  /** \brief .*/
  void setShowLinePrefix(const bool showLinePrefix);

  /** \brief .*/
  void setMaxLenLinePrefix(const bool maxLenLinePrefix);

  /** \brief . */
  void setShowTabCount(const bool showTabCount);

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
  
private:

  streambuf_t	streambuf_;

};

/** \brief Tabbing class for helping to create formated output.
 */
template <typename CharT, typename Traits = std::char_traits<CharT> >
class basic_OSTab
{
public:
  
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
      fancyOStream_->pushTab(tabs_);
      if(linePrefix_.length()) fancyOStream_->pushLinePrefix(linePrefix_);
    }
  
  /** \brief. */
  basic_OSTab( const basic_OSTab &osTab )
    :fancyOStream_(osTab.fancyOStream_)
    ,tabs_(osTab.tabs_)
    {
      fancyOStream_->pushTab(tabs_);
      if(linePrefix_.length()) fancyOStream_->pushLinePrefix(linePrefix_);
    }
  
  /** \brief. */
  ~basic_OSTab()
    {
      fancyOStream_->popTab();
      if(linePrefix_.length()) fancyOStream_->popLinePrefix();
    }
  
  /** \brief. */
  basic_OSTab<CharT,Traits>& operator=( const basic_OSTab &osTab )
    {
      fancyOStream_ = osTab.fancyOStream_;
      tabs_ = osTab.tabs_;
      fancyOStream_->pushTab(tabs_);
      if(linePrefix_.length()) fancyOStream_->pushLinePrefix(linePrefix_);
      return *this;
    }
  
private:
  
  RefCountPtr<basic_FancyOStream<CharT,Traits> >  fancyOStream_;
  int                                             tabs_;
  std::basic_string<CharT,Traits>                 linePrefix_;

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
  )
{
  this->initialize(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::initialize(
  const RefCountPtr<std::basic_ostream<char_type,traits_type> >  &oStream
  ,const std::basic_string<char_type,traits_type>                &tabIndentStr
  ,const int                                                     startingTab
  ,const bool                                                    showLinePrefix
  ,const int                                                     maxLenLinePrefix
  ,const bool                                                    showTabCount
  )
{
  oStream_ = oStream;
  tabIndentStr_ = tabIndentStr;
  showLinePrefix_ = showLinePrefix;
  maxLenLinePrefix_ = maxLenLinePrefix;
  showTabCount_ = showTabCount;
  tabIndent_ = startingTab;
  tabIndentStack_.resize(0);
  linePrefixStack_.resize(0);
  wroteNewline_ = true;
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
      writeBeginningOfLine();
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
void basic_FancyOStream_buf<CharT,Traits>::writeBeginningOfLine()
{
  bool didOutput = false;
  if(showLinePrefix_) {
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
  for( int i = 0; i < tabIndent_; ++i )
    *oStream_ << tabIndentStr_;
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
 )
  : ostream_t(NULL), streambuf_(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount)
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
  )
{
  streambuf_.initialize(oStream,tabIndentStr,startingTab,showLinePrefix,maxLenLinePrefix,showTabCount);
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
void basic_FancyOStream<CharT,Traits>::setShowLinePrefix(const bool showLinePrefix)
{
  streambuf_.setShowLinePrefix(showLinePrefix);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::setMaxLenLinePrefix(const bool maxLenLinePrefix)
{
  streambuf_.setMaxLenLinePrefix(maxLenLinePrefix);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::setShowTabCount(const bool showTabCount)
{
  streambuf_.setShowTabCount(showTabCount);
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

//
// OSTab
//

} // namespace Teuchos

#endif // TEUCHOS_FANCY_O_STREAM_HPP
