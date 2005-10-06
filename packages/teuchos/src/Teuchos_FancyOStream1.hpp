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
#include "Teuchos_OSTab1.hpp"

namespace Teuchos {

/** \brief Stream buffering class that performs the magic of indenting
 * data sent to an std::ostream object.
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
    const Teuchos::RefCountPtr<std::basic_ostream<char_type,traits_type> >  &out
    ,const std::basic_string<char_type,traits_type>                         &tabIndent
    ,const int                                                              tabTag
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr<std::basic_ostream<char_type,traits_type> >  &out
    ,const std::basic_string<char_type,traits_type>                         &tabIndent
    ,const int                                                              tabTag
    );
  
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
  
  Teuchos::RefCountPtr<std::basic_ostream<char_type,traits_type> >  out_;
  std::basic_string<char_type,traits_type>                          tabIndent_;
  int                                                               tabTag_;

  bool                                                              wroteNewline_;

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
    const Teuchos::RefCountPtr< std::basic_ostream<char_type,traits_type> > &out
    ,const std::basic_string<char_type,traits_type>                         &tabIndent = "\t"
    ,const int                                                              tabTag     = 0
    );

  /** \brief . */
  void initialize(
    const Teuchos::RefCountPtr< std::basic_ostream<char_type,traits_type> > &out
    ,const std::basic_string<char_type,traits_type>                         &tabIndent = "\t"
    ,const int                                                              tabTag     = 0
    );

private:

  streambuf_t	streambuf_;

};

/** \brief . */
typedef basic_FancyOStream<char> FancyOStream;

// ////////////////////////////////
// Defintions

// basic_FancyOStream_buf

template<typename CharT, typename Traits>
basic_FancyOStream_buf<CharT,Traits>::basic_FancyOStream_buf(
  const Teuchos::RefCountPtr<std::basic_ostream<char_type,traits_type> >  &out
  ,const std::basic_string<char_type,traits_type>                         &tabIndent
  ,const int                                                              tabTag
  )
{
  this->initialize(out,tabIndent,tabTag);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::initialize(
  const Teuchos::RefCountPtr<std::basic_ostream<char_type,traits_type> >  &out
  ,const std::basic_string<char_type,traits_type>                         &tabIndent
  ,const int                                                              tabTag
  )
{
  out_ = out;
  tabIndent_ = tabIndent;
  tabTag_ = tabTag;
  wroteNewline_ = true;
}

template<typename CharT, typename Traits>
std::streamsize basic_FancyOStream_buf<CharT,Traits>::xsputn(const char_type* s, std::streamsize n)
{
#ifdef TEUCHOS_FANCY_OSTREAM_SHOW_ALL_CALLS
  std::cerr << "\ncalled xsputn()\n";
#endif
  if(n == 0) return 0;
  std::streamsize p = 0, first_p = 0;
  bool done_outputting = false;
  while( !done_outputting ) {
    // Find the next newline
    for( p = first_p; p < n; ++p ) {
      if(s[p] == '\n') {
        break;
      }
    }
    if(p == n) {
      // We did not find a newline at the end!
      --p;
      done_outputting = true;
    }
    else if( p == n-1 && s[p] == '\n' ) {
      // The last character in the string is a newline
      done_outputting = true;
    }
    // Write the beginning of the line if we need to
    if(wroteNewline_) {
      writeBeginningOfLine();
      wroteNewline_ = false;
    }
    // Write up to the newline or the end of the string
    out_->write(s+first_p,p-first_p+1);
    if(s[p] == '\n') {
      wroteNewline_ = true;
    }
    // Update for next search
    if(!done_outputting)
      first_p = p+1; 
  }
  return n; // Is this correct???
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
    this->xsputn(cc,1);
  }
  return traits_type::not_eof(c);
  //return std::basic_streambuf<CharT,Traits>::overflow(c);
}

template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::writeBeginningOfLine()
{
  // ToDo: Add the Prefix name if asked
  // ToDo: Add the processor number if asked
  // ToDo: Add the number of indents if asked
  const int indent = OSTab::getCurrIndent(tabTag_);
  for( int i = 0; i < indent; ++i )
    *out_ << tabIndent_;
}

// basic_FancyOStream

template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>::basic_FancyOStream(
  const Teuchos::RefCountPtr< std::basic_ostream<char_type,traits_type> > &out
  ,const std::basic_string<char_type,traits_type>                         &tabIndent
  ,const int                                                              tabTag
  )
  : ostream_t(NULL), streambuf_(out,tabIndent,tabTag)
{
  this->init(&streambuf_);
}

template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::initialize(
  const Teuchos::RefCountPtr< std::basic_ostream<char_type,traits_type> > &out
  ,const std::basic_string<char_type,traits_type>                         &tabIndent
  ,const int                                                              tabTag
  )
{
  streambuf_.initialize(out,tabIndent,tabTag);
  this->init(&streambuf_);
}

} // namespace Teuchos

#endif // TEUCHOS_FANCY_O_STREAM_HPP
