// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_FANCY_O_STREAM_DEF_HPP
#define TEUCHOS_FANCY_O_STREAM_DEF_HPP

#include "Teuchos_FancyOStream.hpp"

namespace Teuchos {

RCP<basic_FancyOStream<char> >
fancyOStream(
  const RCP< std::basic_ostream<char> >& oStream,
  const std::basic_string<char>& tabIndentStr,
  const int startingTab,
  const bool showLinePrefix,
  const int maxLenLinePrefix,
  const bool showTabCount,
  const bool showProcRank
  )
{
  if (nonnull(oStream)) {
    return rcp(
      new basic_FancyOStream<char>(
        oStream,tabIndentStr,startingTab,showLinePrefix,
        maxLenLinePrefix,showTabCount,showProcRank
        )
      );
  }
  return null;
}


RCP<basic_FancyOStream<char> >
getFancyOStream( const RCP<std::basic_ostream<char> > &out )
{
  if (is_null(out))
    return Teuchos::null;
  RCP<basic_FancyOStream<char> >
    fancyOut = rcp_dynamic_cast<basic_FancyOStream<char> >(out);
  if(nonnull(fancyOut))
    return fancyOut;
  return rcp(new basic_FancyOStream<char>(out));
}

template <typename CharT, typename Traits>
RCP<basic_FancyOStream<CharT,Traits> >
tab(
  const RCP<basic_FancyOStream<CharT,Traits> > &out,
  const int tabs,
  const std::basic_string<CharT,Traits> linePrefix
  )
{
  if(out.get()==NULL)
    return Teuchos::null;
  RCP<basic_FancyOStream<CharT,Traits> > fancyOut = out;
  set_extra_data(
    rcp(new basic_OSTab<CharT,Traits>(out,tabs,linePrefix)),
    "OSTab",
    inOutArg(fancyOut),
    PRE_DESTROY,
    false
    );
  return fancyOut;
}


template <typename CharT, typename Traits>
RCP<basic_FancyOStream<CharT,Traits> >
tab(
  const RCP<std::basic_ostream<CharT,Traits> > &out
  ,const int tabs
  ,const std::basic_string<CharT,Traits> linePrefix
  )
{
  return tab(getFancyOStream(out),tabs,linePrefix);
}


// ////////////////////////////////
// Defintions


//
// basic_FancyOStream_buf
//


template<typename CharT, typename Traits>
basic_FancyOStream_buf<CharT,Traits>::basic_FancyOStream_buf(
  const RCP<std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type> &tabIndentStr
  ,const int startingTab
  ,const bool showLinePrefix
  ,const int maxLenLinePrefix
  ,const bool showTabCount
  ,const bool showProcRank
  )
{
  this->initialize(oStream,tabIndentStr,startingTab,showLinePrefix,
    maxLenLinePrefix,showTabCount,showProcRank);
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::initialize(
  const RCP<std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type> &tabIndentStr
  ,const int startingTab
  ,const bool showLinePrefix
  ,const int maxLenLinePrefix
  ,const bool showTabCount
  ,const bool showProcRank
  )
{
  oStreamSet_ = oStream;
  oStream_ = oStream;
  tabIndentStr_ = tabIndentStr;
  showLinePrefix_ = showLinePrefix;
  maxLenLinePrefix_ = maxLenLinePrefix;
  showTabCount_ = showTabCount;
  showProcRank_ = showProcRank;
  rootRank_ = -1;
  procRank_ = GlobalMPISession::getRank();
  numProcs_ = GlobalMPISession::getNProc();
  rankPrintWidth_ = int(std::log10(float(numProcs_)))+1;
  tabIndent_ = startingTab;
  tabIndentStack_.clear();
  linePrefixStack_.clear();
  wroteNewline_ = true;
  enableTabbingStack_ = 0;
}


template<typename CharT, typename Traits>
RCP<std::basic_ostream<CharT,Traits> >
basic_FancyOStream_buf<CharT,Traits>::getOStream()
{
  return oStreamSet_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setTabIndentStr(
  const std::basic_string<char_type,traits_type> &tabIndentStr
  )
{
  tabIndentStr_ = tabIndentStr;
}


template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream_buf<CharT,Traits>::getTabIndentStr() const
{
  return tabIndentStr_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowLinePrefix(const bool showLinePrefix)
{
  showLinePrefix_ = showLinePrefix;
}


template<typename CharT, typename Traits>
bool basic_FancyOStream_buf<CharT,Traits>::getShowLinePrefix() const
{
  return showLinePrefix_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setMaxLenLinePrefix(const int maxLenLinePrefix)
{
  TEUCHOS_TEST_FOR_EXCEPT( !(maxLenLinePrefix>=5) );
  maxLenLinePrefix_ = maxLenLinePrefix;
}


template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getMaxLenLinePrefix() const
{
  return maxLenLinePrefix_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowTabCount(const bool showTabCount)
{
  showTabCount_ = showTabCount;
}


template<typename CharT, typename Traits>
bool basic_FancyOStream_buf<CharT,Traits>::getShowTabCount() const
{
  return showTabCount_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setShowProcRank(const bool showProcRank)
{
  showProcRank_ = showProcRank;
}


template<typename CharT, typename Traits>
bool basic_FancyOStream_buf<CharT,Traits>::getShowProcRank() const
{
  return showProcRank_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setProcRankAndSize(
  const int procRank, const int numProcs
  )
{
  procRank_ = procRank;
  numProcs_ = numProcs;
}


template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getProcRank() const
{
  return procRank_;
}


template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getNumProcs() const
{
  return numProcs_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::setOutputToRootOnly(
  const int rootRank
  )
{
  rootRank_ = rootRank;
  if(rootRank >= 0) {
    if(rootRank == procRank_)
      oStream_ = oStreamSet_;
    else
      oStream_ = rcp(new oblackholestream());
    // Only processor is being output to so there is no need for line
    // batching!
    lineOut_ = null;
  }
  else {
    oStream_ = oStreamSet_;
    // Output is being sent to all processors so we need line batching to
    // insure that each line will be printed all together!
    lineOut_ = rcp(new std::ostringstream());
  }
}


template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getOutputToRootOnly() const
{
  return rootRank_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushTab(const int tabs)
{
  if( tabIndent_ + tabs < 0 ) {
    tabIndentStack_.push_back(-tabIndent_);
    tabIndent_ = 0;
  }
  else {
    tabIndentStack_.push_back(tabs);
    tabIndent_ += tabs;
  }
}


template<typename CharT, typename Traits>
int basic_FancyOStream_buf<CharT,Traits>::getNumCurrTabs() const
{
  return tabIndent_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::popTab()
{
  tabIndent_ -= tabIndentStack_.back();
  tabIndentStack_.pop_back();
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::pushLinePrefix(
  const std::basic_string<char_type,traits_type> &linePrefix
  )
{
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
std::streamsize basic_FancyOStream_buf<CharT,Traits>::xsputn(
  const char_type* s, std::streamsize n
  )
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
std::ostream& basic_FancyOStream_buf<CharT,Traits>::out()
{
  if(lineOut_.get())
    return *lineOut_;
  return *oStream_;
}


template<typename CharT, typename Traits>
void basic_FancyOStream_buf<CharT,Traits>::writeChars(
  const char_type s[], std::streamsize n
  )
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
      // The last character in the std::string is a newline
      done_outputting = true;
    }
    // Write the beginning of the line if we need to
    if(wroteNewline_) {
      writeFrontMatter();
      wroteNewline_ = false;
    }
    // Write up to the newline or the end of the std::string
    out().write(s+first_p,p-first_p+1);
    if(s[p] == newline) {
      wroteNewline_ = true;
      if(lineOut_.get()) {
        *oStream_ << lineOut_->str() << std::flush;
        lineOut_->str("");
      }
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
  std::ostream &o = this->out();
  if(showProcRank_) {
    o << "p=" << std::right << std::setw(rankPrintWidth_) << procRank_;
    didOutput = true;
  }
  if(showLinePrefix_) {
    if(didOutput)
      o << ", ";
    std::string currLinePrefix = "";
    if ( linePrefixStack_.size() )
      currLinePrefix = this->getTopLinePrefix();
    const int localMaxLenLinePrefix =
      TEUCHOS_MAX( as<int>(currLinePrefix.length()), maxLenLinePrefix_ );
    o << std::left << std::setw(localMaxLenLinePrefix);
    o << currLinePrefix;
    didOutput = true;
  }
  if(showTabCount_) {
    if(didOutput)
      o << ", ";
    o << "tabs=" << std::right << std::setw(2) << tabIndent_;
    didOutput = true;
  }
  // ToDo: Add the Prefix name if asked
  // ToDo: Add the processor number if asked
  // ToDo: Add the number of indents if asked
  if(didOutput) {
    o << " |" << tabIndentStr_;
  }
  if(enableTabbingStack_==0) {
    for( int i = 0; i < tabIndent_; ++i )
      o << tabIndentStr_;
  }
}


//
// basic_FancyOStream
//


template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>::basic_FancyOStream(
  const RCP< std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type> &tabIndentStr
  ,const int startingTab
  ,const bool showLinePrefix
  ,const int maxLenLinePrefix
  ,const bool showTabCount
  ,const bool showProcRank
  )
  :ostream_t(NULL),
   streambuf_(oStream,tabIndentStr,startingTab,showLinePrefix,
     maxLenLinePrefix,showTabCount,showProcRank)
{
  this->init(&streambuf_);
}


template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::initialize(
  const RCP< std::basic_ostream<char_type,traits_type> > &oStream
  ,const std::basic_string<char_type,traits_type> &tabIndentStr
  ,const int startingTab
  ,const bool showLinePrefix
  ,const int maxLenLinePrefix
  ,const bool showTabCount
  ,const bool showProcRank
  )
{
  streambuf_.initialize(oStream,tabIndentStr,startingTab,
    showLinePrefix,maxLenLinePrefix,showTabCount,showProcRank);
  this->init(&streambuf_);
}


template<typename CharT, typename Traits>
RCP<std::basic_ostream<CharT,Traits> >
basic_FancyOStream<CharT,Traits>::getOStream()
{
  return streambuf_.getOStream();
}


template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setTabIndentStr(
  const std::basic_string<char_type,traits_type> &tabIndentStr
  )
{
  streambuf_.setTabIndentStr(tabIndentStr);
  return *this;
}


template<typename CharT, typename Traits>
const std::basic_string<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::getTabIndentStr() const
{
  return streambuf_.getTabIndentStr();
}


template<typename CharT, typename Traits>
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setShowAllFrontMatter(
  const bool showAllFrontMatter
  )
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
basic_FancyOStream<CharT,Traits>::setMaxLenLinePrefix(const int maxLenLinePrefix)
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
basic_FancyOStream<CharT,Traits>&
basic_FancyOStream<CharT,Traits>::setOutputToRootOnly( const int rootRank )
{
  streambuf_.setOutputToRootOnly(rootRank);
  return *this;
}


template<typename CharT, typename Traits>
int basic_FancyOStream<CharT,Traits>::getOutputToRootOnly() const
{
  return streambuf_.getOutputToRootOnly();
}


template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::copyAllOutputOptions(
  const basic_FancyOStream<CharT,Traits> &oStream )
{
  //streambuf_.setTabIndentStr(oStream.streambuf_.getTabIndentStr());
  streambuf_.setShowLinePrefix(oStream.streambuf_.getShowLinePrefix());
  streambuf_.setMaxLenLinePrefix(oStream.streambuf_.getMaxLenLinePrefix());
  streambuf_.setShowTabCount(oStream.streambuf_.getShowTabCount());
  streambuf_.setShowProcRank(oStream.streambuf_.getShowProcRank());
  streambuf_.setProcRankAndSize(oStream.streambuf_.getProcRank(),
    oStream.streambuf_.getNumProcs());
  streambuf_.setOutputToRootOnly(oStream.streambuf_.getOutputToRootOnly());
}


template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushTab(const int tabs)
{
  streambuf_.pushTab(tabs);
}


template<typename CharT, typename Traits>
int basic_FancyOStream<CharT,Traits>::getNumCurrTabs() const
{
  return streambuf_.getNumCurrTabs();
}


template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::popTab()
{
  streambuf_.popTab();
}


template<typename CharT, typename Traits>
void basic_FancyOStream<CharT,Traits>::pushLinePrefix(
  const std::basic_string<char_type,traits_type> &linePrefix
  )
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


} // namespace Teuchos


#endif // TEUCHOS_FANCY_O_STREAM_DEF_HPP
