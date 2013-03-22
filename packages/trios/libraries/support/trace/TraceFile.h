/**
//@HEADER
// ************************************************************************
//
//                   Trios: Trilinos I/O Support
//                 Copyright 2011 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//Questions? Contact Ron A. Oldfield (raoldfi@sandia.gov)
//
// *************************************************************************
//@HEADER
 */
#ifndef TRACE_FILE_H_
#define TRACE_FILE_H_


/**
 * @brief An abstract base class for trace files.
 */
class TraceFile {

    public:
        const char *fname;

    public:
        TraceFile(const char *f): fname(f) { }

        virtual ~TraceFile() {}

        /**
         * @brief Return the file name.
         */
        const char *get_fname() {return fname;}

        virtual int output_generic_event(
                const int eventID,
                const int pid,
                const char *data) = 0;

        virtual int output_interval_event(
                const int eventID,
                const int pid,
                const int level,
                const char *data,
                double duration) = 0;

        virtual int output_tput_event(
                const int eventID,
                const int pid,
                const int level,
                const char *data,
                double duration,
                const long num_processed) = 0;

        virtual int output_count_event(
                const int eventID,
                const int pid,
                const char *data,
                const int count) = 0;

        virtual int set_buffer_size(
                const unsigned long int bufsize) = 0;
};



#endif /*TRACE_FILE_H*/
