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
#ifndef _SDDF_H_
#define _SDDF_H_

#include "Trios_config.h"

#ifdef HAVE_TRIOS_PABLO

#include "TraceFile.h"
#include <AsciiPipeWriter.h>
#include <BinaryPipeWriter.h>
#include <RecordDossier.h>
#include <OutputFileStreamPipe.h>

#include "Trios_threads.h"


/**
 * @brief A class for outputting SDDF trace files
 *        using Pablo.
 */
class SDDF : public TraceFile {

    private:
        int ftype;

        enum RecordTags {
            GENERIC_RECORD=1,
            COUNT_RECORD,
            INTERVAL_RECORD,
            THROUGHPUT_RECORD
        };

        /* The RecordDossier holds information about the
         * structure of a record.  */
        RecordDossier *genericRecord;
        RecordDossier *countRecord;
        RecordDossier *intervalRecord;
        RecordDossier *throughputRecord;

        /* need mutexes to protect the RecordDossier */
        nthread_lock_t genericMutex;
        nthread_lock_t countMutex;
        nthread_lock_t intervalMutex;
        nthread_lock_t throughputMutex;

        /* Where to put the data */
        OutputFileStreamPipe *outFile;
        PipeWriter *pipeWriter;

        /* need mutex to protect the output stream */
        nthread_lock_t outputMutex;

        /* When the trace file is initialized */
        double starttime;

    public:
        SDDF(const char *f, const int t);

        ~SDDF();

        int output_record(RecordDossier *rec);

        int output_generic_event(
                const int eventID,
                const int pid,
                const char *data);

        int output_interval_event(
                const int eventID,
                const int pid,
                const int level,
                const char *data,
                double duration);

        int output_tput_event(
                const int eventID,
                const int pid,
                const int level,
                const char *data,
                double duration,
                const long num_processed);

        int output_count_event(
                const int eventID,
                const int pid,
                const char *data,
                const int count);

        int set_buffer_size(
                const unsigned long int bufsize);

    private:
        int define_generic_event(
                const int tag,
                PipeWriter *pipeWriter);

        int define_count_event(
                const int tag,
                PipeWriter *pipeWriter);

        int define_interval_event(
                const int tag,
                PipeWriter *pipeWriter);

        int define_throughput_event(
                const int tag,
                PipeWriter *pipeWriter);

};


#endif /* HAVE_TRIOS_PABLO */

#endif /* _SDDF_H_*/
