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
/**  @file pablo_interface.C
 *
 *   @brief Makes calls to the Pablo SDDF library to generate
 *          trace files.
 *
 *   @author Ron Oldfield (raoldfi@sandia.gov)
 *   @version $Revision: 406 $
 *   @date $Date: 2005-10-07 15:08:29 -0600 (Fri, 07 Oct 2005) $
 *
 */

#include "Trios_config.h"

#ifdef HAVE_TRIOS_PABLO

#include "SDDF.h"
#include <AsciiPipeWriter.h>
#include <BinaryPipeWriter.h>
#include <OutputFileStreamPipe.h>
#include <RecordDossier.h>
#include <StructureDescriptor.h>

#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "Trios_timer.h"
#include "Trios_logger.h"
#include "Trios_threads.h"


/*------------ SUPPORT FOR INTERVAL EVENTS WITH PABLO ------------*/

/**
  * A generic event has the following fields
  *
  *	double timestamp;
  *	int id;
  *	int pid;
  *	char *data;
  */


int SDDF::define_generic_event(
        const int tag,
        PipeWriter *writer)
{
    Attributes *attributes = new Attributes();
    FieldDescriptor *fieldP;
    StructureDescriptor *structureP;

    attributes->clearEntries();
    attributes->insert("event", "Generic trace event");
    structureP = new StructureDescriptor("event", *attributes);

    /* double timestamp */
    attributes->clearEntries();
    attributes->insert("timestamp", "Time of event");
    attributes->insert("units", "sec");
    fieldP = new FieldDescriptor("timestamp", *attributes, DOUBLE, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int id */
    attributes->clearEntries();
    attributes->insert("id", "Identifier");
    fieldP = new FieldDescriptor("id", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int pid */
    attributes->clearEntries();
    attributes->insert("pid", "Process/Thread identifier");
    fieldP = new FieldDescriptor("pid", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* char *data */
    attributes->clearEntries();
    attributes->insert("data", "App-specific character data");
    fieldP = new FieldDescriptor("data", *attributes, CHARACTER, 1);
    structureP->insert(*fieldP);
    delete fieldP;

    /* Now we can write the structure descriptor to the output pipe
     * and create the RecordDossier for a generic event.  */
    writer->putDescriptor(*structureP, tag);

    /* create the global genericRecord structure */
    genericRecord = new RecordDossier(tag, *structureP);
    delete structureP;
    delete attributes;

    return 0;
}


/**
  * A count event has is a generic event with an extra count field.
  *
  *	double timestamp;
  *	long count;
  *	int id;
  *	int pid;
  *	char *data;
  */

int SDDF::define_count_event(
        const int tag,
        PipeWriter *writer)
{
    Attributes *attributes = new Attributes();
    FieldDescriptor *fieldP;
    StructureDescriptor *structureP;

    attributes->clearEntries();
    attributes->insert("count", "Generic trace event");
    structureP = new StructureDescriptor("count", *attributes);

    /* double timestamp */
    attributes->clearEntries();
    attributes->insert("timestamp", "Time of event");
    attributes->insert("units", "sec");
    fieldP = new FieldDescriptor("timestamp", *attributes, DOUBLE, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* long count */
    attributes->clearEntries();
    attributes->insert("count", "Count");
    fieldP = new FieldDescriptor("count", *attributes, LONG, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int id */
    attributes->clearEntries();
    attributes->insert("id", "Identifier");
    fieldP = new FieldDescriptor("id", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int pid */
    attributes->clearEntries();
    attributes->insert("pid", "Process/Thread identifier");
    fieldP = new FieldDescriptor("pid", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* char *data */
    attributes->clearEntries();
    attributes->insert("data", "App-specific character data");
    fieldP = new FieldDescriptor("data", *attributes, CHARACTER, 1);
    structureP->insert(*fieldP);
    delete fieldP;

    /* Now we can write the structure descriptor to the output pipe
     * and create the RecordDossier for a generic event.  */
    writer->putDescriptor(*structureP, tag);

    /* create the global countRecord structure */
    countRecord = new RecordDossier(tag, *structureP);

    delete structureP;
    delete attributes;

    return 0;
}


/**
  * A count event has is a generic event with an extra count field.
  *
  *	double timestamp;
  *	int id;
  *	int pid;
  *	int level;
  *     double duration;
  *	char *data;
  */
int SDDF::define_interval_event(
        const int tag,
        PipeWriter *writer)
{
    Attributes *attributes = new Attributes();
    FieldDescriptor *fieldP;
    StructureDescriptor *structureP;

    attributes->clearEntries();
    attributes->insert("interval", "Interval event");
    structureP = new StructureDescriptor("interval", *attributes);

    /* double timestamp */
    attributes->clearEntries();
    attributes->insert("timestamp", "Time of event");
    attributes->insert("units", "sec");
    fieldP = new FieldDescriptor("timestamp", *attributes, DOUBLE, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int id */
    attributes->clearEntries();
    attributes->insert("id", "Identifier");
    fieldP = new FieldDescriptor("id", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int pid */
    attributes->clearEntries();
    attributes->insert("pid", "Process/Thread identifier");
    fieldP = new FieldDescriptor("pid", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int level */
    attributes->clearEntries();
    attributes->insert("level", "Level/Scope");
    fieldP = new FieldDescriptor("level", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* double duration */
    attributes->clearEntries();
    attributes->insert("duration", "Length of interval");
    attributes->insert("units", "sec");
    fieldP = new FieldDescriptor("duration", *attributes, DOUBLE, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* char *data */
    attributes->clearEntries();
    attributes->insert("data", "App-specific character data");
    fieldP = new FieldDescriptor("data", *attributes, CHARACTER, 1);
    structureP->insert(*fieldP);
    delete fieldP;

    /* Now we can write the structure descriptor to the output pipe
     * and create the RecordDossier for a generic event.  */
    writer->putDescriptor(*structureP, tag);

    /* create the global intervalRecord structure */
    intervalRecord = new RecordDossier(tag, *structureP);

    delete structureP;
    delete attributes;

    return 0;
}


/**
  * A throughput event has the following fields
  *
  *	double timestamp;
  *	int id;
  *	int pid;
  *	int level;
  *	char *data;
  *	double duration;
  *	long count;
  */

int SDDF::define_throughput_event(
        const int tag,
        PipeWriter *writer)
{
    Attributes *attributes = new Attributes();
    FieldDescriptor *fieldP;
    StructureDescriptor *structureP;

    attributes->clearEntries();
    attributes->insert("throughput", "Throughput Event");
    structureP = new StructureDescriptor("throughput", *attributes);

    /* double timestamp */
    attributes->clearEntries();
    attributes->insert("timestamp", "Time of event");
    attributes->insert("units", "sec");
    fieldP = new FieldDescriptor("timestamp", *attributes, DOUBLE, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int id */
    attributes->clearEntries();
    attributes->insert("id", "Identifier");
    fieldP = new FieldDescriptor("id", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int pid */
    attributes->clearEntries();
    attributes->insert("pid", "Process/Thread identifier");
    fieldP = new FieldDescriptor("pid", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* int level */
    attributes->clearEntries();
    attributes->insert("level", "Level/Scope");
    fieldP = new FieldDescriptor("level", *attributes, INTEGER, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* double duration */
    attributes->clearEntries();
    attributes->insert("duration", "Length of interval");
    attributes->insert("units", "sec");
    fieldP = new FieldDescriptor("duration", *attributes, DOUBLE, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* long count */
    attributes->clearEntries();
    attributes->insert("count", "Number of objects processed");
    fieldP = new FieldDescriptor("count", *attributes, LONG, 0);
    structureP->insert(*fieldP);
    delete fieldP;

    /* char *data */
    attributes->clearEntries();
    attributes->insert("data", "App-specific character data");
    fieldP = new FieldDescriptor("data", *attributes, CHARACTER, 1);
    structureP->insert(*fieldP);
    delete fieldP;

    /* Now we can write the structure descriptor to the output pipe
     * and create the RecordDossier for a generic event.  */
    writer->putDescriptor(*structureP, tag);

    /* create the global throughputRecord structure */
    throughputRecord = new RecordDossier(tag, *structureP);

    delete structureP;
    delete attributes;

    return 0;
}


/**
 * @brief Initialize the pablo output interface.
 *
 * @param fname @input_type  The name of the outputfile.
 * @param type  @input_type  Type of file (0=binary, 1=ascii)
 *
 */
SDDF::SDDF(const char *f, const int t):
    TraceFile(f), ftype(t), genericRecord(0), countRecord(0),
    intervalRecord(0), throughputRecord(0)
{
    static const int bufsize = 204800;

    /* get the current time */
    time_t now = time(NULL);

    /* store the start time */
    starttime = trios_get_time();

    Attributes attributes;


    /* initialize the mutex for the genericRecord */
    nthread_lock_init(&genericMutex);
    /* initialize the mutex for the countRecord  */
    nthread_lock_init(&countMutex);
    /* initialize the mutex for the intervalRecord  */
    nthread_lock_init(&intervalMutex);
    /* initialize the mutex for the throughputRecord  */
    nthread_lock_init(&throughputMutex);

    /* initialize the mutex for the output stream  */
    nthread_lock_init(&outputMutex);


    /* Open file */
    outFile = new OutputFileStreamPipe(fname, bufsize);
    if (ftype) {
        pipeWriter = new AsciiPipeWriter(outFile);
    }
    else {
        pipeWriter = new BinaryPipeWriter(outFile);
    }

    /* Stream Attribute */
    attributes.clearEntries();
    attributes.insert("run date", ctime(&now));
    /* ... what else goes in the header? */
    pipeWriter->putAttributes(attributes);

    //output_header(pipeWriter);

    /* ---- Describe the types of records we expect ---- */
    define_generic_event(GENERIC_RECORD, pipeWriter);
    define_count_event(COUNT_RECORD, pipeWriter);
    define_interval_event(INTERVAL_RECORD, pipeWriter);
    define_throughput_event(THROUGHPUT_RECORD, pipeWriter);

    output_generic_event(0,0,"init");
    output_count_event(0,0,"init",0);
    output_interval_event(0,0,0,"init",0);
    output_tput_event(0,0,0,"init",0,0);
}

SDDF::~SDDF(void)
{
    fprintf(stderr, "delete outfile\n");
    delete outFile;
    fprintf(stderr, "delete pipewriter\n");
    delete pipeWriter;
    fprintf(stderr, "delete genericRec\n");
    delete genericRecord;
    fprintf(stderr, "delete countRecord\n");
    delete countRecord;
    fprintf(stderr, "delete throughputRecord\n");
    delete throughputRecord;
    fprintf(stderr, "delete intervalRecord\n");
    delete intervalRecord;

    /* cleanup the mutex for the genericRecord */
    nthread_lock_fini(&genericMutex);
    /* cleanup the mutex for the countRecord  */
    nthread_lock_fini(&countMutex);
    /* cleanup the mutex for the intervalRecord  */
    nthread_lock_fini(&intervalMutex);
    /* cleanup the mutex for the throughputRecord  */
    nthread_lock_fini(&throughputMutex);

    /* cleanup the mutex for the output stream  */
    nthread_lock_fini(&outputMutex);


    fprintf(stderr, "finished ~SDDF()\n");
}

int SDDF::output_record(RecordDossier *rec)
{
    int rc = 0;

    if (!rec)
        return rc;
    /* need to protect the putData function */
    nthread_lock(&mutex);
    pipeWriter->putData(*rec);
    nthread_unlock(&mutex);

    return rc;
}

/**
 * @brief Output a generic trace event.
 *
 * @param eventID @input_type  The ID of the trace.
 * @param pid     @input_type  Process ID.
 * @param data    @input_type  User-defined data passed in a character string.
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int SDDF::output_generic_event(
        const int eventID,
        const int pid,
        const char *data)
{
    int rc = 1;

    nthread_lock(&genericMutex);

    genericRecord->setValue("timestamp", trios_get_time() - starttime);
    genericRecord->setValue("id", eventID);
    genericRecord->setValue("pid", pid);

    /*fprintf(stderr, "setting data to \"%s\"\n",data);*/
    genericRecord->setCString("data", data);

    nthread_unlock(&genericMutex);

    output_record(genericRecord);

    return rc;
}



/**
 * @brief Output an interval event.
 *
 * We use the generic pablo trace and encode the additional information
 * we want in the data field.  The new data field will be, "interval:$name:duration".
 *
 * Pablo has its own interval records, but they are inadequate because it is
 * difficult to measure concurrent intervals (e.g., in threads).
 * A better (and more efficient) way to do this would be to create our own
 * Pablo record type, but this is a quick "hack" to address our needs.
 *
 * @param eventID @input_type  The ID of the trace.
 * @param pid     @input_type  Process ID.
 * @param level   @input_type
 * @param data    @input_type  User-defined data passed in a character
 *                             string.
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int SDDF::output_interval_event(
        const int eventID,
        const int pid,
        const int level,
        const char *data,
        double duration)
{
    int rc = 1;

    nthread_lock(&intervalMutex);
    intervalRecord->setValue("timestamp", trios_get_time() - starttime);
    intervalRecord->setValue("id", eventID);
    intervalRecord->setValue("pid", pid);
    intervalRecord->setValue("level", level);
    intervalRecord->setValue("duration", duration);

    /*fprintf(stderr, "setting data to \"%s\"\n",data);*/
    intervalRecord->setCString("data", data);
    nthread_unlock(&intervalMutex);

    output_record(intervalRecord);

    return rc;
}

int SDDF::output_tput_event(
        const int eventID,
        const int pid,
        const int level,
        const char *data,
        double duration,
        const long count)
{
    int rc = 1;

    nthread_lock(&throughputMutex);
    throughputRecord->setValue("timestamp", trios_get_time() - starttime);
    throughputRecord->setValue("id", eventID);
    throughputRecord->setValue("pid", pid);
    throughputRecord->setValue("level", level);
    throughputRecord->setValue("duration", duration);
    throughputRecord->setValue("count", count);
    nthread_unlock(&throughputMutex);

    /*fprintf(stderr, "setting data to \"%s\"\n",data);*/
    throughputRecord->setCString("data", data);

    output_record(throughputRecord);

    return rc;
}


/**
 * @brief Output a count event.
 *
 * We use the generic pablo trace and encode the additional information
 * we want in the data field.  The new data field will be, "interval:$name:duration".
 *
 * Pablo has its own count records, but they are inadequate because they only
 * increment values.  We want to increment, decrement, and set count events.
 * A better (and more efficient) way to do this would be to create our own
 * Pablo record type, but this is a quick "hack" that will still work.
 *
 * @param intervalID @input_type The interval ID (unique for each interval).
 * @param eventID @input_type  The ID of the trace.
 * @param data    @input_type  User-defined data passed in a character
 *                             string.
 *
 * @return non-zero if successfull
 * @return 0 if failure
 */
int SDDF::output_count_event(
        const int eventID,
        const int pid,
        const char *data,
        const int count)
{
    int rc = 1;

    nthread_lock(&countMutex);
    countRecord->setValue("timestamp", trios_get_time() - starttime);
    countRecord->setValue("id", eventID);
    countRecord->setValue("pid", pid);
    countRecord->setValue("count", count);

    /*fprintf(stdout, "setting count data to \"%s\"\n",data);*/
    countRecord->setCString("data", data);
    nthread_unlock(&countMutex);

    output_record(countRecord);

    return rc;
}

int SDDF::set_buffer_size(const unsigned long int bufsize)
{
    return -1;

}

#endif // HAVE_TRIOS_PABLO
