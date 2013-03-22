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

#include "Trios_timer.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <ctime>


double Trios::GetTime()
{
    return trios_get_time();
}

long Trios::GetTimeMS()
{
    return trios_get_time_ms();
}

long Trios::GetTimeNS()
{
    return trios_get_time_ns();
}


long Trios::GetTimeUS()
{
    return trios_get_time_us();
}


int Trios::WriteTimings(
        const std::string &fname,
        const std::string &header,
        const std::vector<std::string> &timings_desc,
        const std::vector<double> &timings)
{
    if (fname.empty())
        return WriteTimings(std::cout, header, timings_desc, timings, true);

    // If the file does not exist, put the complete header
    bool write_header = true;
    std::fstream fin;
    fin.open(fname.c_str(), std::ios::in);
    if (fin.is_open()) {
        write_header = false;
    }
    fin.close();

    // Open the file for writing
    std::ofstream fout;
    fout.open(fname.c_str(), std::ios::out | std::ios::app);
    WriteTimings(fout, header, timings_desc, timings, write_header);
    fout.close();
    return 0;
}

int Trios::WriteTimings(std::ostream &out,
        const std::string &header,
        const std::vector<std::string> &timing_desc,
        const std::vector<double> &timings,
        const bool write_header)
{
    using namespace std;

    string prefix("%");
    string separator(",\t");

    if (write_header) {
        time_t rawtime = time(NULL);
        out << prefix << setfill('-') << setw(52) << " " << std::endl;
        out << prefix << " Timings: " << ctime(&rawtime);
        if (!header.empty())
            out << prefix << header << endl;

        out << prefix << " Column IDS:" << endl;
        for (int i=0; i<(int)timing_desc.size(); i++) {
            out << prefix << setfill(' ') << setw(6) << (i+1) << ": " << timing_desc[i] << endl;
        }
        out << prefix << setfill('-') << setw(52) << " " << endl;
    }

    // just write timing information
    for (int i=0; i<(int)timings.size()-1; i++) {
        out.precision(5);
        out << scientific << timings[i] << separator;
    }
    out << timings[timings.size()-1] << endl;
    return 0;
}
