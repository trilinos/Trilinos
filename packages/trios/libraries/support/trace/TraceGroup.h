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
#ifndef TRACE_GROUP_H_
#define TRACE_GROUP_H_

#include <algorithm>
#include <list>

using namespace std;
using namespace __gnu_cxx;




/**
 * @brief The data structure used to manage trace groups.
 */
class TraceGroup {

    public:
        TraceGroup *parent;
        const char *name;
        const int gid;
        bool enabled;

        /* Children of this group */
        list<TraceGroup *> children;

    public:
        TraceGroup(const char *s, const int id, bool enable_flag):
            parent(0), name(s), gid(id), enabled(enable_flag), children() { }

        ~TraceGroup() {}

        /**
         * @brief Set the enabled flag on this group and all its children.
         */
        int enable() {
            enabled = true;

            /* enable each of the children */
            //for_each(children.begin(), children.end(), enable_group);

            return 0;
        }

        /**
         * @brief Enable the specified TraceGroup.
         */
        static void enable_group(TraceGroup *grp) {
            grp->enable();
        }

        static void disable_group(TraceGroup *g) {
            g->disable();
        }

        struct disabler {
            void operator() (TraceGroup *g) {
                g->disable();
            }
        };

        /**
         * @brief Set the disabled flag on this group and all its children.
         */
        int disable() {
            enabled = false;

            /* enable each of the children */
            for_each(children.begin(), children.end(), disabler());
            return 0;
        }


        /**
         * @brief Add a sub-group to this group.
         */
        int add_child(TraceGroup *child) {
            children.push_back(child);
            if (!child->parent) {
                child->parent = this;
            }
            return 0;
        }
};



#endif /*TRACE_GROUP_H*/
