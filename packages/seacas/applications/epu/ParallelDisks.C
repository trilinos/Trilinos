/*
 * Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */

#include <smart_assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <string>

#include <ParallelDisks.h>

#include <to_string.h>

/*****************************************************************************/
Excn::ParallelDisks::ParallelDisks()
  : disk_names(NULL), number_of_raids(0), raid_offset(0)
{}

/*****************************************************************************/
Excn::ParallelDisks::~ParallelDisks()
{
  delete [] disk_names;
}

/*****************************************************************************/
void Excn::ParallelDisks::Number_of_Raids(int i)
{
  number_of_raids = i;
  create_disk_names();
}

/*****************************************************************************/
void Excn::ParallelDisks::Raid_Offset(int i)
{
  raid_offset = i;
  create_disk_names();
}

/*****************************************************************************/
int Excn::ParallelDisks::Number_of_Raids() const
{
  return number_of_raids;
}
 
/*****************************************************************************/
int Excn::ParallelDisks::Raid_Offset() const
{
  return raid_offset;
}


/*****************************************************************************/
void Excn::ParallelDisks::rename_single_file_for_mp(const std::string& curdir,
					      std::string& name) const
{
  // Possible to have node layout without parallel disks 
  if(curdir.length() && name[0] != '/') {
    name = curdir + "/" + name;
  }
  return;
}

/*****************************************************************************/
void Excn::ParallelDisks::rename_file_for_mp(const std::string& rootdir,
				       const std::string& subdir, 
				       std::string& name,
				       int node, int numproc) const
{
  // Possible to have node layout without parallel disks 

  std::string prepend;
  if (rootdir.length()) {
    prepend = rootdir + "/";
  } else if (name[0]=='/') {
    prepend = "";
  } else {
    prepend = "./";
  }
  
  int lnn = node;
  if(number_of_raids) {
    int diskn = lnn % number_of_raids;
    Create_IO_Filename(name, lnn, numproc);
    name = disk_names[diskn] + "/" + subdir +  "/" + name;
  }
  else {
    Create_IO_Filename(name, lnn, numproc);
    if (subdir.length()) {
      name = subdir + "/" + name;
    }
  }
  name = prepend + name;
  return;
}

/*****************************************************************************/
void Excn::ParallelDisks::create_disk_names()
{

  if(!number_of_raids) return;

  delete [] disk_names;
  disk_names = new std::string[number_of_raids];

  char str[3];

  for(int i = 0; i < number_of_raids; i++) {
    int num = i + raid_offset;
    if(num < 10) {
#ifdef COUGAR
      sprintf(str, "%d", num);
#else
      sprintf(str, "0%d", num);
#endif
      disk_names[i] = std::string(str);
    }
    else {
      sprintf(str, "%d", num); 
      disk_names[i] = std::string(str);
    } 
  }
}

/*****************************************************************************/
void Excn::ParallelDisks::Create_IO_Filename(std::string &name,
				       int lnn, int num_proc )
{
  if (num_proc >= 10000000) {
    std::cerr << "ERROR: EPU is currently limited to less than 10,000,000 processors.\n"
	      << "       Please contact gdsjaar@sandia.gov to increase the limit.\n";
    exit(EXIT_FAILURE);
  }
  
  name = name + "." + ToString(num_proc);
  if(num_proc >= 1000000) {
    if      (lnn <       10) { name = name + ".000000" + ToString(lnn); }
    else if (lnn <      100) { name = name + ".00000" + ToString(lnn); }
    else if (lnn <     1000) { name = name + ".0000" + ToString(lnn); }
    else if (lnn <    10000) { name = name + ".000" + ToString(lnn); }
    else if (lnn <   100000) { name = name + ".00" + ToString(lnn); }
    else if (lnn <  1000000) { name = name + ".0" + ToString(lnn); }
    else if (lnn < 10000000) { name = name + "." + ToString(lnn); }
  }
  if(num_proc >= 100000) {
    if      (lnn <      10) { name = name + ".00000" + ToString(lnn); }
    else if (lnn <     100) { name = name + ".0000" + ToString(lnn); }
    else if (lnn <    1000) { name = name + ".000" + ToString(lnn); }
    else if (lnn <   10000) { name = name + ".00" + ToString(lnn); }
    else if (lnn <  100000) { name = name + ".0" + ToString(lnn); }
    else if (lnn < 1000000) { name = name + "." + ToString(lnn); }
  }
  if(num_proc >= 10000) {
    if      (lnn <     10) { name = name + ".0000" + ToString(lnn); }
    else if (lnn <    100) { name = name + ".000" + ToString(lnn); }
    else if (lnn <   1000) { name = name + ".00" + ToString(lnn); }
    else if (lnn <  10000) { name = name + ".0" + ToString(lnn); }
    else if (lnn < 100000) { name = name + "." + ToString(lnn); }
  }
  else if(num_proc >= 1000) {
    if      (lnn <    10)  { name = name + ".000" + ToString(lnn); }
    else if (lnn <   100)  { name = name + ".00" + ToString(lnn); }
    else if (lnn <  1000)  { name = name + ".0" + ToString(lnn); }
    else if (lnn < 10000)  { name = name + "." + ToString(lnn); }
  }
  else if(num_proc >= 100) {
    if      (lnn <   10)   { name = name + ".00" + ToString(lnn); }
    else if (lnn <  100)   { name = name + ".0" + ToString(lnn); }
    else if (lnn < 1000)   { name = name + "." + ToString(lnn); }
  }
  else if(num_proc >= 10) {
    if      (lnn <  10)    { name = name + ".0" + ToString(lnn); }
    else if (lnn < 100)    { name = name + "." + ToString(lnn); }
  }
  else if(num_proc >= 1) {
    name=name + "." + ToString(lnn);
  }
  return;
}
