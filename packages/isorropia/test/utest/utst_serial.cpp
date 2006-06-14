//@HEADER
/*
************************************************************************

              Isorropia: Partitioning and Load Balancing Package
                Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
USA
Questions? Contact Alan Williams (william@sandia.gov)
                or Erik Boman    (egboman@sandia.gov)

************************************************************************
*/
//@HEADER

//This file is organized as follows:
//- include-directives
//- unit-test functions, each function unit-tests a single capability
//  in either Isorropia itself, or the Isorropia test-utilities.
//- function 'run_serial_utests', which calls the above unit-test functions.

//---------------------------------------------------------------
//Include-directives, to include relevant Isorropia declarations 
//
#include <Isorropia_configdefs.hpp>
#include <Isorropia_Exception.hpp>
#include <Isorropia_Utils.hpp>

#include <ispatest_read_distribution.hpp>

//---------------------------------------------------------------
//Unit-test functions, one function per unit-test.
//
void test_create_comm_plan(bool verbose)
{
  if (verbose) {
    std::cout << "testing Isorropia::Utils::create_comm_plan... ";
  }

  std::vector<int> all_old_offsets(4);
  std::vector<int> all_new_offsets(4);

  all_old_offsets[0] = 0;
  all_old_offsets[1] = 5;
  all_old_offsets[2] = 9;
  all_old_offsets[3] = 15;

  all_new_offsets[0] = 0;
  all_new_offsets[1] = 7;
  all_new_offsets[2] = 12;
  all_new_offsets[3] = 15;

  std::vector<int> send_info;
  std::vector<int> recv_info;

  //this is a serial test, but imagine we were proc 2...
  int myPID = 2;
  Isorropia::Utils::create_comm_plan(myPID, all_old_offsets, all_new_offsets,
                                     send_info, recv_info);

  //from the above old/new offset data, we expect that proc 2 would
  //send 3 elements, starting at offset 0, to proc 1, and that's all.
  //So we expect recv_info.size()==0 and send_info.size()==3.
  if (recv_info.size() != 0) {
    throw Isorropia::Exception("create_comm_plan test 1a failed.");
  }
  if (send_info.size() != 3) {
    throw Isorropia::Exception("create_comm_plan test 1b failed.");
  }

  //We expect send_info's contents to be the following:
  if (send_info[0] != 1 || send_info[1] != 0 || send_info[2] != 3) {
    throw Isorropia::Exception("create_comm_plan test 1c failed.");
  }

  //now imagine we were proc 1...
  myPID = 1;
  Isorropia::Utils::create_comm_plan(myPID, all_old_offsets, all_new_offsets,
                                     send_info, recv_info);

  //from the old/new offset data, we expect that proc 1 would send
  //2 elements, starting at offset 0, to proc 0.
  //we also expect that proc 1 would recv 3 elements from proc 2,
  //to be stored at offset 2 in the new local element list.
  if (recv_info.size() != 3 || send_info.size() != 3) {
    throw Isorropia::Exception("create_comm_plan test 2a failed.");
  }

  if (recv_info[0] != 2 || recv_info[1] != 2 || recv_info[2] != 3) {
    throw Isorropia::Exception("create_comm_plan test 2b failed.");
  }

  if (send_info[0] != 0 || send_info[1] != 0 || send_info[2] != 2) {
    throw Isorropia::Exception("create_comm_plan test 2c failed.");
  }

  all_old_offsets[0] = 0;
  all_old_offsets[1] = 2;
  all_old_offsets[2] = 4;
  all_old_offsets[3] = 15;

  all_new_offsets[0] = 0;
  all_new_offsets[1] = 5;
  all_new_offsets[2] = 10;
  all_new_offsets[3] = 15;

  //now imagine we were proc 1...
  myPID = 1;
  Isorropia::Utils::create_comm_plan(myPID, all_old_offsets, all_new_offsets,
                                     send_info, recv_info);

  //from the old/new offset data, we expect that proc 1 would send
  //2 elements, starting at offset 0, to proc 0.
  //we also expect that proc 1 would recv 5 elements from proc 2,
  //to be stored at offset 0 in the new local element list.
  if (recv_info.size() != 3 || send_info.size() != 3) {
    throw Isorropia::Exception("create_comm_plan test 3a failed.");
  }

  if (recv_info[0] != 2 || recv_info[1] != 0 || recv_info[2] != 5) {
    throw Isorropia::Exception("create_comm_plan test 3b failed.");
  }

  if (send_info[0] != 0 || send_info[1] != 0 || send_info[2] != 2) {
    throw Isorropia::Exception("create_comm_plan test 3c failed.");
  }

  //now imagine we were proc 2...
  myPID = 2;
  Isorropia::Utils::create_comm_plan(myPID, all_old_offsets, all_new_offsets,
                                     send_info, recv_info);

  //from the old/new offset data, we expect that proc 2 would send
  //1 element, starting at offset 0, to proc 0.
  //we also expect that proc 2 would send 5 elements, starting at
  //offset 1, to proc 1.
  //we expect that proc 2 would not recv any elements.
  if (recv_info.size() != 0 || send_info.size() != 6) {
    throw Isorropia::Exception("create_comm_plan test 4a failed.");
  }

  if (send_info[0] != 0 || send_info[1] != 0 || send_info[2] != 1) {
    throw Isorropia::Exception("create_comm_plan test 4b failed.");
  }

  if (send_info[3] != 1 || send_info[4] != 1 || send_info[5] != 5) {
    throw Isorropia::Exception("create_comm_plan test 4c failed.");
  }

  if (verbose) {
    std::cout << "ok" << std::endl;
  }
}

void test_read_distribution(bool verbose, const char* path)
{
  if (verbose) {
    std::cout << "testing ispatest::read_distribution... ";
  }

  std::vector<int> rows1;
  std::vector<int> cols1;
  std::vector<int> partitions1;

  std::string fullname(path);
  fullname += "/utest_dist1";
  ispatest::read_distribution(fullname.c_str(), rows1, cols1, partitions1);

  std::vector<int> rows2;
  std::vector<int> cols2;
  std::vector<int> partitions2;

  fullname = path;
  fullname += "/utest_dist2";
  ispatest::read_distribution(fullname.c_str(), rows2, cols2, partitions2);

  //We know that the utest_dist1 file should contain 6
  //row-col-partition triplets. So we'll check the vectors to see if
  //they all have size 6.
  if (rows1.size() != 6 || cols1.size() != 6 || partitions1.size() != 6) {
    throw Isorropia::Exception("test_read_distribution check1 failed");
  }

  //We also know that the contents of utest_dist1 and
  //utest_dist2 should be the same.
  if (rows1 != rows2 || cols1 != cols2 ||
      partitions1 != partitions2) {
    throw Isorropia::Exception("test_read_distribution check2 failed");
  }

  std::vector<int> rows3;
  std::vector<int> cols3;
  std::vector<int> partitions3;

  fullname = path;
  fullname += "/utest_dist3";
  ispatest::read_distribution(fullname.c_str(), rows3, cols3, partitions3);

  //We know that utest_dist3 should contain the same rows and cols as
  //utest_dist1 and utest_dist2, but different partitions.
  if (rows3 != rows2 || cols3 != cols1 ||
      partitions3 == partitions1) {
    throw Isorropia::Exception("test_read_distribution check3 failed");
  }

  if (verbose) {
    std::cout << "ok"<<std::endl;
  }
}

//---------------------------------------------------------------
//The following function, run_serial_utests, is kind of like a 'main',
//it calls each of the above unit-test functions.
//
void run_serial_utests(bool verbose)
{
  if (verbose) {
#ifdef HAVE_EPETRA
    std::cout << "run_serial_utests: HAVE_EPETRA defined."<<std::endl;
#else
    std::cout << "run_serial_utests: HAVE_EPETRA not defined."<<std::endl;
#endif
  }

  test_create_comm_plan(verbose);

  try {
    test_read_distribution(verbose, ".");
  }
  catch(...) {
    test_read_distribution(verbose, "./utest");
  }
}

