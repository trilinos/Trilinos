// Copyright(C) 1999-2010
// Sandia Corporation. Under the terms of Contract
// DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
// certain rights in this software.
//         
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <Ioss_Field.h>
#include <Ioss_VariableType.h>
#include <Ioss_Transform.h>
#include <Ioss_Utils.h>

#include <assert.h>
#include <iostream>
#include <string>

namespace {
  size_t internal_get_size(Ioss::Field::BasicType type, size_t count,
			   const Ioss::VariableType* storage);
  
  std::string type_string(Ioss::Field::BasicType type)
  {
    switch (type) {
    case Ioss::Field::REAL:
      return std::string("real");
    case Ioss::Field::INTEGER:
      return std::string("integer");
    case Ioss::Field::COMPLEX:
      return std::string("complex");
    case Ioss::Field::STRING:
      return std::string("string");
    case Ioss::Field::CHARACTER:
      return std::string("char");
    case Ioss::Field::INVALID:
      return std::string("invalid");
    default:
      return std::string("internal error");
    }
  }
  
  void error_message(const Ioss::Field &field,
		     Ioss::Field::BasicType requested_type)
  {
    std::ostringstream errmsg;
    errmsg << "ERROR: For field named '" << field.get_name()
	   << "', code requested value of type '" << type_string(requested_type)
	   << "', but field type is '" << type_string(field.get_type())
	   << "'. Types must match\n";
    IOSS_ERROR(errmsg);
  }
}

Ioss::Field::Field() :
  name_(""), rawCount_(0), transCount_(0), size_(0), index_(0),
  type_(INVALID), role_(INTERNAL), rawStorage_(NULL), transStorage_(NULL)
{ rawStorage_ = transStorage_ = Ioss::VariableType::factory("invalid"); }

Ioss::Field::Field(const std::string &name,
		   const Ioss::Field::BasicType type,
		   const std::string &storage,
		   const Ioss::Field::RoleType role,
		   size_t value_count, size_t index) :
  name_(name), rawCount_(value_count), transCount_(value_count), size_(0), index_(index),
  type_(type), role_(role), rawStorage_(NULL), transStorage_(NULL)
{
  rawStorage_ = transStorage_ = Ioss::VariableType::factory(storage);
  size_ = internal_get_size(type_, rawCount_, rawStorage_);
}

Ioss::Field::Field(const std::string &name,
		   const Ioss::Field::BasicType type,
		   const std::string &storage,
		   int copies, 
		   const Ioss::Field::RoleType role,
		   size_t value_count, size_t index) :
  name_(name), rawCount_(value_count), transCount_(value_count), size_(0), index_(index),
  type_(type), role_(role), rawStorage_(NULL), transStorage_(NULL)
{
  rawStorage_ = transStorage_ = Ioss::VariableType::factory(storage, copies);
  size_ = internal_get_size(type_, rawCount_, rawStorage_);
}

Ioss::Field::Field(const std::string &name,
		   const Ioss::Field::BasicType type,
		   const Ioss::VariableType *storage,
		   const Ioss::Field::RoleType role,
		   size_t value_count, size_t index) :
  name_(name), rawCount_(value_count), transCount_(value_count), size_(0), index_(index),
  type_(type), role_(role), rawStorage_(storage), transStorage_(storage)
{
  size_ = internal_get_size(type_, rawCount_, rawStorage_);
}

Ioss::Field::Field(const Ioss::Field& from)
  : name_(from.name_), rawCount_(from.rawCount_), transCount_(from.transCount_),
    size_(from.size_), index_(from.index_),
    type_(from.type_), role_(from.role_),
    rawStorage_(from.rawStorage_), transStorage_(from.transStorage_), 
    transforms_(from.transforms_)
{}

Ioss::Field& Ioss::Field::operator=(const Field& from)
{
  name_ = from.name_;
  rawCount_ = from.rawCount_;
  transCount_ = from.transCount_;
  size_ = from.size_;
  index_ = from.index_;
  type_  = from.type_;
  role_ = from.role_;
  rawStorage_ = from.rawStorage_;
  transStorage_ = from.transStorage_;
  transforms_ = from.transforms_;
  return *this;
}

Ioss::Field::~Field()
{
}

// Verify that data_size is valid.
// If return value >= 0, then it is the maximum number of
// entities to get 'count'
// Throws runtime error if data_size too small.
size_t Ioss::Field::verify(size_t data_size) const
{
  if (data_size > 0) {
    // Check sufficient storage
    size_t required = get_size();
    if (required > data_size) {
	std::ostringstream errmsg;
	errmsg << "Field " << name_ << " requires "
	       << required  << " bytes to store its data. Only "
	       << data_size << " bytes were provided." << std::endl;
	IOSS_ERROR(errmsg);
    }
  }
  return rawCount_;
}

void Ioss::Field::check_type(BasicType the_type) const
{
  if (type_ != the_type)
    error_message(*this, the_type);
}

void Ioss::Field::reset_count(size_t new_count)
{
  if (transCount_ == rawCount_) transCount_ = new_count;
  rawCount_ = new_count;
  size_ = 0;
}

void Ioss::Field::reset_type(Ioss::Field::BasicType new_type)
{
  type_ = new_type;
  size_ = 0;
}

// Return number of bytes required to store entire field
size_t Ioss::Field::get_size() const
{
  if (size_ == 0) {
    Ioss::Field *new_this = const_cast<Ioss::Field*>(this);
    new_this->size_ = internal_get_size(type_, rawCount_, rawStorage_);

    new_this->transCount_   = rawCount_;
    new_this->transStorage_ = rawStorage_;
    std::vector<Transform*>::const_iterator I = transforms_.begin();
    while (I != transforms_.end()) {
      Transform* my_transform = *I++;
      new_this->transCount_   = my_transform->output_count(transCount_);
      new_this->transStorage_ = my_transform->output_storage(transStorage_);
      size_t size = internal_get_size(type_, transCount_, transStorage_);
      if (size > size_)
	new_this->size_ = size;
    }
  }
  return size_;
}

bool Ioss::Field::add_transform(Transform* my_transform)
{
  const Ioss::VariableType *new_storage = my_transform->output_storage(transStorage_);
  size_t new_count = my_transform->output_count(transCount_);

  if (new_storage != NULL && new_count > 0) {
    transStorage_ = new_storage;
    transCount_   = new_count;
  } else {
    return false;
  }

  if (transCount_ < rawCount_)
    role_ = REDUCTION;

  size_t size = internal_get_size(type_, transCount_, transStorage_);
  if (size > size_)
    size_ = size;

  transforms_.push_back(my_transform);
  return true;
}

bool Ioss::Field::transform(void *data)
{
  std::vector<Transform*>::const_iterator I = transforms_.begin();
  transStorage_ = rawStorage_;
  transCount_   = rawCount_;

  while (I != transforms_.end()) {
    Transform* my_transform = *I++;
    my_transform->execute(*this, data);

    transStorage_ = my_transform->output_storage(transStorage_);
    transCount_   = my_transform->output_count(transCount_);
  }
  return true;
}

namespace {
  size_t internal_get_size(Ioss::Field::BasicType type, size_t count,
			   const Ioss::VariableType* storage)
  {
    // Calculate size of the low-level data type
    size_t basic_size = 0;
    switch (type) {
    case Ioss::Field::REAL:
      basic_size = sizeof(double);
      break;
    case Ioss::Field::INTEGER:
      basic_size = sizeof(int);
      break;
    case Ioss::Field::COMPLEX:
      basic_size = sizeof(Complex);
      break;
    case Ioss::Field::STRING:
      basic_size = sizeof(std::string*);
      break;
    case Ioss::Field::CHARACTER:
      basic_size = sizeof(char);
      break;
    case Ioss::Field::INVALID:
      basic_size = 0;
      break;
    }
    // Calculate size of the storage type
    size_t storage_size = storage->component_count();

    return basic_size * storage_size * count;
  }
}
