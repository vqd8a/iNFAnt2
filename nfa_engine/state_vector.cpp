// Original work:
// Copyright (C) 2010  
// Pierluigi Rolando (pierluigi.rolando@polito.it)
// Netgroup - DAUIN - Politecnico di Torino
//
// Niccolo' Cascarano (niccolo.cascarano@polito.it)
// Netgroup - DAUIN - Politecnico di Torino
//
// Modified work:
// Copyright (C) 2017  
// Vinh Dang (vqd8a@virginia.edu)
// University of Virginia
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#include "cuda_allocator.h"
#include "state_vector.h"

#include <sstream>

using namespace std;

StateVector::StateVector(const StateVector &other) :
	all_(other.all_), state_count_(other.state_count_), size_(other.size_),
	device_(0) {
	
	vector_ = all_.alloc_host<ST_BLOCK>(size_);
	memmove(vector_, other.get_host(), size_);
}

StateVector::StateVector(unsigned int state_c, CudaAllocator &a) :
	all_(a), state_count_(state_c), vector_(0), device_(0)
{
	size_ = BYTE_ROUND_UP(state_c, bit_sizeof(ST_BLOCK))*sizeof(ST_BLOCK);

	vector_ = all_.alloc_host<ST_BLOCK>(size_);
	assert(vector_);
	memset(vector_, 0, size_);

	return;
}

StateVector::~StateVector() {
}

size_t StateVector::get_size() const {
	return size_;
}

size_t StateVector::get_count() const {
	return state_count_;
}

void StateVector::set_bit(unsigned which) {
	const unsigned int chunk = which/bit_sizeof(ST_BLOCK);
	assert(chunk < size_);
	const unsigned int off   = which%bit_sizeof(ST_BLOCK);
	vector_[chunk] |= 1 << off;
	return;
}

void StateVector::set_bits(const StateVector &bits, unsigned base)
{
	assert(bits.get_count() <= get_count() - base);
	memmove(vector_ + base / bit_sizeof(*vector_), bits.get_host(), 
			bits.get_size());

	return;
}


ST_BLOCK *StateVector::get_device() const {
	if(!device_)
		device_ = all_.alloc_device<ST_BLOCK>(get_size());
	assert(device_);

	cudaError_t retval;
	retval = cudaMemcpy(device_, vector_, size_, cudaMemcpyHostToDevice);

	CUDA_CHECK(retval, "Error while copying state vector to device memory");

	return device_;
}

ST_BLOCK *StateVector::get_host(bool sync) const {
	if(sync) {
		cudaError_t retval;
		retval = cudaMemcpy(vector_, device_, size_, cudaMemcpyDeviceToHost);
		CUDA_CHECK(retval, "Error while synchronizing host state vector");
	}	
	return vector_;
}

ST_BLOCK *StateVector::get_mutable_host() {
	return vector_;
}

string StateVector::toString() const {
	stringstream ss;

	ss << "{" << size_ << "} ";
	for (unsigned i = 0; i < size_ * bit_sizeof(unsigned char); i++) {
		ST_BLOCK mask = 1;
		mask <<= i % bit_sizeof(ST_BLOCK);
		if (vector_[i/bit_sizeof(ST_BLOCK)] & mask)
			ss << i << " ";
	}

	return ss.str();
}

string StateVector::toString(unsigned int split_num) const {
	stringstream ss;
	const unsigned int sub_vector_count = get_count() / split_num;

	vector<StateVector *> final_vectors = split(split_num);
	assert(final_vectors.size() == split_num);

	for (unsigned cur_v = 0; cur_v < split_num; ++cur_v) {
		ss << "Vector " << cur_v << ": ";
		ss << final_vectors[cur_v]->toString() << endl;
		delete final_vectors[cur_v];
	}

	return ss.str();
}

vector<StateVector *> StateVector::split(unsigned int split_num) const {
	vector<StateVector *>final_vectors;
	const unsigned int sub_vector_size = get_size() / split_num;
	const unsigned int sub_vector_count = get_count() / split_num;

#ifdef DEBUG_1
	cout << "Long vector: " << toString() << endl;
#endif
	for (unsigned i = 0; i < split_num; ++i) {
		StateVector *sv = new StateVector(sub_vector_count, all_);
#ifdef DEBUG_1
		cout << "Memmove from " << i * sub_vector_size << " of: " 
			<< sv->get_size()<< endl;
#endif
		memmove(sv->get_mutable_host(), (char *)get_host() + i * sub_vector_size,
				sv->get_size());
		final_vectors.push_back(sv);

#ifdef DEBUG_1
		cout << "Splitted: " << i << " " <<  sv->toString() << endl;
#endif
	}

	return final_vectors;
}

void StateVector::free_device() {
	all_.dealloc_device(device_);
}

void StateVector::free_host() {
	all_.dealloc_host(vector_);
}