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

#ifndef STATE_VECTOR_H_
#define STATE_VECTOR_H_

#include <fstream>
#include <vector>

#include <stdio.h>

#include "common.h"
#include "half_trie.h"

class CudaAllocator;

class StateVector {
	private:
		CudaAllocator &all_;

		unsigned int state_count_;
		size_t size_;
		mutable ST_BLOCK *vector_;
		mutable ST_BLOCK *device_;

	public:
		StateVector(const StateVector &state_vector);
		StateVector(unsigned int state_c, CudaAllocator &a);
		~StateVector();

		void set_bit(unsigned which);
		void set_bits(const StateVector &bits, unsigned base = 0);
		size_t get_size() const;
		size_t get_count() const;
		std::vector<StateVector *> split(unsigned int split_num) const;

		ST_BLOCK *get_device() const;
		
		ST_BLOCK *get_host(bool sync=0) const;
		
		ST_BLOCK *get_mutable_host();

		std::string toString() const;
		std::string toString(unsigned int split) const;
		
		void free_device();
		void free_host();
};

#endif


