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

#ifndef TRANSITION_GRAPH_H_
#define TRANSITION_GRAPH_H_

#include <fstream>
#include <map>
#include <set>
#include <vector>

#include <stdio.h>

#include "common.h"
#include "half_trie.h"
#include "sparse.h"
#include "state_vector.h"
#include "transition_graph.h"

class TransitionGraph {
	private:
		CudaAllocator &all_;

		size_t nfa_table_size_;
		size_t offset_table_size_;
		size_t state_vector_size_;

		st_t *nfa_table_; // transition sequence
		st_t *src_table_; // transition source table
		unsigned int *offset_table_;
		
		st_t *nfa_table_optim_; // transition sequence
		st_t *src_table_optim_; // transition source table
        unsigned int *group_offset_;
		///////
		
		ST_BLOCK *d_nfa_table_; // transition sequence
		ST_BLOCK *d_src_table_; // transition source table
		unsigned int *d_offset_table_;
		StateVector persistent_states_;
		StateVector initial_states_;
        StateVector accept_states_;
		
		unsigned int transition_count_;

		std::map<unsigned, std::set<unsigned> > accepting_states_;

		void add(Graph & tnsn_graph_, unsigned src, unsigned dst, const Ranges r);

	public:
		TransitionGraph(std::istream &, CudaAllocator &, unsigned int, unsigned int *trans_per_sym);//Changed

		void copy_to_device();
		void accepting_rules(const StateVector &, 
				std::set<unsigned int> &results) const;
		
		//Added for matching operation		
		//void mapping_states2rules(unsigned int *match_count, unsigned int *match_offset, unsigned int *match_states, 
		//                          unsigned int match_vec_size, std::vector<unsigned long> cur_size_vec, std::ofstream &fp, int *rulestartvec, unsigned int gid) const;
		void mapping_states2rules(unsigned int *match_count, match_type *match_array, unsigned int match_vec_size, std::vector<unsigned long> cur_size_vec, std::vector<unsigned long> pad_size_vec, std::ofstream &fp
#ifdef DEBUG
		                                                           , int *rulestartvec, unsigned int gid
#endif
		                                                                                                 ) const;//multi-byte fetching
								  
		ST_BLOCK *get_d_nfa_table() const;
		ST_BLOCK *get_d_src_table() const;
		unsigned int *get_d_offset_table() const;

		st_t *get_nfa_table();
		st_t *get_src_table();
		unsigned int *get_offset_table();
				
		StateVector &get_mutable_persistent_states();
		StateVector &get_mutable_initial_states();
		StateVector &get_accept_states();
		
		std::map<unsigned, std::set<unsigned> > &get_accepting_states_();
		//Returning a reference means that the calling code can modify the value of your member variable after you return. 

		size_t get_nfa_table_size() const;
		size_t get_offset_table_size() const;
		
		void free_devmem();
		StateVector *tmp_states_;
		unsigned int get_transition_count() const;//For testing only
		
		void set_nfa_table_size(size_t new_size);
		void set_offset_table_size(size_t new_size);
		unsigned int *get_group_offset();
		st_t *get_nfa_table_optim();
		st_t *get_src_table_optim();
		void free_hostmem();
};

TransitionGraph *load_nfa_file(const char *pattern_name, unsigned int gid, unsigned int *trans_per_sym);

#endif /* HOST_FUNCTIONS_H_ */
