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

#include <boost/algorithm/string/regex.hpp>
#include <boost/foreach.hpp>
#include <boost/iterator.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/multi_array.hpp>
#include <boost/range.hpp>
#include <boost/regex.hpp>
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <pthread.h>//parallelize transition duplicating section
#include <time.h>
#include <sys/time.h>

#include "config_options.h"
#include "cuda_allocator.h"
#include "globals.h"
#include "half_trie.h"
#include "transition_graph.h"

#include <algorithm>//for "find" function

using namespace std;
using namespace half_trie;

extern ConfigOptions cfg;

int print_en_transition_graph=0;//debug purpose: 0-don't print; 1-print
/*------------------------------------------------------------------------------------*/
struct MyComparator {
    const vector<unsigned short> & value_vector;

    MyComparator(const vector<unsigned short> & val_vec): value_vector(val_vec) {}

    bool operator()(int i1, int i2)
    {
        return value_vector[i1] < value_vector[i2];
    }
};
/*------------------------------------------------------------------------------------*/
void set_bit_(unsigned int which, ST_BLOCK *vector) {
	unsigned int chunk = which/bit_sizeof(ST_BLOCK);
	unsigned int off   = which%bit_sizeof(ST_BLOCK);
	vector[chunk] |= 1 << off;
}
/*------------------------------------------------------------------------------------*/
struct thread_data{
		std::vector< std::vector<unsigned int> > dst_state_groups;
		std::vector< std::vector<st_t> > src_table_vec;
		std::vector< std::vector<st_t> > nfa_table_vec;
		int symstart;
        int Symspercpu;
     	st_t *nfa_table;
		st_t *src_table;
		unsigned int *offset_table;
        int trans_cnt;
		int threadID;
};
/*------------------------------------------------------------------------------------*/
void *SearchAndDup (void *p) {
	struct thread_data *my_data;  
	vector<st_t> src_state_ids;
	vector<st_t> dst_state_ids;
	
	my_data = (struct thread_data *) p;
	
	unsigned int trans_cnt_=0;
	for (unsigned short j = 0; j < my_data->Symspercpu; j++) {
		unsigned base = my_data->offset_table[j   + my_data->symstart];
		unsigned end  = my_data->offset_table[j+1 + my_data->symstart];
		cout << "Thread " << my_data->threadID <<", Duplicating: symbol "<< j << endl;
		for (unsigned short prev_j = 0; prev_j < cfg.get_alphabet_size(); prev_j++) {
			unsigned int counters_=0;
			for (unsigned jj = base; jj < end; jj++){
				if (my_data->src_table[jj] != 0xFFFF) {
					std::vector<unsigned int>::iterator it;
					it = std::find(my_data->dst_state_groups[prev_j].begin(), my_data->dst_state_groups[prev_j].end(), my_data->src_table[jj]);
					if ( it != my_data->dst_state_groups[prev_j].end() ){
						src_state_ids.push_back(my_data->src_table[jj]);
						dst_state_ids.push_back(my_data->nfa_table[jj]);
						counters_++;
						trans_cnt_++;
					}
				}
			}
			if(counters_ % 2) {
				src_state_ids.push_back(0xFFFF);//???necessary???
				dst_state_ids.push_back(0xFFFF);//???necessary???
				trans_cnt_++;
			}
			my_data->src_table_vec.push_back(src_state_ids);
			my_data->nfa_table_vec.push_back(dst_state_ids);
			src_state_ids.clear();
			dst_state_ids.clear();
		}
	}
	my_data->trans_cnt = trans_cnt_;
	pthread_exit(NULL);
	return NULL;
}
/*------------------------------------------------------------------------------------*/
TransitionGraph::TransitionGraph(istream &file, CudaAllocator &allocator, unsigned int gid, unsigned int *trans_per_sym) ://Changed
	all_(allocator), persistent_states_(cfg.get_state_count(gid), all_), initial_states_(cfg.get_state_count(gid), all_), accept_states_(cfg.get_state_count(gid), all_)
{
	using boost::algorithm::split_regex;
	using boost::regex;
	using boost::lexical_cast;

	transition_count_ = 0;

	Graph tnsn_graph;
	string line;
    
	unsigned int persistent_count = 0;
	
	while(getline(file, line)) {
		list<string> parts;

		// Skip comments
		if (line[0] == '#')
			continue;

		split_regex(parts, line, regex("[[:blank:]]*:[[:blank:]]*"));

		// Skip empty lines
		if(parts.size() == 0)
			continue;

		vector<string> srcdst;
		split_regex(srcdst, parts.front(), regex("[[:blank:]]*->[[:blank:]]*"));

		parts.pop_front();

		// Skip empty atoms
		if(parts.size() == 0)
			continue;

		assert(srcdst.size() >= 1);
		assert(parts.front().size());

		unsigned src = lexical_cast<unsigned>(srcdst.front());

		if(srcdst.size() == 2) {
			/* Handle a transition */
			unsigned dst = lexical_cast<unsigned>(srcdst.back());
			
			vector<string> atoms;
			split_regex(atoms, parts.front(), regex("[[:blank:]]+(?!\\|)[[:blank:]]*"));
            
			//cout << "      MY TEST: src= " << src << " dst= " << dst << ", atoms.size()= " << atoms.size() << endl;
			//cout << "      ";
			//for(int i = 0; i < atoms.size(); ++i) cout << atoms[i] << " ";
			//cout << endl;
	
			BOOST_FOREACH(string s, atoms) {
				if(s.size() > 0) {
					vector<string> ssrcend;
					split_regex(ssrcend, s, regex("\\|"));
					assert(ssrcend.size() <= 2);

					unsigned srange, erange;
					srange = erange = lexical_cast<unsigned>(ssrcend.front());
					if(ssrcend.size() > 1) {
						erange = lexical_cast<unsigned>(ssrcend.back());
					}
					assert(erange < cfg.get_alphabet_size()); //cout << "      srange= "<< srange << ", erange= " << erange << ", cfg.get_alphabet_size()= " << cfg.get_alphabet_size() << endl;
					
					if ((src==dst)&&(srange==0)&&(erange==255)) {
						//persistent_states_.set_bit(src);//test: states with a self-transition on every character of the alphabet
						persistent_count ++;
					}
					//else {
						for(unsigned i = srange; i <= erange; ++i)
							add(tnsn_graph, src, dst, Ranges(Range(srange, erange)));
					//}
				}
			}
		} else if(parts.front().find("accepting") != string::npos) {
			/* Handle an accepting state and its related rules */
			vector<string> atoms;
			split_regex(atoms, parts.front(), regex("[[:blank:]]+"));

			BOOST_FOREACH(string s, atoms) {
				if(s.size() > 0 && isdigit(s[0])) {
					accepting_states_[src].insert(lexical_cast<unsigned>(s)); //insert a new value of rule to the set of this state
					//cout << "ACCE_states_ " << src << ", Rule: " << s << endl;
				}
			}
			//accepting_states_[23].insert(100);//TEST
			
            //printf("persistent_states_[%d]\n",src);//Check content of "persistent_states_" vector?? It seems to me that this vector is the same as "accepting_states_" vector
			
			//persistent_states_.set_bit(src);//Commented
			
			//Commenting the line aboved seems working fine. Matchings are kept track by other data structures: match_count, match_offset, match_states
			//The purpose of commenting is for the simplicity of coding in the optimization part: searching for compatible groups; 
			//and for having the same output as vasim: in some NFAs, accepting states are not the last states, but they transition to other states, so keeping these states as persistent states is not a good idea 
			
			accept_states_.set_bit(src);
		} else if(parts.front().find("initial") != string::npos) {
			/* Enable the initial state in all the required state vector */
			initial_states_.set_bit(src);
		}  else {
			cerr << "[warning] cannot parse line '" << line << "'\n";
		}
	}

	//clog << "Read " << transition_count_ << " transitions.\n";
	//clog << "Have " << persistent_count << " persistent states.\n";

	/* Allocate the iNFAnt NFA data structure in host memory and fill it */
#if DEBUG
	cout << "Alphabet size: " << cfg.get_alphabet_size() << endl;
#endif
	
	offset_table_size_ = (cfg.get_alphabet_size()+1)*sizeof(*offset_table_);
	
	// padding is not required as in the original iNFAnt 
	nfa_table_size_ = transition_count_*sizeof(*nfa_table_); 
	nfa_table_      = all_.alloc_host<st_t>(nfa_table_size_);
	src_table_      = all_.alloc_host<st_t>(nfa_table_size_);
	offset_table_   = all_.alloc_host<unsigned int>(offset_table_size_);
	
	//cout<< "sizeof(*offset_table_) = " << sizeof(*offset_table_) << " and sizeof(*nfa_table_) = " << sizeof(*nfa_table_) << endl;
	//cout << " sizeof(*nfa_table_) = " << sizeof(*nfa_table_) << ", nfa_table_ size (bytes)= " << nfa_table_size_ << endl;
	
	hash_map<symbol_t, map<Graph::state_t, set<Graph::state_t> > > hm;
	tnsn_graph.project(hm);

	typedef boost::multi_array<unsigned, 2> array_t;
	typedef array_t::index index_t;
	array_t tnsn_per_char(
			boost::extents[cfg.get_alphabet_size()][cfg.get_state_count(gid)]);
	for(index_t i = 0; i < cfg.get_alphabet_size(); ++i)
		for(index_t j = 0; j < cfg.get_state_count(gid); ++j)
			tnsn_per_char[i][j] = 0;

	unsigned nfa_current = 0;
	for (symbol_t s = 0; s < cfg.get_alphabet_size(); ++s) {
		offset_table_[s] = nfa_current;
		
		hash_map<symbol_t, map<Graph::state_t, set<Graph::state_t> > >::iterator ss;
		if((ss = hm.find(s)) == hm.end())
			continue;

		map<Graph::state_t, set<Graph::state_t> >::iterator ii;
		for(ii = ss->second.begin(); ii != ss->second.end(); ++ii) {
			set<Graph::state_t>::iterator jj;
			for(jj = ii->second.begin(); jj != ii->second.end(); ++jj) {
				tnsn_per_char[ss->first][ii->first]++;
				nfa_table_[nfa_current] = *jj;
				src_table_[nfa_current++] = ii->first;
			}
		}
	}
	offset_table_[cfg.get_alphabet_size()] = nfa_current;

	//iNFAnt2: Represent accepting states in the NFA state table as negative numbers 
	for (unsigned int i = 0; i < transition_count_; i++) {
		unsigned int state_bit   =  1 << (((unsigned int)nfa_table_[i]) % bit_sizeof(ST_BLOCK));
		unsigned int state_chunk =        ((unsigned int)nfa_table_[i]) / bit_sizeof(ST_BLOCK);
		ST_BLOCK *tmp_st_vector = accept_states_.get_host(false);
		ST_BLOCK match_check = state_bit & tmp_st_vector[state_chunk];
		if(match_check)
			nfa_table_[i] = -nfa_table_[i];	
	}
	//
	//cout << "ORIGINAL iNFAnt CODE" << endl;
	//Keep track of the number of transitions of each symbol
	unsigned int cnt_=0;
	for (unsigned short j = 0; j < cfg.get_alphabet_size(); j++) {
		trans_per_sym[j]=offset_table_[j+1] - offset_table_[j];
		cnt_ += offset_table_[j+1] - offset_table_[j];
	}
	cout << "Total transitions (including paddings): " << cnt_<< endl;
	
	//clog << "NFA loading done.\n";
	return;
}
/*------------------------------------------------------------------------------------*/
void TransitionGraph::copy_to_device(){
	cudaError_t retval;
			
	d_nfa_table_ = all_.alloc_device<ST_BLOCK>(nfa_table_size_);
	retval = cudaMemcpy(d_nfa_table_, nfa_table_, nfa_table_size_,
			cudaMemcpyHostToDevice);
	CUDA_CHECK(retval, "Error while copying NFA table to device memory");

	d_src_table_ = all_.alloc_device<ST_BLOCK>(nfa_table_size_);
	retval = cudaMemcpy(d_src_table_, src_table_, nfa_table_size_,
			cudaMemcpyHostToDevice);
	CUDA_CHECK(retval, "Error while copying NFA src table to device memory");

	d_offset_table_ = all_.alloc_device<unsigned int>(offset_table_size_);
	retval = cudaMemcpy(d_offset_table_, offset_table_,
			offset_table_size_, cudaMemcpyHostToDevice);
	CUDA_CHECK(retval,
			"Error while copying character start offset table to device memory");

#ifdef DEBUG
	cout << "Automata:" << endl;
	for (unsigned short i = 0; i < cfg.get_alphabet_size(); i++) {
		
		unsigned base = offset_table_[i];
		unsigned end = offset_table_[i+1];
		//cout << "Symbol " << i << ", base " << base << ", end " << end << endl;
		printf("Symbol: %d, base: %d, end: %d, number of transition (end-base): %d\n", i, base, end, end-base);
		/*for (unsigned j = base; j < end; ++j) {
			//cout << src_table_[j] << " -> " << nfa_table_[j] << endl;
			printf("%d -> %d\n",src_table_[j],nfa_table_[j]);
		}*/
	}

#endif
#ifdef DEBUG
	cout << "Automata by night" << endl;
	for (unsigned i = 0; i < nfa_table_size_/2; i++) {
		cout << src_table_[i] << " -> " << nfa_table_[i] << endl;
	}
#endif

	all_.dealloc_host(nfa_table_);
	all_.dealloc_host(src_table_);
	all_.dealloc_host(offset_table_);
	
	return;
}
/*------------------------------------------------------------------------------------*/
namespace boost {
	template<> struct range_const_iterator<Ranges> {
		typedef Ranges::const_iterator type;
	};
}
/*------------------------------------------------------------------------------------*/
void TransitionGraph::add(Graph &tnsn_graph, unsigned src, unsigned dst,
		const Ranges r)
{
	for (Ranges::const_iterator ii = r.begin(); ii != r.end(); ++ii) {
		if(!tnsn_graph.get(src, dst).includes(*ii))
			++transition_count_;
	}

	tnsn_graph.add(src, dst, r);
	return;
}
/*------------------------------------------------------------------------------------*/
void TransitionGraph::accepting_rules(const StateVector &final_vector, std::set<unsigned int> &results) const {
	
	vector<unsigned>accepted_rules;
	//cout << "size of final_vector " << final_vector.get_size() * bit_sizeof(unsigned char) << endl;
	for (unsigned i = 0; i < final_vector.get_size() * bit_sizeof(unsigned char); i++) {
		ST_BLOCK mask = 1;
		mask <<= i % bit_sizeof(ST_BLOCK);
		if (final_vector.get_host()[i/bit_sizeof(ST_BLOCK)] & mask) {
			map<unsigned, set<unsigned> >::const_iterator it = accepting_states_.find(i);
			if (it != accepting_states_.end()) {
				set<unsigned>::iterator iitt;
				for (iitt = it->second.begin();	iitt != it->second.end(); ++iitt)
					results.insert(*iitt);
			}
		}
	}
#ifdef DEBUG	
	//TEST ONLY
	for (unsigned i = 0; i < final_vector.get_size() * bit_sizeof(unsigned char); i++) {
		map<unsigned, set<unsigned> >::const_iterator it = accepting_states_.find(i);
		if (it != accepting_states_.end()) {
			set<unsigned>::iterator iitt;
			cout << "State: " << i << ", Rules: ";
			for (iitt = it->second.begin();	iitt != it->second.end(); ++iitt) cout << *iitt << " ";
			cout << endl;
		}		
	}	
	//END TEST
#endif
}
/*------------------------------------------------------------------------------------*/
//Added for matching operation
/*void TransitionGraph::mapping_states2rules(unsigned int *match_count, unsigned int *match_offset, unsigned int *match_states, 
                                           unsigned int match_vec_size, std::vector<unsigned long> cur_size_vec, std::ofstream &fp, int *rulestartvec, unsigned int gid) const {
	unsigned int total_matches=0;	
	for (int j = 0; j < cur_size_vec.size(); j++)	total_matches += match_count[j];
	fp   << "REPORTS: Total matches: " << total_matches << endl;
	
	for (int j = 0; j < cur_size_vec.size(); j++) {
		for (unsigned i = 0; i < match_count[j]; i++) {
			map<unsigned, set<unsigned> >::const_iterator it = accepting_states_.find(match_states[match_vec_size*j + i]);		
			if (j==0) fp   << match_offset[match_vec_size*j + i] << "::" << endl;
			else      fp   << match_offset[match_vec_size*j + i] + cur_size_vec[j-1] << "::" << endl;
			if (it != accepting_states_.end()) {
				set<unsigned>::iterator iitt;
				for (iitt = it->second.begin();	iitt != it->second.end(); ++iitt) {
					fp   << "    Rule: " << *iitt + rulestartvec[gid] << endl;
				}
			}
		}
	}
}*/
void TransitionGraph::mapping_states2rules(unsigned int *match_count, match_type *match_array, unsigned int match_vec_size, std::vector<unsigned long> cur_size_vec, std::vector<unsigned long> pad_size_vec, std::ofstream &fp
#ifdef DEBUG
                                                                         , int *rulestartvec, unsigned int gid
#endif
                                                                                                             ) const {//version 2: multi-byte fetching
	unsigned int total_matches=0;	
	for (int j = 0; j < cur_size_vec.size(); j++)	total_matches += match_count[j];
	fp   << "REPORTS: Total matches: " << total_matches << endl;
	
	for (int j = 0; j < cur_size_vec.size(); j++) {
		for (unsigned i = 0; i < match_count[j]; i++) {
			map<unsigned, set<unsigned> >::const_iterator it = accepting_states_.find(match_array[match_vec_size*j + i].stat);		
			if (j==0) fp   << match_array[match_vec_size*j + i].off  << "::" << endl;
			else      fp   << match_array[match_vec_size*j + i].off + cur_size_vec[j-1] - (pad_size_vec.empty() ? 0 : pad_size_vec[j-1]) << "::" << endl;
			if (it != accepting_states_.end()) {
				set<unsigned>::iterator iitt;
				for (iitt = it->second.begin();	iitt != it->second.end(); ++iitt) {
#ifdef DEBUG
					fp   << "    Rule: " << *iitt + rulestartvec[gid] << endl;
#else
                    fp   << "    Rule: " << *iitt << endl;
#endif
				}
			}
		}
	}
}
//END
/*------------------------------------------------------------------------------------*/
TransitionGraph *load_nfa_file(const char *pattern_name, unsigned int gid, unsigned int *trans_per_sym) {//Changed
	int trn_good = 0;
	vector<ifstream *> trn_file;
	ifstream file;
	
	// Setting for no multistride
	cfg.set_max_out(256);

	/* Load the translation tries from disk */
	if (pattern_name) {
		vector<string> prefixes;
		prefixes.push_back("/2");
		prefixes.push_back("/4");
		prefixes.push_back("/8");
		prefixes.push_back("/F");

		pair<string, string> p = utils::dir_base_name(pattern_name);
		
		//cout << "cfg.get_good_tries() = " << cfg.get_good_tries() << endl;
		
		string dump_filename;
		if (cfg.get_good_tries()) {
			dump_filename = p.first + prefixes[cfg.get_good_tries() - 1] + p.second 
				+ ".nfa";
		}
		else {
			dump_filename = p.first + "/" + p.second + ".nfa";//orig
			//dump_filename = p.first + "/" + p.second + ".dumpdfa";//testing
			//cout << "pattern_name = " << pattern_name << ", p.first = " << p.first << ", p.second = " << p.second << ", dump_filename = " << dump_filename << endl;
		}

		file.open(dump_filename.c_str());
		if (!file.good()) {
			cout << "Impossibile aprire il file " << dump_filename << endl;
			return NULL;
		}

		for(unsigned int i = 0; i < prefixes.size() &&	i < cfg.get_good_tries(); ++i) {
			string cur_name = p.first + prefixes[i] + p.second + ".translation";
			trn_file.push_back(new ifstream(cur_name.c_str()));
			cfg.get_mutable_max_outs()[i+1] = 0;
			if(trn_file[i]->good()) {
				cfg.add_trns(new HalfTrie<symbol_t, out_symbol_t>(*trn_file[i],
							cfg.get_mutable_max_outs()[i+1]));
				++trn_good;
			} else
				clog << "No good " << i << '\n';
		}	
		
		if (cfg.get_good_tries() > 0) {
#ifdef DEBUG
			cout << "Good tries: " << trn_good << endl;
#endif
			assert(trn_good);
			cfg.set_good_tries(trn_good);
			cfg.set_max_out(cfg.get_mutable_max_outs()[trn_good]);
		}
	} else {
		cout << "Ruleset name is invalid" << endl;
		return NULL;
	}
	
	/* Read the state count from the first transition graph line */
	{
		string line("#");
		
		unsigned int t;
		
		while (line[0] == '#') {
			if(!getline(file, line))
				error("reading NFA size", NFA_ERROR);

			istringstream iss(line);
			iss >> t;
		}
		cfg.set_state_count(t);
	}
		
	/* Read the transition graph from disk */
	TransitionGraph *tg = new TransitionGraph(file, cfg.get_allocator(), gid, trans_per_sym);
	
	/* Create the state vectors */
	tg->tmp_states_ = new StateVector(tg->get_mutable_initial_states());//this might cause memory overwritten, conflicts the resizing of memory footprint of group_offset_optim
	cfg.set_state_vector(tg->tmp_states_);
	
	return tg;
}
/*------------------------------------------------------------------------------------*/
ST_BLOCK *TransitionGraph::get_d_nfa_table() const {
	return d_nfa_table_;
}

ST_BLOCK *TransitionGraph::get_d_src_table() const {
	return d_src_table_;
}

unsigned int *TransitionGraph::get_d_offset_table() const {
	return d_offset_table_;
}

st_t *TransitionGraph::get_nfa_table() {
	return nfa_table_;
}

st_t *TransitionGraph::get_src_table() {
	return src_table_;
}

unsigned int *TransitionGraph::get_offset_table() {
	return offset_table_;
}

StateVector &TransitionGraph::get_mutable_persistent_states() {
	return persistent_states_;
}

StateVector &TransitionGraph::get_mutable_initial_states() {
	return initial_states_;
}

StateVector &TransitionGraph::get_accept_states() {
	return accept_states_;
}

std::map<unsigned, std::set<unsigned> > &TransitionGraph::get_accepting_states_() {
//Returning a reference means that the calling code can modify the value of your member variable after you return. 
	return accepting_states_;
}

size_t TransitionGraph::get_nfa_table_size() const {
	return nfa_table_size_;
}

size_t TransitionGraph::get_offset_table_size() const {
	return offset_table_size_;
}

void TransitionGraph::free_devmem(){
	all_.dealloc_device(d_nfa_table_);
	all_.dealloc_device(d_src_table_);
	all_.dealloc_device(d_offset_table_);	
}

unsigned int TransitionGraph::get_transition_count() const {
	return transition_count_;
}//For testing only

void TransitionGraph::set_nfa_table_size(size_t new_size){
	nfa_table_size_ = new_size;
}

void TransitionGraph::set_offset_table_size(size_t new_size){
	offset_table_size_ = new_size;	
}

unsigned int *TransitionGraph::get_group_offset() {
	return group_offset_;
}

st_t *TransitionGraph::get_nfa_table_optim() {
	return nfa_table_optim_;
}

st_t *TransitionGraph::get_src_table_optim() {
	return src_table_optim_;
}

void TransitionGraph::free_hostmem(){
	all_.dealloc_host(nfa_table_);
	all_.dealloc_host(src_table_);
	all_.dealloc_host(offset_table_);
}
