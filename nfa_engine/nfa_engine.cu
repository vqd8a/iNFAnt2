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

#include <iostream>
#include <fstream>
#include <string>

#include <stdio.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <algorithm>

#include "burst.h"
#include "globals.h"
#include "host_functions.h"

using namespace std;

size_t getFilesize(const char* filename);
void Usage(void);
bool ParseCommandLine(int argc, char *argv[]);

const char *pattern_name = NULL;

#ifdef DEBUG
const char *timing_filename = NULL;
const char *blksiz_filename = NULL;
#endif

int total_rules   = 0;
int blksiz_tuning = 0;

extern ConfigOptions cfg;

int main(int argc, char* argv[]){

	int retval;
	std::vector<TransitionGraph *> nfa_vec;

	char char_temp;
    char filename[1500], bufftmp[10];

	struct timeval c1, c2, c3, c4, c5, c6;
	long seconds, useconds;
	double t_alloc, t_kernel, t_collect, t_NFAload, t_in;	
	string str_line;

#ifdef DEBUG
	double t_exec, t_out;
	int rulespergroup, *rulestartvec;
	//ifstream fp_tmp;
	ofstream fp_timing;
	ofstream fp_trans;
	ofstream fp_blksiz;
#endif
	
	unsigned int *trans_per_sym;
	int blockSize;

	if (argc == 1) {
		Usage();
		return 0;
	}

	// Load the NFAs
	gettimeofday(&c1, NULL);
	
	retval = ParseCommandLine(argc, argv);
    
	if(!retval)
		return 0;

    unsigned int total_bytes = getFilesize(cfg.get_trace_file_name());

	cout<< "-----------------User input info--------------------" << endl;
	cout<< "Total number of rules (subgraphs): " << total_rules << endl;
	cout<< "Total input bytes: "   << total_bytes << endl;
	unsigned int n_subsets   = cfg.get_blocks_y();         
	unsigned int n_packets   = cfg.get_parallel_packets();
	unsigned int packet_size = ((total_bytes%n_packets)==0)?(total_bytes/n_packets):(total_bytes/n_packets+1);
	cout<< "Graph(s) (NFA(s)) combined: "   << n_subsets << endl;
	cout<< "Packet(s): "   << n_packets << endl;
	cout<< "Packet size (bytes): " << packet_size << endl;

	trans_per_sym    = (unsigned int*)malloc (n_subsets * (256+1) * sizeof(unsigned int));

#ifdef DEBUG
	rulestartvec = (int*)malloc (n_subsets * sizeof(int));
	
    if ((total_rules%n_subsets)==0)
		rulespergroup = total_rules/n_subsets;
	else
		rulespergroup = total_rules/n_subsets + 1;
	//printf("rulespergroup=%d\n",rulespergroup);
	for (unsigned int i=0; i<n_subsets; i++) {
		rulestartvec[i]=i*rulespergroup;
		//printf("rulestartvec[%d]=%d\n",i,rulestartvec[i]);
	}
#endif

	cout << endl;
	cout<< "-----------------Loading NFA(s) from text file(s)--------------------" << endl;
	cout << "Loading..." << endl;
	
	for (unsigned int i = 0; i < n_subsets; ++i) {		
		
		strcpy (filename,pattern_name);
		
		strcat (filename,"_");
		snprintf(bufftmp, sizeof(bufftmp),"%d",n_subsets);
		strcat (filename,bufftmp);
		strcat (filename,"/");
		snprintf(bufftmp, sizeof(bufftmp),"%d",i+1);
		strcat (filename,bufftmp); //cout<< "NFA " << i + 1 << ":"<< filename << endl;		
		
		TransitionGraph *nfa_tmp = NULL;
		
		nfa_tmp = load_nfa_file(filename, i, &trans_per_sym[i*(256+1)]);
				
		if(!nfa_tmp){
			printf("Error while loading NFA on the device\n");
			return 0;
		}
		nfa_vec.push_back(nfa_tmp);
	}
	
	cout << "\nNFA loading done!!!\n\n";
	
#ifdef DEBUG	
	for (unsigned int i = 0; i < n_subsets; ++i) {
		if (i!=n_subsets-1) cout << "Sub-ruleset "<< i + 1 << ": Rules: " << rulestartvec[i+1] - rulestartvec[i] <<", States: "<< cfg.get_state_count(i) 
		                         << ", Transitions: " << nfa_vec[i]->get_transition_count() << endl;	
	    else cout << "Sub-ruleset "<< i + 1 << ": Rules: " << total_rules - rulestartvec[i] <<", States: "<< cfg.get_state_count(i) 
		          << ", Transitions: " << nfa_vec[i]->get_transition_count() << endl;	
	}
	
	//Keep track of the number of transitions of each symbol
	strcpy (filename,"Trans_per_sym_");
	snprintf(bufftmp, sizeof(bufftmp),"%d",n_subsets);
	strcat (filename,bufftmp);
	strcat (filename,".bin");	
	fp_trans.open(filename,ios::binary | ios::out);
	fp_trans.write((char *)trans_per_sym, n_subsets * (256+1) * sizeof(unsigned int));
	fp_trans.close();	
#endif

	gettimeofday(&c2, NULL);
		
	printf("-----------------Starting nfa execution--------------------\n");
    	
	// open input stream file    
	ifstream fp(cfg.get_trace_file_name());
#ifdef DEBUG
	if (timing_filename != NULL)//and timing file
		fp_timing.open(timing_filename,ios::binary | ios::out);
	if (blksiz_filename != NULL)//and timing file
		fp_blksiz.open(blksiz_filename,ios::binary | ios::out);
#endif
	
	unsigned int processed_packets = 0;
{	   
    gettimeofday(&c3, NULL);
	
	Burst burst;
	vector<unsigned char> payload;
	vector<set<unsigned> > accepted_rules;
    vector<unsigned int> payload_count;

	burst.set_required_translations(cfg.get_good_tries());
	
	// Read input stream file
	unsigned int cnt2=0;
	unsigned int cnt=0;
	if (fp){								
		while ( fp.get(char_temp) ){
			cnt2++;	
			payload.push_back(char_temp);
			cnt++;					
			if (cnt==packet_size){
				//note: padding to each packet if packet_size is not evenly divided by fetch_bytes (e.g. 4, 8)
				if ( (cnt%fetch_bytes) != 0 ) {
					for (unsigned int i = 0; i < (fetch_bytes-(cnt%fetch_bytes)); i++)
						payload.push_back(0);
					burst.save_n_padded_bytes(fetch_bytes-(cnt%fetch_bytes));
				}
				burst.append_payload(payload);
				payload_count.push_back(payload.size());//cout << payload.size() << endl;
				processed_packets++;
				payload.clear();
				cnt=0;
			}
		}
		if ((cnt>0)&&(cnt<packet_size)){
			//-- note: padding to each packet if packet_size is not evenly divided by fetch_bytes (e.g. 4, 8)
			if ( (cnt%fetch_bytes) != 0 ) {
				for (unsigned int i = 0; i < (fetch_bytes-(cnt%fetch_bytes)); i++)
					payload.push_back(0);
				burst.save_n_padded_bytes(fetch_bytes-(cnt%fetch_bytes));
			}
			burst.append_payload(payload);
			payload_count.push_back(payload.size());//cout << payload.size() << endl;
			processed_packets++;
			payload.clear();
			cnt=0;
		}
	}
	else{
		cout<< "Cannot open input file" << endl;				
	}
	cout << "Number of processed packets: "<< processed_packets << " and total number of bytes: "<< cnt2 << endl;
	for (unsigned int i = 0; i < processed_packets; i++){
		cout << "Packet "<< i+1 << ": " << payload_count[i] << endl;				
	}
	
	for (unsigned int i = 0; i < n_subsets; i++) {//Changed
		burst.init_state_vector(cfg.get_state_vector(i), i);
	}
	
	gettimeofday(&c4, NULL);

	accepted_rules = nfa_execute(nfa_vec, burst, n_subsets, 
#ifdef DEBUG
	                             rulestartvec,
#endif	
	                             &t_alloc, &t_kernel, &t_collect, &blockSize, trans_per_sym, blksiz_tuning);

	gettimeofday(&c5, NULL);
		
	burst.free_device();
    
#ifdef DEBUG		
    for (unsigned i = 0; i < accepted_rules.size(); ++i) {
		cout << "Accepted Rules for Packet " << i << ": ";
		set<unsigned>::iterator it;
		for (it = accepted_rules[i].begin(); it != accepted_rules[i].end(); ++it)
			cout << *it << " ";
		cout << endl;
	}		
#endif		
	gettimeofday(&c6, NULL);
}	
	cout << "----------------- Kernel execution done -----------------" << endl;
	cout << endl;
    
	seconds  = c2.tv_sec  - c1.tv_sec;
	useconds = c2.tv_usec - c1.tv_usec;
    t_NFAload   = ((double)seconds * 1000 + (double)useconds/1000.0);
	
	seconds  = c4.tv_sec  - c3.tv_sec;
	useconds = c4.tv_usec - c3.tv_usec;
    t_in     = ((double)seconds * 1000 + (double)useconds/1000.0);
#ifdef DEBUG	
	seconds  = c5.tv_sec  - c4.tv_sec;
	useconds = c5.tv_usec - c4.tv_usec;
    t_exec   = ((double)seconds * 1000 + (double)useconds/1000.0);
	
	seconds  = c6.tv_sec  - c5.tv_sec;
	useconds = c6.tv_usec - c5.tv_usec;
    t_out   = ((double)seconds * 1000 + (double)useconds/1000.0);
#endif	
    printf("Execution times: NFA loading (from text): %lf(ms), Input stream loading: %lf(ms), GPU mem alloc: %lf(ms), GPU kernel execution: %lf(ms), Result collecting: %lf(ms)\n", t_NFAload, t_in, t_alloc, t_kernel, t_collect);
	
#ifdef DEBUG
	//Write timing result and blocksize to file
	double t_NFAs[7];
	t_NFAs[0] = t_alloc;
	t_NFAs[1] = t_kernel;
	t_NFAs[2] = t_collect;
	t_NFAs[3] = t_NFAload;
	t_NFAs[4] = t_in;
	t_NFAs[5] = t_exec;
	t_NFAs[6] = t_out;
	
	if (timing_filename != NULL)
		fp_timing.write((char *)t_NFAs, 7*sizeof(double));
	if (blksiz_filename != NULL)
		fp_blksiz.write((char *)&blockSize, sizeof(int));
#endif
	
	// close the file
	fp.close();

#ifdef DEBUG	
	if (timing_filename != NULL)
		fp_timing.close();
	if (blksiz_filename != NULL)
		fp_blksiz.close();
#endif
	
    for (unsigned int i = 0; i < n_subsets; i++) {
		nfa_vec[i]->get_mutable_persistent_states().free_host();
		nfa_vec[i]->get_mutable_initial_states().free_host();
		nfa_vec[i]->get_accept_states().free_host();
		nfa_vec[i]->tmp_states_->free_host();
		nfa_vec[i]->free_hostmem();
		delete nfa_vec[i]->tmp_states_;
		delete nfa_vec[i];
	}
	
	cudaDeviceReset();//Explicitly destroys and cleans up all resources associated with the current device in the current process. Note that this function will reset the device immediately. It is the caller's responsibility to ensure that the device is not being accessed by any other host threads from the process when this function is called.
	//To prevent strange memory leak in some machines (or drivers)
	
#ifdef DEBUG	
	free(rulestartvec);
#endif	
	free(trans_per_sym);
	
	return 0;
}

bool ParseCommandLine(int argc, char *argv[])
{
	int CurrentItem = 1;
	int retVal;

	while (CurrentItem < argc) {

		if (strcmp(argv[CurrentItem], "-a") == 0)
		{
			CurrentItem++;
			pattern_name=argv[CurrentItem];
			CurrentItem++;
			continue;
		}

		if (strcmp(argv[CurrentItem], "-i") == 0)
		{
			CurrentItem++;
			char *trace_filename = NULL;
			trace_filename=argv[CurrentItem];
			cfg.set_trace_file_name(trace_filename);
			CurrentItem++;
			continue;
		}

		if (strcmp(argv[CurrentItem], "-p") == 0)
		{
			CurrentItem++;
			unsigned int parallel_packets;
			retVal = sscanf(argv[CurrentItem],"%d", &parallel_packets);
			cfg.set_parallel_packets(parallel_packets);
			if(retVal!=1 || parallel_packets < 1 ){
				printf("Invalid parallel_packets number: %s\n", argv[CurrentItem]);
				return false;
			}
			CurrentItem++;
			continue;
		}
		
		if (strcmp(argv[CurrentItem], "-T") == 0)
		{
			CurrentItem++;
			unsigned int threads_per_block;
			retVal = sscanf(argv[CurrentItem],"%d", &threads_per_block);
			cfg.set_threads_per_block(threads_per_block);
			if(retVal!=1 || threads_per_block < 1 ){
				printf("Invalid THREADS_PER_BLOCK number: %s\n", argv[CurrentItem]);
				return false;
			}
			CurrentItem++;
			continue;
		}
#ifdef DEBUG		
		if (strcmp(argv[CurrentItem], "-f") == 0)
		{
			CurrentItem++;
			timing_filename=argv[CurrentItem];
			CurrentItem++;
			continue;
		}	
		
		if (strcmp(argv[CurrentItem], "-fg") == 0)
		{
			CurrentItem++;
			blksiz_filename=argv[CurrentItem];
			CurrentItem++;
			continue;
		}	
#endif
        if (strcmp(argv[CurrentItem], "-g") == 0)
		{
			CurrentItem++;
			unsigned int blocks_y;
			retVal = sscanf(argv[CurrentItem],"%d", &blocks_y);
			cfg.set_blocks_y(blocks_y);
			if(retVal!=1 || blocks_y < 1 ){
				printf("Invalid BLOCKS_Y number: %s\n", argv[CurrentItem]);
				return false;
			}
			CurrentItem++;
			continue;
		}
		
		if (strcmp(argv[CurrentItem], "-N") == 0)
			{
				CurrentItem++;
				retVal = sscanf(argv[CurrentItem],"%d", &total_rules);
				if(retVal!=1 || total_rules < 1 ){
					printf("Invalid TOTAL_RULES number: %s\n", argv[CurrentItem]);
					return false;
				}
				CurrentItem++;
				continue;
		}

		if (strcmp(argv[CurrentItem], "-O") == 0)
			{
				CurrentItem++;
				retVal = sscanf(argv[CurrentItem],"%d", &blksiz_tuning);
				if(retVal!=1 || blksiz_tuning > 1 ){
					printf("Invalid blksiz_tuning param: %s\n", argv[CurrentItem]);
					return false;
				}
				CurrentItem++;
				continue;
		}

		if (strcmp(argv[CurrentItem], "-h") == 0 || strcmp(argv[CurrentItem], "-?") == 0)
		{
			CurrentItem++;
			Usage();
			return false;
		}
	}


	return true;
}


void Usage(void) {
	char string[]= "USAGE: ./nfa_engine [OPTIONS] \n" \
					 "\t-a <file> :   transition graph file name with full directory path (must NOT contain the file extension)\n" \
					 "\t-i <file> :   input file to be processed (with file extension)\n"  \
					 "\t-T <n>    :   number of threads per block (overwritten if block size tuning feature is used)\n" \
					 "\t-g <n>    :   number of graphs (or NFAs)(number of thread blocks in grid.y) to be executed (default: 1)\n" \
					 "\t-p <n>    :   number of parallel packets to be examined (number of thread blocks in grid.x)(defaul: 1)\n"\
					 "\t-N <n>    :   total number of rules (subgraphs)\n" \
					 "\t-O <n>    :   0 - block size tuning not enabled; 1 - block size tuned (optional, default: 0 - not tuning)\n" \
#ifdef DEBUG
					 "\t-f <name> :   timing result filename (optional, default: empty)\n" \
					 "\t-fg <name>:   blocksize filename (optional, default: empty)\n" \
#endif
					 "Ex:\t./nfa_engine -a ../data/Sample_NFA/Sample_NFA -i ../data/random_stream_1MB.input -T 1024 -g 2 -p 1 -N 3072 -O 0\n" \
					 "\t./nfa_engine -a ../data/Sample_NFA/Sample_NFA -i ../data/random_stream_1MB.input -g 2 -p 1 -N 3072 -O 1\n";
	fprintf(stderr, "%s", string);
}

/**
 * Get the size of a file.
 * @return The filesize, or 0 if the file does not exist.
 */
 size_t getFilesize(const char* filename) {
    struct stat st;
    if(stat(filename, &st) != 0) {
        return 0;
    }
    return st.st_size;   
}