// Copyright (C) 2010  
// Pierluigi Rolando (pierluigi.rolando@polito.it)
// Netgroup - DAUIN - Politecnico di Torino
//
// Niccolo' Cascarano (niccolo.cascarano@polito.it)
// Netgroup - DAUIN - Politecnico di Torino
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

#ifndef __UTILS_H__
#define __UTILS_H__

#include "common.h"
#include <arpa/inet.h>
#include <libgen.h>
#include <netinet/in.h>

#include <iostream>
#include <string>
#include <cstdio>
#include <algorithm>

namespace utils {

//template<typename in_sym, typename out_sym>
//bool read_line(istream &, in_sym &, in_sym &, out_sym &);

std::pair<std::string, std::string> dir_base_name(std::string path);


template<typename in_sym, typename out_sym>
bool read_line(std::istream &in, in_sym &start, in_sym &end, out_sym &out)
{
	std::string s;
	
	if(!std::getline(in, s))
		return 0;
	
	if(std::sscanf(s.c_str(), "%lu|%lu |> %hu", &start, &end, &out) != 3)
		if(std::sscanf(s.c_str(), "%lu |> %hu", &start, &out) == 2)
			end = start;
		else
			return 0;
		
		return 1;
}
}
#define CUDA_CHECK(a, b) \
	do { \
		if((a) != cudaSuccess) { \
			fprintf(stderr, b ": %s\n", cudaGetErrorString((a))); \
			std::exit(NFA_ERROR); \
		} \
	} while(0)

#define MALLOC_CHECK(a, b) \
	if(a == NULL) { \
			fprintf(stderr, b); \
			return NFA_ERROR; \
	}

#define error(a, b) \
	do { \
		std::cerr << (a) << std::endl; \
		std::exit(b); \
	} while(0)

#endif /* __UTILS_H__ */

