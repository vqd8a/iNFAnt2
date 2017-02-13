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

#ifndef __SPARSE_H__
#define __SPARSE_H__

#include <iterator>
#include <ext/hash_map>
#include <set>
#include <map>

#include "common.h"

namespace std { using namespace __gnu_cxx; }

class Range {
	private:
		std::pair<symbol_t, symbol_t> range;
		
	public:
		Range(symbol_t);
		Range(symbol_t, symbol_t);
		Range(const Range &);
		Range(const std::pair<symbol_t, symbol_t> &);
		
		bool operator<(const Range &other) const;
		
		const std::pair<symbol_t, symbol_t> &get() const;
		std::pair<symbol_t, symbol_t> &set();
		
		bool includes(symbol_t) const;
		std::pair<bool, Range> overlap(const Range &other) const;
		unsigned size() const;
};


class Ranges {
	private:
		std::set<Range> ranges;

	public:
		Ranges();

		Ranges(Range);
		template <class T> Ranges(const T &s) {
			typename T::const_iterator ii;
			for(ii = s.begin(); ii != s.end(); ++ii)
				add(*ii);
			return;
		}

		void add(symbol_t);
		void add(Range);
		void add(const Ranges &);
		std::set<Range> get() const;
		
		bool includes(symbol_t) const;

		void subtract(const Ranges &other);
		bool empty() const;
		unsigned size() const;
		unsigned count() const;

		void check() const;
		void print() const;

		Ranges complement(symbol_t start, symbol_t stop) const;

		friend class const_iterator;

		class const_iterator {
			private:
				const std::set<Range> *original;
				std::set<Range>::const_iterator current;
				symbol_t current_sym;
			public:
				typedef std::forward_iterator_tag iterator_category;
				typedef symbol_t value_type;
				typedef ptrdiff_t difference_type;
				typedef symbol_t *pointer;
				typedef symbol_t &reference;
				const_iterator();
				const_iterator(const Ranges &);
				virtual ~const_iterator() {};
				bool operator !=(const const_iterator &other) const;
				const_iterator &operator=(const const_iterator &other);
				const symbol_t &operator*() const;
				const_iterator &operator++();
				const_iterator &go_to_end();
		};

		const_iterator begin() const;
		const_iterator end() const;
};


Ranges range_intersection(const Ranges & a, const Ranges &b);


class Graph {
	public:
		typedef unsigned int state_t;

	private:
		std::hash_map<state_t, std::hash_map<state_t, Ranges> > incoming;

	public:
		void add(state_t, state_t, Ranges);
		
		Ranges get(state_t, state_t);
		
		void project(std::hash_map<symbol_t, std::map<state_t, std::set<state_t> > > &) const;
};

#endif /* __SPARSE_H__ */

