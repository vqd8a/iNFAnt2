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

#include "sparse.h"
#include <iostream>
#include <cassert>
#include <algorithm>


using namespace std;


Range::Range(symbol_t s) {
	range.first = range.second = s;
	return;
}

Range::Range(symbol_t s, symbol_t e) {
	assert(s <= e);
	range.first = s;
	range.second = e;
}

Range::Range(const Range &other) {
	range = other.range;
}

Range::Range(const pair<symbol_t, symbol_t> &other) {
	Range(other.first, other.second);
}

bool Range::operator<(const Range &other) const {
	if(range.first == other.range.first)
		return range.second < other.range.second;
	return range.first < other.range.first;
}

const pair<symbol_t, symbol_t> &Range::get() const {
	return range;
}

pair<symbol_t, symbol_t> &Range::set() {
	return range;
}

bool Range::includes(symbol_t s) const {
	return s >= range.first && s <= range.second;
}

pair<bool, Range> Range::overlap(const Range &other) const {
//	cerr << "overlap fra " << get().first << ", " << get().second << " e " << other.get().first <<
//		", " << other.get().second << ' ';
//
	symbol_t start = max(get().first, other.get().first);
	symbol_t end = min(get().second, other.get().second);

	if(start < end)
		return make_pair(true, Range(start, end));
	else
		return make_pair(false, Range(0, 0));
}

unsigned Range::size() const {
	assert(range.second >= range.first);
	return range.second - range.first + 1;
}


Ranges::Ranges() {
	// NOP
	return;
}

Ranges::Ranges(Range r) {
	ranges.insert(r);
	return;
}

void Ranges::add(symbol_t s) {
	add(Range(s));
	return;
}

void Ranges::add(Range r) {
	set<Range>::iterator ii, jj;
	bool done(false);

	Range copy = r;


	do {
		done = true;

//		cerr << "aggiungo " << r.get().first << ", " << r.get().second << '\n';
		pair<set<Range>::iterator, bool> rr(ranges.insert(r));
		if(!rr.second) {
//			cerr << "non fatto\n";
			break;
		}
		
		ii = rr.first;

		assert(r.get().first == ii->get().first);
		assert(r.get().second == ii->get().second);
	
		if(ii != ranges.begin()) {
//			cerr << "non sono il primo\n";
			--ii;
		}

		jj = ii;

		if(jj != ranges.end()) {
//			cerr << "non sono l'ultimo\n";
			++jj;
		}
	
		// Merge ranges to maintain the invariant
		unsigned w(0);
		for(; jj != ranges.end(); ++jj, ++ii, ++w) {
			assert(ii->get().first <= jj->get().first);

//			cerr << "finisce a " << ii->get().second << " e riprende a " << jj->get().first << '\n';

			if(jj->get().first <= (ii->get().second + 1)) {
//				cerr << "faccio il merge\n";
				r = Range(ii->get().first, max(ii->get().second, jj->get().second));
				ranges.erase(ii);
				ranges.erase(jj);
				done = false;
				break;
			} //else if(w > 2);
			//	break;
		}
	} while(!done);

//	cerr << "esco\n";

//	cerr << "status:\n";
//	set<Range>::iterator zz;
//	for(zz = ranges.begin(); zz != ranges.end(); ++zz)
//		cerr << '[' << zz->get().first << ", " << zz->get().second << "]\n";

	assert(includes(copy.get().first));
	assert(includes(copy.get().second));
	
	return;
}

void Ranges::add(const Ranges &other) {
	set<Range> o = other.get();
	set<Range>::iterator ii;
	for(ii = o.begin(); ii != o.end(); ++ii)
		add(*ii);
	return;
}

set<Range> Ranges::get() const {
	return ranges;
}

bool Ranges::includes(symbol_t s) const {
	set<Range>::iterator ii;
	for(ii = ranges.begin(); ii != ranges.end(); ++ii)
		if(ii->includes(s))
			return true;
	return false;
}

bool Ranges::empty() const {
	return ranges.empty();
}

Ranges Ranges::complement(symbol_t start, symbol_t end) const {
	Ranges result;

	if(empty())
		return Ranges(Range(start, end));

	symbol_t old_limit(start);
	symbol_t new_limit;

	set<Range>::iterator ii;

	for(ii = ranges.begin(); ii != ranges.end(); ++ii) {
		new_limit = ii->get().first;
		if(old_limit < new_limit)
			result.add(Range(old_limit, new_limit));
		old_limit = ii->get().second+1;
	}

	if(old_limit < end)
		result.add(Range(old_limit, end));


	Ranges intersection(range_intersection(*this, result));
	assert(intersection.size() == 0);
	assert(intersection.empty());

	return result;
}

unsigned Ranges::size() const {
	unsigned total(0);
	set<Range>::iterator ii;
	for(ii = ranges.begin(); ii != ranges.end(); ++ii)
		total += ii->size();
	return total;
}

void Ranges::print() const {
	set<Range>::iterator ii;
	for(ii = ranges.begin(); ii != ranges.end(); ++ii)
		cerr << '[' << ii->get().first << ", " << ii->get().second << "]\n";
	return;
}

unsigned Ranges::count() const {
	return ranges.size();
}

void Ranges::subtract(const Ranges &other) {
	if(empty())
		return;

	if(other.empty())
		return;

	symbol_t initial(0);

	set<Range>::iterator ii, jj;

	ii = ranges.end();
	--ii;

	jj = other.ranges.end();
	--jj;

	symbol_t final(max(ii->get().second, jj->get().second));

	Ranges complement(other.complement(initial, final));
	Ranges intersection(range_intersection(*this, complement));

	ranges = intersection.ranges;
	return;
}

#if 0
void Ranges::subtract(const Ranges &other) {
	cerr << "inizio subtract... ";
	set<Range> o(other.get());
	set<Range>::const_iterator ii, jj;
	set<Range> remaining;

	for(ii = ranges.begin(); ii != ranges.end(); ++ii) {
		bool overlapped = false;
		for(jj = o.begin(); jj != o.end(); ++jj) {
			cerr << "confronto " << ii->get().first << ", " << ii->get().second << "con "
				<< jj->get().first << ", " << jj->get().second << '\n';

			bool include_start(jj->includes(ii->get().first));
			bool include_end(jj->includes(ii->get().second));

			if(include_start || include_end) {
				cerr << "both includes, perfect overlap\n";
				overlapped = true;
			}

			if(!include_start && !include_end) {
				cerr << "no includes\n";
				if(ii->get().first < jj->get().first && ii->get().second > jj->get().second) {
					cerr << "reverse inclusion\n";
					overlapped = true;
					remaining.insert(Range(ii->get().first, jj->get().first - 1));
					cerr << "adding " << ii->get().first << ", " << (jj->get().first -1 ) << '\n';
					remaining.insert(Range(jj->get().second + 1, ii->get().second));
					cerr << "adding " << (jj->get().second + 1) << ", " << ii->get().second << '\n';
				}
			} else if(include_start && !include_end) {
				cerr << "includes start but not end\n";
				remaining.insert(Range(jj->get().second + 1, ii->get().second));
				cerr << "adding " << (jj->get().second + 1) << ", " << ii->get().second << '\n';
			} else if(!include_start && include_end) {
				cerr << "includes end but not start\n";
				remaining.insert(Range(ii->get().first, jj->get().first -1 ));
				cerr << "adding " << ii->get().first << ", " << (jj->get().first -1);
			}
		}

		if(!overlapped) {
			cerr << "didn't overlap, inserting full: " << ii->get().first << ", " << ii->get().second << '\n';
			remaining.insert(*ii);
		}
	}

	Ranges inters_empty(range_intersection(Ranges(remaining), Ranges(o)));
	if(!(inters_empty.size() == 0)) {
		cerr << "inters_empty: " << inters_empty.size() << '\n';
		inters_empty.print();
		cerr << "remaining: " << Ranges(remaining).size() << '\n';
		Ranges(remaining).print();
		cerr << "o: " << Ranges(o).size() << '\n';
		Ranges(o).print();
	}
	assert(inters_empty.size() == 0);
	if(!inters_empty.empty())
		cerr << "inters_empty not empty\n";
	assert(inters_empty.empty());

	Ranges intersection(range_intersection(*this, Ranges(o)));

	cerr << '\n';
	cerr << "remaining: " << Ranges(remaining).size() << '\n';
	cerr << "o:         " << Ranges(o).size() << '\n';
	cerr << "this:      " << size() << '\n';
	cerr << "other:     " << other.size() << '\n';
	cerr << "intersect: " << intersection.size() << '\n';

	assert((Ranges(remaining).size() + Ranges(o).size()) == (size() + other.size() - intersection.size()));



	if((Ranges(remaining).size() + intersection.size()) != size()) {
		cerr << "subtract size mismatch: remaining: " << Ranges(remaining).size() << ", size: " << size()
			<< ", intersection: " << intersection.size() << ", delta: " << (signed)(Ranges(remaining).size() + intersection.size() - size())
			<< '\n';
	}

	assert((Ranges(remaining).size() + intersection.size()) == size());

	ranges = remaining;
	cerr << "fine\n";
	return;
}
#endif

Ranges::const_iterator Ranges::begin() const {
	return const_iterator(*this);
}

Ranges::const_iterator Ranges::end() const {
	return begin().go_to_end();
}

Ranges::const_iterator::const_iterator(const Ranges &r) : original(&r.ranges) {
	current = original->begin();
	if(!original->empty())
		current_sym = current->get().first;
	else
		current_sym = 0;
}

Ranges::const_iterator::const_iterator() : original(0) {
	// Nop
}

Ranges::const_iterator &Ranges::const_iterator::go_to_end() {
	current = original->end();
	if(!original->empty()) {
		--current;
		current_sym = current->get().second + 1;
		++current;
	} else
		current_sym = 0;
	return *this;
}

bool Ranges::const_iterator::operator!=(const Ranges::const_iterator &other) const {
	if(!(current != other.current) && current_sym == other.current_sym)
		return false;
	return true;
}

const symbol_t &Ranges::const_iterator::operator*() const {
	return current_sym;
}

Ranges::const_iterator &Ranges::const_iterator::operator++() {
	if((current != original->end())) {
		if(((++current_sym) > current->get().second)) {
			if(++current != original->end())
				current_sym = current->get().first;
		}
	}
	return *this;
}

Ranges::const_iterator &Ranges::const_iterator::operator=(const const_iterator &other) {
	original = other.original;
	current = other.current;
	current_sym = other.current_sym;
	return *this;
}


Ranges range_intersection(const Ranges &a, const Ranges &b)
{
//	cerr << "intersection, sizes: " << a.size() << '/' << a.count()
//		<< ", " << b.size() << '/' << b.count() << "... ";

	set<Range> first(a.get());
	set<Range> second(b.get());

	set<Range>::const_iterator ii, jj;

	Ranges result;

	for(ii = first.begin(); ii != first.end(); ++ii) {
		bool found(false);
		for(jj = second.begin(); jj != second.end(); ++jj) {
			pair<bool, Range> overlap(ii->overlap(*jj));
			if(overlap.first) {
				found = true;
				result.add(overlap.second);
			} //else if(found)
			//	break;
		}

//		cerr << "status result: " << result.size() << '/' << result.count() << '\n';
	}

//	cerr << "done\n";

	result.check();
//	cerr << "status result: " << result.size() << '/' << result.count() << '\n';

	return result;
}

void Ranges::check() const {
	set<Range>::const_iterator ii, jj;
	unsigned k = 0;

//	cerr << "inizio check\n";

	jj = ranges.begin();
	ii = ranges.begin();

	if(jj != ranges.end())
		++jj;

	for(; jj != ranges.end(); ++jj, ++ii, ++k) {
//		cerr << "check: " << k << '/' << ranges.size() << '\n';

		if(ii->get().second >= jj->get().first)
			cerr << "Invariant 1 violation: " << ii->get().second << ", " << jj->get().first << '\n';
		if(ii->get().second + 1 >= jj->get().first)
			cerr << "Invariant 2 violation: " << ii->get().second << ", " << jj->get().first << '\n';
		assert(ii->get().second < jj->get().first);
		assert(ii->get().second + 1 < jj->get().first);
	}
	
	return;
}


void Graph::add(state_t src, state_t dst, Ranges r) {
	hash_map<state_t, hash_map<state_t, Ranges> >::iterator ii;
	hash_map<state_t, Ranges>::iterator jj;
	
	if((ii = incoming.find(dst)) == incoming.end()) {
		hash_map<state_t, Ranges> hm;
		ii = incoming.insert(make_pair(dst, hm)).first;
	}
	
	if((jj = ii->second.find(src)) == ii->second.end()) {
		Ranges rr;
		jj = ii->second.insert(make_pair(src, rr)).first;
	}
	
	jj->second.add(r);
	
	return;
}
		
Ranges Graph::get(state_t src, state_t dst) {
	Ranges empty;
	
	hash_map<state_t, hash_map<state_t, Ranges> >::iterator ii;
	hash_map<state_t, Ranges>::iterator jj;
	
	if((ii = incoming.find(dst)) != incoming.end() && (jj = ii->second.find(src)) != ii->second.end())
		return jj->second;
	else
		return empty;
}

void Graph::project(hash_map<symbol_t, map<state_t, set<state_t> > > &hm) const {
	hm.clear();

	hash_map<state_t, hash_map<state_t, Ranges> >::const_iterator ii;
	for(ii = incoming.begin(); ii != incoming.end(); ++ii) {
		state_t dst(ii->first);
		hash_map<state_t, Ranges>::const_iterator jj;
		for(jj = ii->second.begin(); jj != ii->second.end(); ++jj) {
			state_t src(jj->first);
			
			set<Range> sr(jj->second.get());
			set<Range>::iterator kk;
			for(kk = sr.begin(); kk != sr.end(); ++kk)
				for(symbol_t i = kk->get().first; i <= kk->get().second; ++i)
					hm[i][src].insert(dst);
		}
	}
	
	return;
}


#ifdef TEST_MAIN

#include <cstdlib>
#include <ctime>
#include <vector>


int main(void)
{
	const unsigned transitions = 100;

	Graph g;
	
	vector<Graph::state_t> srcs, dsts;
	vector<symbol_t> syms;
	
	srand(time(0));
	
	for(unsigned i = 0; i < transitions; ++i) {
		srcs.push_back(rand());
		dsts.push_back(rand());
		syms.push_back(rand());
		
		g.add(srcs.back(), dsts.back(), Ranges(syms.back()));
	}
	
	for(unsigned i = 0; i < transitions; ++i)
		assert(g.get(srcs.at(i), dsts.at(i)).includes(syms.at(i)));
		
	return EXIT_SUCCESS;
}
#endif /* TEST_MAIN */
