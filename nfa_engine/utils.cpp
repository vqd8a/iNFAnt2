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

#include "utils.h"

#include <libgen.h>
#include <netinet/in.h>

#include <iostream>
#include <string>
#include <string.h>
#include <cstdio>
#include <algorithm>

std::pair<std::string, std::string> utils::dir_base_name(std::string path)
{
	std::string base, dir;
	char cpath[path.size()+2];
	strcpy(cpath, path.c_str());
	base = basename(cpath);
	strcpy(cpath, path.c_str());
	dir = dirname(cpath);
	return make_pair(dir, base);
}
