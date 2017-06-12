/*
 *  Copyright (C) 2010  Regents of the University of Michigan
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ERROR_H_
#define _ERROR_H_

#include <string>

extern "C" {
  size_t hts_realloc_or_die(unsigned long, unsigned long, unsigned long, unsigned long, int, void**, char const*);
}


// #ifdef __cplusplus
// extern "C" {
// #endif

void error(const char * msg, ...);
void warning(const char * msg, ...);
void numerror(const char * msg, ...);
void notice(const char * msg, ...);
void catprintf(std::string &s, const char * msg, ...);

// #ifdef __cplusplus
//   };
// #endif


#endif
