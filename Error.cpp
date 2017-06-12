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

#include "Error.h"
#include "pException.h"

#include <string>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <stdarg.h>

void error(const char * msg, ...)
{
  va_list  ap;

  va_start(ap, msg);

  fprintf(stderr, "\nFATAL ERROR - \n");
  vfprintf(stderr, msg, ap);
  fprintf(stderr, "\n\n");

  va_end(ap);

  throw pexception;
  //exit(EXIT_FAILURE);
}

void warning(const char * msg, ...)
{
  va_list  ap;

  va_start(ap, msg);

  fprintf(stderr,"\n\aWARNING - \n");
  vfprintf(stderr, msg, ap);
  fprintf(stderr,"\n");

  va_end(ap);
}

void numerror(const char * msg , ...)
{
  va_list  ap;

  va_start(ap, msg);

  fprintf(stderr,"\nFATAL NUMERIC ERROR - ");
  vfprintf(stderr, msg, ap);
  fprintf(stderr,"\n\n");

  va_end(ap);

  exit(EXIT_FAILURE);
}

void notice(const char * msg, ...) {
  va_list ap;
  va_start(ap, msg);

  time_t current_time;
  char buff[255];
  current_time = time(NULL);

  strftime(buff, 120, "%Y/%m/%d %H:%M:%S", localtime(&current_time));

  fprintf(stderr,"NOTICE [%s] - ", buff);
  vfprintf(stderr, msg, ap);
  fprintf(stderr,"\n");

  va_end(ap);
}

void catprintf(std::string &s, const char * msg, ...)
{
  va_list ap;

  va_start(ap, msg);

  char buf[1000];
  vsprintf(buf, msg, ap);

  s += buf;
  va_end(ap);
}
