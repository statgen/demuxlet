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

#ifndef __PARAM_H__
#define __PARAM_H__

#include <ctype.h>
#include <stddef.h>
#include <vector>
#include <map>
#include <string>

typedef std::map<std::string,void*> StringMap;
typedef std::map<std::string,void*>::iterator StringMapIterator;

class paramList;  // list of parameters

class param  // represent each parameter - an abstract class
{
 protected:
  char ch;                  // single character to represent the option
  std::string description;  // full name of the option
  void * var;               // pointer to store actual value
  std::string helpstring;   // detailed description of the option

  static int nameCol;       // length of name column
  static int statusCol;     // length of status column
  static int helpCol;       // length of the help column

  // get option by full string
  virtual int TranslateExtras(const char * value, const char * extras);

  static bool CheckInteger(const char * value); 
  static bool CheckDouble(const char * value);

  std::string * errors;
  std::string * messages;

 public:

  // constructor 
  param(char c, const char * desc, void * v, const char * help = NULL);

  // destructor
  virtual ~param() {}

  // Read argn-th argument and assing values
  //virtual bool Read(int argc, char ** argv, int argn) {}

  // virtual function which prints the name and values
  virtual void Status() = 0;
  virtual void HelpMessage() = 0;

  // modify nameCol
  static void SetNameLen(int len)
  {
    nameCol = len;
  }

  // modify statusCol
  static void SetStatusLen(int len)
  {
    statusCol = len;
  }

  // set error buffer
  void SetErrorBuffer(std::string & buffer)
  {
    errors = &buffer;
  }

  // set error buffer
  void SetMessageBuffer(std::string & buffer)
  {
    messages = &buffer;
  }

  // function printing warning
  void error(const char * format, ...);
  void message(const char * format, ...);

  // give full access to paramList class
  friend class paramList;
};

// contents for long parameter list
struct longParamList
{
  const char * desc;
  void * value;
  bool exclusive;
  int  type;
  bool touched;
  const char * help;
};

#define LP_BOOL_PARAM            1
#define LP_INT_PARAM             2
#define LP_DOUBLE_PARAM          3
#define LP_STRING_PARAM          4
#define LP_MULTI_INT_PARAM      12
#define LP_MULTI_DOUBLE_PARAM   13
#define LP_MULTI_STRING_PARAM   14

#define BEGIN_LONG_PARAMS(array)   longParamList array[] = {	\
    { NULL,  NULL,      false,  0, 0, NULL},

#define LONG_PARAM_GROUP(label, help)                 { label, NULL,      false,  0, 0, help},

#define LONG_PARAM(label,boolptr, help)               { label, boolptr,   false,  LP_BOOL_PARAM, 0, help},
#define EXCLUSIVE_PARAM(label,boolptr, help)          { label, boolptr,   true,   LP_BOOL_PARAM, 0, help},
#define LONG_INT_PARAM(label,intptr, help)             { label, intptr,    false,  LP_INT_PARAM, 0, help},
#define LONG_SMARTINT_PARAM(label,intptr, help)        { label, intptr,    true,   LP_INT_PARAM, 0, help},
#define LONG_DOUBLE_PARAM(label,doubleptr, help)       { label, doubleptr, false,  LP_DOUBLE_PARAM, 0, help},
#define LONG_STRING_PARAM(label,stringptr, help)       { label, stringptr, false,  LP_STRING_PARAM, 0, help},

#define LONG_MULTI_INT_PARAM(label, intvecptr, help)    { label, intvecptr, false, LP_MULTI_INT_PARAM, 0, help == NULL ? label : help},
#define LONG_MULTI_DOUBLE_PARAM(label, dblvecptr, help) { label, dblvecptr, false, LP_MULTI_DOUBLE_PARAM, 0, help == NULL ? label : help},
#define LONG_MULTI_STRING_PARAM(label, strvecptr, help) { label, strvecptr, false, LP_MULTI_STRING_PARAM, 0, help == NULL ? label : help},

#define END_LONG_PARAMS()                    { NULL,  NULL,      false,  0, 0, 0 }};

// long parameter class
class longParams : public param
{
 public:
  longParams(const char * desc, longParamList * list);  // constructor

  virtual void Status();                     // Print the status
  virtual void HelpMessage();                     // Print the status

  longParams * SetPrecision(int precision)   // Set precision for output
  {
    this->precision = precision;

    return this;
  }

 protected:
  std::map<std::string, longParamList*> index;
  std::map<std::string, longParamList*> legacyIndex;

  longParamList * list;
  int group_len;
  int name_len;
  int precision;

  //virtual void Translate(const char * value);
  virtual int TranslateExtras(const char * value, const char * extras);

  void Status(longParamList * ptr, int & line_len, bool & need_a_comma);
  void HelpMessage(longParamList * ptr);
};

// List of parameters
class paramList
{
 protected:
  bool help;
  std::vector<param*> pl; // vector of pointers;

 public:

 paramList() : help(false) {}

  virtual ~paramList();

  void Add(param * p);

  // Tries to process all command line arguments
  virtual void Read(int argc, char ** argv, int start = 1);

  // Allows for trailing, unprocessed, filenames in the command line
  // The number of translated argv[] items is returned
  virtual int ReadWithTrailer(int argc, char ** argv, int start = 1);

  // Outputs summary of parameter switches and settings
  virtual void Status();
  virtual void HelpMessage();

  // Keeps track of warnings generated during parameter processing
  std::string  errors;
  std::string  messages;
};

#endif // Params.h
