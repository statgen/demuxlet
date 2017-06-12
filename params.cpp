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

#include "params.h"
#include "Constant.h"
#include "Error.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctype.h>
#include <stdarg.h>

int param::nameCol = 30;
int param::statusCol = 15;
int param::helpCol = 20;

// constructor -- assign initial values
param::param(char c, const char * desc, void * v, const char* help) 
  : ch((char)tolower(c)),
    description(desc == NULL ? "" : desc),
    var(v),
    helpstring(help == NULL ? "" : help),
    errors(NULL),
    messages(NULL)
{
}

// ???
int param::TranslateExtras(const char * , const char *)
{
  return -1;
}

// set parameter speficific errors
void param::error(const char * format, ...)
{
  va_list ap;
  va_start(ap, format);
  char buf[65535];
  vsprintf(buf, format, ap);
  va_end(ap);

  if (errors == NULL) {
    ::error(buf);
  }
  else {
    (*errors) += buf;
  }
}

// set parameter speficific errors
void param::message(const char * format, ...)
{
  va_list ap;
  va_start(ap, format);
  char buf[65535];
  vsprintf(buf, format, ap);
  va_end(ap);

  if (messages == NULL) {
    ::printf(buf);
  }
  else {
    (*messages) += buf;
  }
}

// constructor of long parameter
longParams::longParams(const char * desc, longParamList * lst)
  : param('-', desc, NULL)
{
  list = lst;

  index.clear();
  group_len = 0;
  name_len = 0;

  longParamList * ptr = list + 1;  // first parameter

  while (ptr->desc != NULL) {  // unless it is non-informative 
    if (ptr->value != NULL) {  // if it is a real parameter
      index[ptr->desc] = ptr;  // map string to param
      int tmp = strlen(ptr->desc);   // if it is just a group
      if (tmp > name_len) name_len = tmp; // modify group_len      
    }
    else {
      int tmp = strlen(ptr->desc);   // if it is just a group
      if (tmp > group_len) group_len = tmp; // modify group_len
    }
    
    ptr++;
  }

  precision = 2;  // ???
}

// cstr is [option_name] and extra is value
int longParams::TranslateExtras(const char * cstr, const char * extras)
{
  std::map<std::string, longParamList*>::iterator it = index.find(cstr);

  if ( it != index.end() ) { // keyword was found
    longParamList* ptr = it->second;
    switch(ptr->type) {
    case LP_BOOL_PARAM:
      if ( ptr->touched )          error("[E:%s:%d %s] Redundant use of option --%s or its exclusive neighbor is not allowed",__FILE__,__LINE__,__FUNCTION__,cstr);
      if ( *(bool*)ptr->value == true )
	message("Option --%s specified when it is [ON] by default. The argument has no effect and the variable will be still [ON]",cstr);
      *(bool *) ptr->value = true;
      ptr->touched = true;
      if ( ptr->exclusive ) {
	for (int i = -1; ptr[i].exclusive; i--) {
	  *(bool *)ptr[i].value = false;
	  ptr[i].touched = true;
	}
	for (int i =  1; ptr[i].exclusive; i++) {
	  *(bool *)ptr[i].value = false;
	  ptr[i].touched = true;
	}
      }
      return 0;
    case LP_INT_PARAM:
      if ( !CheckInteger(extras) ) error("[E:%s:%d %s] Invalid argument --%s %s. Integer was expected",__FILE__,__LINE__,__FUNCTION__, cstr, extras);
      else if ( ptr->touched )          error("[E:%s:%d %s] Redundant use of option --%s is not allowed",__FILE__,__LINE__,__FUNCTION__,cstr);
      else {
	*(int *) ptr->value = atoi(extras);
	ptr->touched = true;
      }
      return 1;
    case LP_DOUBLE_PARAM:
      if ( !CheckDouble(extras) ) error("[E:%s:%d %s] Invalid argument --%s %s. Double was expected",__FILE__,__LINE__,__FUNCTION__, cstr, extras);
      else if ( ptr->touched )         error("[E:%s:%d %s] Redundant use of option --%s is not allowed",__FILE__,__LINE__,__FUNCTION__,cstr);
      else {
	*(double *) ptr->value = atof(extras);
	ptr->touched = true;      
      }
      return 1;
    case LP_STRING_PARAM:
      if ( extras == NULL ) error("[E:%s:%d %s] Invalid argument --%s %s. String was expected",__FILE__,__LINE__,__FUNCTION__, cstr, extras);
      else if ( ptr->touched )         error("[E:%s:%d %s] Redundant use of option --%s is not allowed",__FILE__,__LINE__,__FUNCTION__,cstr);
      else {
	*(std::string *) ptr->value = extras;
	ptr->touched = true;      
      }
      return 1;
    case LP_MULTI_INT_PARAM:
      if ( !CheckInteger(extras) ) error("[E:%s:%d %s] Invalid argument --%s %s. Integer was expected",__FILE__,__LINE__,__FUNCTION__, cstr, extras);
      else
	((std::vector<int> *) ptr->value)->push_back(atoi(extras));
      return 1;
    case LP_MULTI_DOUBLE_PARAM:
      if ( !CheckDouble(extras) ) error("[E:%s:%d %s] Invalid argument --%s %s. Double was expected",__FILE__,__LINE__,__FUNCTION__, cstr, extras);
      else 
	((std::vector<double> *) ptr->value)->push_back(atof(extras));
      return 1;
    case LP_MULTI_STRING_PARAM:
      if ( extras == NULL ) 
	error("[E:%s:%d %s] Invalid argument --%s %s. String was expected",__FILE__,__LINE__,__FUNCTION__, cstr, extras);
      else 
      ((std::vector<std::string> *) ptr->value)->push_back(extras);
      return 1;
    default:
      return -1; // ignore?
    }
  }
  else {
    return -1; // ignore?
  }
}

// Print the status of parameter List
void longParams::Status(longParamList * ptr, int & line_len, bool & need_a_comma)
{
  std::string state;
  int line_start = group_len ? group_len + 5 : 0;
  int i;

  if (ptr->value == NULL) {  // if group parameter, end previous group (if exists) 
    fprintf(stderr, "%s %*s :", need_a_comma ? "\n" : "", group_len + 2, ptr->desc);
    need_a_comma = false;
    line_len = line_start;
  }
  else {                     // otherwise, print argument name and 
    switch(ptr->type) {
    case LP_BOOL_PARAM:
      state = * (bool *) ptr->value ? " [ON]" : "";
      break;
    case LP_INT_PARAM:
      if (((* (int *) ptr->value == 1) && (ptr->exclusive)) || (* (int *) ptr->value == 0))
	state = *(int *) ptr->value ? " [ON]" : "";
      else
	catprintf(state," [%d]",*(int*)ptr->value);
      break;
    case LP_DOUBLE_PARAM:
      if (* (double *) ptr->value != _NAN_) {
	  double value = * (double *) ptr->value;
	  
	  state = " [";
	  if (value == 0.0 || value >= 0.01)
	    catprintf(state, "%.*f", precision, value);
	  else
	    catprintf(state, "%.1e", value);
	  state += ']';
	}
      else
	state = "";
      break;
    case LP_STRING_PARAM:
      if ( ((std::string*) ptr->value)->empty() ) 
	state = "";
      else
	state = " [" + * (std::string *) ptr->value + "]";
      break;
    case LP_MULTI_INT_PARAM:
      {
	std::vector<int>* v = (std::vector<int>*) ptr->value;
	if ( v->empty() ) 
	  state = "";
	else {
	  state = " [";
	  for(i=0; i < (int)v->size(); ++i) {
	    if ( i > 0 ) 
	      catprintf(state, ", %d", v->at(i));
	    else
	      catprintf(state, "%d", v->at(i));
	  }
	  state += ']';
	}
	break;
      }
    case LP_MULTI_DOUBLE_PARAM:
      {
	std::vector<double>* v = (std::vector<double>*) ptr->value;
	if ( v->empty() ) 
	  state = "";
	else {
	  state = " [";
	  for(i=0; i < (int)v->size(); ++i) {
	    if ( i > 0 ) 
	      state += ", ";
	    
	    if (v->at(i) == 0.0 || v->at(i) >= 0.01)
	      catprintf(state, "%.*f", precision, v->at(i));
	    else
	      catprintf(state, "%.1e", v->at(i));
	  }
	  state += ']';
	}
	break;
      }
    case LP_MULTI_STRING_PARAM:
      {
	std::vector<std::string>* v = (std::vector<std::string>*) ptr->value;
	if ( v->empty() ) 
	  state = "";
	else {
	  state = " [";
	  for(i=0; i < (int)v->size(); ++i) {
	    if ( i > 0 ) 
	      state += ", ";
	    state += v->at(i);
	  }
	  state += ']';
	}
	break;
      }
    default:
      error("[E:%s:%d %s] Cannot recognize the parameter type %d",__FILE__,__LINE__,__FUNCTION__,ptr->type);
    }
    
    int item_len = 3 + strlen(ptr->desc) + need_a_comma + state.size();
    
    if (item_len + line_len > 78 && line_len > line_start)
	{
	  line_len = line_start;
	  fprintf(stderr, "%s\n%*s", need_a_comma ? "," : "", line_len,  "");
	  need_a_comma = 0;
	  item_len -= 1;
	}

      fprintf(stderr, "%s --%s%s", need_a_comma ? "," : (need_a_comma = true, ""),
	      ptr->desc, state.c_str());

      need_a_comma = true;
      line_len += item_len;
  }
}

// Print the status of parameter List
void longParams::HelpMessage(longParamList * ptr)
{
  std::string state;
  int i;

  if (ptr->value == NULL) {  // if group parameter, end previous group (if exists) 
    fprintf(stderr, "\n%s%s%s\n", ptr->desc, ptr->help ? " - " : "", ptr->help ? ptr->help : "");
  }
  else {                     // otherwise, print argument name and 
    switch(ptr->type) {
    case LP_BOOL_PARAM:
      state = * (bool *) ptr->value ? " [FLG: ON]" : " [FLG: OFF]";
      break;
    case LP_INT_PARAM:
      if (((* (int *) ptr->value == 1) && (ptr->exclusive)) || (* (int *) ptr->value == 0))
	state = *(int *) ptr->value ? " [INT: ON]" : " [INT: 0]";
      else
	catprintf(state," [INT: %d]",*(int*)ptr->value);
      break;
    case LP_DOUBLE_PARAM:
      if (* (double *) ptr->value != _NAN_) {
	  double value = * (double *) ptr->value;
	  
	  state = " [FLT: ";
	  if (value == 0.0 || value >= 0.01)
	    catprintf(state, "%.*f", precision, value);
	  else
	    catprintf(state, "%.1e", value);
	  state += ']';
      }
      else
	state = " [FLT: NaN]";
      break;
    case LP_STRING_PARAM:
      if ( ((std::string*) ptr->value)->empty() ) 
	state = " [STR: ]";
      else
	state = " [STR: " + * (std::string *) ptr->value + "]";
      break;
    case LP_MULTI_INT_PARAM:
      {
	std::vector<int>* v = (std::vector<int>*) ptr->value;
	if ( v->empty() ) 
	  state = " [V_INT: ]";
	else {
	  state = " [V_INT: ";
	  for(i=0; i < (int)v->size(); ++i) {
	    if ( i > 0 ) 
	      catprintf(state, ", %d", v->at(i));
	    else
	      catprintf(state, "%d", v->at(i));
	  }
	  state += ']';
	}
	break;
      }
    case LP_MULTI_DOUBLE_PARAM:
      {
	std::vector<double>* v = (std::vector<double>*) ptr->value;
	if ( v->empty() ) 
	  state = " [V_FLT: ]";
	else {
	  state = " [V_FLT: ";
	  for(i=0; i < (int)v->size(); ++i) {
	    if ( i > 0 ) 
	      state += ", ";
	    
	    if (v->at(i) == 0.0 || v->at(i) >= 0.01)
	      catprintf(state, "%.*f", precision, v->at(i));
	    else
	      catprintf(state, "%.1e", v->at(i));
	  }
	  state += ']';
	}
	break;
      }
    case LP_MULTI_STRING_PARAM:
      {
	std::vector<std::string>* v = (std::vector<std::string>*) ptr->value;
	if ( v->empty() ) 
	  state = " [V_STR: ]";
	else {
	  state = " [V_STR: ";
	  for(i=0; i < (int)v->size(); ++i) {
	    if ( i > 0 ) 
	      state += ", ";
	    state += v->at(i);
	  }
	  state += ']';
	}
	break;
      }
    default:
      error("[E:%s:%d %s] Cannot recognize the parameter type %d",__FILE__,__LINE__,__FUNCTION__,ptr->type);
    }
    

    fprintf(stderr, "  --%-*s%-*s%s%s%s\n", name_len, ptr->desc, param::helpCol, state.c_str(), " : ", ptr->help ? ptr->help : "", ptr->exclusive ? (ptr->help ? " (EXCLUSIVE PARAMETER)" : "(EXCLUSIVE PARAMETER)") : "");
  }
}

// print the status of the parameter
void longParams::HelpMessage()
{
  if (!description.empty() && description[0] != 0)  // group option
    fprintf(stderr, "\n%s\n\n", description.c_str());    
    //fprintf(stderr, "\n%s - %s\n", description.c_str(), helpstring.c_str());

  // for the rest of the group, print parameters
  for (longParamList * ptr = list + 1; ptr->desc != NULL; ptr++)  
    HelpMessage(ptr);

  fprintf(stderr, "\n");
}

// print the status of the parameter
void longParams::Status()
{
  if (!description.empty() && description[0] != 0)  { // group option
    fprintf(stderr, "\n%s\n\n", description.c_str());

    fprintf(stderr, "The following parameters are available. Ones with \"[]\" are in effect:\n");
  }

  bool need_a_comma = false;
  int  line_len = 0;

  // for the rest of the group, print parameters
  for (longParamList * ptr = list + 1; ptr->desc != NULL; ptr++)  
    Status(ptr, line_len, need_a_comma);

  fprintf(stderr, "\n");
}

// Add parameter
void paramList::Add(param * p)
{
  p->SetErrorBuffer(errors);
  p->SetMessageBuffer(messages);
  pl.push_back(p);
};

// Read parameters from argument
void paramList::Read(int argc, char ** argv, int start)
{
  int i, j;
  // iterate from first argument
  for (i=start; i < argc; i++) {
    bool success = false;

    if (argv[i][0] == '-' && argv[i][1]) { // first is -, second is non-null
      if ( ( argv[i][1] == 'h' ) || ( strcmp(argv[i],"--help") == 0 ) ) { // printing help requested
	HelpMessage();
	if ( messages.empty() ) 
	  fprintf(stderr,"NOTES:\n");
	fprintf(stderr, "When --help was included in the argument. The program prints the help message but do not actually run\n");
	exit(1);
      }

      for (j=0; j<(int)pl.size(); j++) {      // compare with all available parameters
	success = tolower(argv[i][1]) == pl[j]->ch;  // option should match
	if (success) {
	  // see if it can be parsed using two consecutive arguments
	  //if ((i+1 < argc) && pl[j]->TranslateExtras(argv[i]+2, argv[i+1]))
	  int ret = pl[j]->TranslateExtras(argv[i]+2, argv[i+1]);
	  if ( ret >= 0 ) {
	    if ( ret == 1 ) i++;
	    break;
	  }
	}
      }

      if ( j == (int)pl.size() ) {
	catprintf(errors, "Command line parameter %s (#%d) not recognized\n", argv[i], i);
      }
    }
    else {
      catprintf(errors, "Cannot correspond command line parameter %s (#%d) to any of the options\n", argv[i], i);
    }
  }
}



int paramList::ReadWithTrailer(int argc, char ** argv, int start)
{
  int last_success = start - 1;
  bool split = false;

  for (int i=start; i < argc; i++)
    {
      bool success = false;

      if (argv[i][0] == '-' && argv[i][1])
	for (int j=0; j<(int)pl.size(); j++)
	  {
	    success = tolower(argv[i][1]) == pl[j]->ch;

	    if (success)
	      {
		if ((i+1 < argc) && pl[j]->TranslateExtras(argv[i]+2, argv[i+1]))
		  split = true;
		break;
	      }
	  }

      if (success)
	for (last_success++; last_success < i; last_success++)
	  catprintf(errors,"Command line parameter %s (#%d) ignored\n",
			  argv[last_success], last_success);

      if (split)
	{
	  split = false;
	  last_success++;
	  i++;
	}
    }

  return last_success;
};

void paramList::HelpMessage()
{
  fprintf(stderr, "\nDetailed instructions of parameters are available. Ones with \"[]\" are in effect:\n");

  for (int i=0; i<(int)pl.size(); i++)
    pl[i]->HelpMessage();

  fprintf(stderr, "\n");

  if (errors.size())
    {
      error("[E:%s:%d %s] Problems encountered parsing command line:\n\n%s",__FILE__,__LINE__,__FUNCTION__,
		errors.c_str());
      errors.clear();
    }

  if (messages.size())
    {
      ::printf("NOTES:\n%s\n", messages.c_str());
      //messages.clear();
    }
}

// Print the total parameter list
void paramList::Status()
{
  //fprintf(stderr, "\nThe following parameters are available. Ones with \"[]\" are in effect:\n");

  for (int i=0; i<(int)pl.size(); i++)
    pl[i]->Status();

  fprintf(stderr, "\nRun with --help for more detailed help messages of each argument.\n");  
  fprintf(stderr, "\n");

  if (errors.size())
    {
      error("[E:%s:%d %s] Problems encountered parsing command line:\n\n%s",__FILE__,__LINE__,__FUNCTION__,
		errors.c_str());
      errors.clear();
    }

  if (messages.size())
    {
      ::printf("NOTES:\n%s\n", messages.c_str());
      messages.clear();
    }
}

//

paramList::~paramList()
{
  for (int i = 0; i < (int)pl.size(); i++)
    delete pl[i];
};

bool param::CheckInteger(const char * value)
{
  if ( value == NULL ) return false;
  if (value[0] != '+' && value[0] != '-' &&
      (value[0] < '0' || value[0] > '9'))
    return false;

  int pos = 1;
  while (value[pos] != 0)
    if (value[pos] < '0' || value[pos] > '9')
      return false;
    else
      pos++;

  return true;
}

bool param::CheckDouble(const char * value)
{
  if ( value == NULL ) return false;
  if (value[0] != '+' && value[0] != '-' && value[0] != '.' &&
      (value[0] < '0'  || value[0] > '9'))
    {
      return false;
    }

  bool decimal = value[0] == '.';

  for (int pos = 1; value[pos] != 0; pos++)
    {
      if (value[pos] < '0' || value[pos] > '9')
	{
	  if (!decimal && value[pos] == '.')
	    {
	      decimal = true;
	    }
	  else if (value[pos] == 'e' || value[pos] == 'E')
	    {
	      return CheckInteger(value + pos + 1);
	    }
	}
    }

  return true;
}
