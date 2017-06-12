#ifndef __PLUG_EXCEPTION_H
#define __PLUG_EXCEPTION_H

#include <exception>

class pException : public std::exception {
public:
  virtual const char* what() const throw() {
    return "Exception was thrown";
  }
} pexception;

#endif // __PLUG_EXCEPTION_H
