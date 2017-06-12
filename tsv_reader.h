#ifndef __TSV_READER_H
#define __TSV_READER_H

#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include "htslib/kstring.h"
#include "htslib/kseq.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "Error.h"

class tsv_reader {
public:
  std::string filename;
  htsFile* hp;
  tbx_t* tbx;
  hts_itr_t* itr;
  kstring_t str;
  int32_t lstr;
  int32_t nfields;
  int32_t* fields;
  int32_t nlines;

  
  bool open(const char* filename);
  int32_t read_line();
  const char* str_field_at(int32_t idx);
  int32_t int_field_at(int32_t idx);
  double double_field_at(int32_t idx);
  int32_t store_to_vector(std::vector<std::string>& v);
  bool jump_to(const char* reg);
  bool jump_to(const char* chr, int32_t beg, int32_t end = INT_MAX);

 tsv_reader() : hp(NULL), tbx(NULL), itr(NULL), lstr(0), nfields(0), fields(NULL), nlines(0) {
    str.l = str.m = 0; str.s = NULL;
  }

 tsv_reader(const char* filename) : hp(NULL), tbx(NULL), itr(NULL), lstr(0), nfields(0), fields(NULL), nlines(0) {
    str.l = str.m = 0; str.s = NULL;    
    open(filename);
  }  

  ~tsv_reader() {
    if ( str.s ) free(str.s);
    if ( fields ) free(fields);
  }
};
#endif
