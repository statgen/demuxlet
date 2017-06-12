#ifndef __SAM_FILTERED_READER_H
#define __SAM_FILTERED_READER_H

#include <cstdlib>
#include <cstring>
#include <vector>

#include "bam_ordered_reader.h"
#include "Error.h"
#include "tsv_reader.h"
#include "genomeLoci.h"

struct _sam_filter_arg {
  uint32_t include_flag;
  uint32_t exclude_flag;
  int32_t minMQ;    
  double probThin;

  _sam_filter_arg() {
    include_flag = 0xffff;
    exclude_flag = 0x0000;
    minMQ = 0;
    probThin = 1.0;
  }
};

typedef struct _sam_filter_arg sam_filter_arg;

class SAMFilteredReader {
public:
  // arguments that can be controlled externally
  std::string sam_file_name;
  std::string ref_file_name;
  std::string target_region;
  std::string target_interval_list;

  sam_filter_arg filt;
  int32_t verbose;
  int32_t n_read;
  int32_t n_skip;

  int32_t ridx;
  std::vector<bam1_t*> rbufs;
  bool unlimited_buffer;
  int32_t nbuf;
  int32_t max_jumping_distance;
  bool ignore_index;
  bool eof;

  genomeLoci target_loci;

  samFile*   file;
  bam_hdr_t* hdr;
  hts_idx_t* idx;
  hts_itr_t* itr;
  htsExactFormat  ftype;
  kstring_t  str;

 SAMFilteredReader() : verbose(1000000), n_read(0), n_skip(0), ridx(-1), unlimited_buffer(false), nbuf(0), max_jumping_distance(0), ignore_index(false), eof(false), file(NULL), hdr(NULL), idx(NULL), itr(NULL), ftype(unknown_format) {
    str.l=str.m=0; str.s = NULL;
    set_buffer_size(1);
  }
  void init_params();
  void set_buffer_size(int32_t buffer_size);
  inline bam1_t* cursor() { return rbufs[ridx]; }
  bam1_t* cursor(int32_t idx) {
    if ( idx < nbuf )
      return rbufs[(ridx + rbufs.size() - idx) % rbufs.size()];
    else {
      error("[E:%s:%d %s] Cannot move cursor to %d-th variant, which is out of bound from %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, idx, nbuf);
      return 0;
    }
  }
  bam1_t* read();
  void pop() { --nbuf; }
  int32_t clear_buffer_before(const char* chr = NULL, int32_t pos = INT_MAX);

  bool initialize_current_interval();
  bool load_index(bool continue_on_fail = false);  
  
  bool passed_filter(bam_hdr_t* _hdr=NULL, bam1_t* b=NULL);
  bool jump_to(const char* chr = NULL, int32_t pos=INT_MAX);
  void close();
};

#endif
