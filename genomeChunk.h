#ifndef __GENOME_CHUNK_H
#define __GENOME_CHUNK_H

#define GENOME_CHUNK_TYPE_CHR 1
#define GENOME_CHUNK_TYPE_BEG 2
#define GENOME_CHUNK_TYPE_END 3

#define GENOME_CHUNK_SEP_CHR "-_CHR_-"
#define GENOME_CHUNK_SEP_BEG "-_BEG_-"
#define GENOME_CHUNK_SEP_END "-_END_-"
#define GENOME_CHUNK_SEP_LEN 7

#include <cstdlib>
#include <cstring>
#include <climits>
#include <vector>
#include <algorithm>
#include "genomeLoci.h"
#include "Error.h"
#include "hts_utils.h"
#include "reference_sequence.h"

class genomeChunk {
 public:
  bool separated_by_chromosome;
  bool chromosome_is_chunked;
  bool custom_chunk_used;
  bool is_eof;
  int32_t chunk_unit;
  genomeLoci chunk_intervals;

  std::string current_file_name;
  std::vector<int32_t> v_types;
  std::vector<std::string> v_substrs;

  void init(const char* patternOrFileName, const char* refFile, const char* intervalFile,
	    int32_t unit, genomeLoci* pTarget);


  genomeChunk() {}
  genomeChunk(const char* patternName, const char* refName, const char* intervalFile, int32_t unit, genomeLoci* pTarget) {
    init(patternName, refName, intervalFile, unit, pTarget);
  }

  genomeChunk(const char* patternName, const char* refName, int32_t unit, genomeLoci* pTarget = NULL) {
    init(patternName, refName, NULL, unit, pTarget);
  }
  genomeChunk(const char* patternName, const char* intervalFile, genomeLoci* pTarget = NULL) {
    init(patternName, NULL, intervalFile, INT_MAX, pTarget);    
  }
  genomeChunk(const char* fileName, genomeLoci* pTarget = NULL) {
    init(fileName, NULL, NULL, INT_MAX, pTarget);
  }

  void setFirstFileName();
  bool setNextFileName();
  std::pair<bool,bool> jumpTo(const char* chr = NULL, int32_t pos1 = INT_MAX);
  void setFileName();
};
#endif
