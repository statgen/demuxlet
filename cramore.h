#ifndef __CRAMORE_H
#define __CRAMORE_H

#include <getopt.h>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <utility>
#include <string>
#include <set>

// Hyun's codes
#include "params.h"
#include "Error.h"
#include "PhredHelper.h"

// Adrian's codes
//#include "genome_interval.h"
//#include "bcf_ordered_reader.h"
//#include "bcf_chunked_reader.h"
//#include "bcf_ordered_writer.h"
#include "hts_utils.h"

extern "C" {
  size_t hts_realloc_or_die(unsigned long, unsigned long, unsigned long, unsigned long, int, void**, char const*);
}

// bcftools's code
//#include "filter.h"
//#include "genomeChunk.h"
//#include "genomeLoci.h"

// htslib's code

#endif
