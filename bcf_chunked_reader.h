/* The MIT License

   Copyright (c) 2013 Adrian Tan <atks@umich.edu>

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef BCF_CHUNKED_READER_H
#define BCF_CHUNKED_READER_H

#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "hts_utils.h"
#include "genome_interval.h"
#include "interval_tree.h"
#include "Error.h"
#include "genomeChunk.h"
#include "genomeLoci.h"

/**
 * A class for reading chunked VCF/BCF files.
 *
 * Basically a record iterator that hides the
 * htslib interface from the programs in vt.
 *
 * The main impetus for this is that htslib
 * is currently incorporated as a very early
 * stage and is thus lacking many feature that
 * is useful to us at this point in time. This
 * allows us to isolate the changes required in
 * future to simply methods in hts_utils.
 *
 * The following cases are supported.
 *
 * 1) Input is an unindexed file which is not necessarily ordered.
 * 2) Input is an indexed file
 *
 * This class hides the handling of indices from
 * the user and also allows for the selection of
 * records in intervals in both cases 1 and 2.
 */

class BCFChunkedReader
{
    public:

    ///////
    //i/o//
    ///////

    // information needed if pattern is used
    genomeChunk chunk;
    vcfFile *file;
    bcf_hdr_t *hdr;
    hts_idx_t *idx;
    tbx_t *tbx;
    hts_itr_t *itr;
    bcf1_t *v;
    //bcf1_t *buffer_v;

    //for control
    htsFormat ftype;
    //bool intervals_present;
    //bool index_loaded;
    //bool random_access_enabled;

    //list of intervals
    genomeLoci target_intervals;

    //for storing unused bcf records
    //std::list<bcf1_t*> pool;

    //shared objects for string manipulation
    kstring_t s;

    /**
     * Initialize files and intervals.
     *
     * @input_vcf_file_name     name of the input VCF file
     * @intervals          list of intervals, if empty, all records are selected.
     */
    void init(genomeLoci* pIntervals);

    void init(const char* input_pattern_name, const char* input_ref_file, const char* input_interval_list, int32_t unit, genomeLoci* pIntervals);

    bool open_current_file();

    BCFChunkedReader() : file(NULL), hdr(NULL), idx(NULL), tbx(NULL), itr(NULL), v(NULL) {}

    BCFChunkedReader(const char* input_vcf_file_name, genomeLoci* pIntervals = NULL);

    BCFChunkedReader(const char* input_pattern_name, const char* input_ref_file, int32_t unit_chunk, genomeLoci* pIntervals = NULL);

    BCFChunkedReader(const char* input_pattern_name, const char* input_interval_list, genomeLoci* pIntervals = NULL);

    BCFChunkedReader(const char* input_pattern_name, const char* input_ref_file, const char* input_interval_list, int32_t unit, genomeLoci* pIntervals);
    
    
    /**
     * Returns next vcf record.
     */
    bool read(bcf1_t *v);

    /**
    * Initialize next interval.
    * Returns false only if all intervals are accessed.
    */
    bool initialize_current_interval();

//    /**
//     * Returns next set of vcf records at a start position.
//     * Note that this function should never be used in conjunction with read(bcf1_t *v)
//     */
//    bool read_next_position(std::vector<bcf1_t *>& vs);

    /**
     * Gets sequence name of a record.
     */
    const char* get_seqname(bcf1_t *v);

    /**
     * Checks if index is loaded.
     */
    //bool is_index_loaded();
    
    /**
     * Gets bcf header.
     */
    bcf_hdr_t* get_hdr();

    /**
     * Closes the file.
     */
    void close(bool destroy_hdr = true);

    bool load_index();

    bool jump_to(const char* chr = NULL, int32_t pos1 = INT_MAX);
};

#endif
