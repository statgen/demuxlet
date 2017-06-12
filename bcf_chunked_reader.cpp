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

#include "bcf_chunked_reader.h"

BCFChunkedReader::BCFChunkedReader(const char* input_pattern_name, const char* input_ref_file, const char* input_interval_list, int32_t unit, genomeLoci* pIntervals) : chunk(input_pattern_name, input_ref_file, input_interval_list, unit, pIntervals) {
  init(pIntervals);
}

BCFChunkedReader::BCFChunkedReader(const char* input_vcf_file_name, genomeLoci* pIntervals) : chunk(input_vcf_file_name, pIntervals) {
  init(pIntervals);
}

BCFChunkedReader::BCFChunkedReader(const char* input_pattern_name, const char* input_ref_file, int32_t unit_chunk, genomeLoci* pIntervals) : chunk(input_pattern_name, input_ref_file, unit_chunk, pIntervals) {
  init(pIntervals);
}

BCFChunkedReader::BCFChunkedReader(const char* input_pattern_name, const char* input_interval_list, genomeLoci* pIntervals) : chunk(input_pattern_name, input_interval_list, pIntervals) {
  init(pIntervals);
}

void BCFChunkedReader::init(const char* input_pattern_name, const char* input_ref_file, const char* input_interval_list, int32_t unit, genomeLoci* pIntervals) {
  chunk.init(input_pattern_name, input_ref_file, input_interval_list, unit, pIntervals);
  init(pIntervals);
}


void BCFChunkedReader::init(genomeLoci* pIntervals) {
  //notice("[%s]", __PRETTY_FUNCTION__);
  file = NULL;
  hdr = NULL;
  idx = NULL;
  tbx = NULL;
  itr = NULL;
  
  if ( pIntervals ) {
    this->target_intervals = *pIntervals;
    this->target_intervals.resolveOverlaps();
  }
  
  //index_loaded = false;

  //intervals_present = this->target_intervals.empty() ? false : true;
  
  if ( !open_current_file() ) {
    error("[%s:%d %s] Cannot open any of the input file. Last file name is %s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
  }
}

bool BCFChunkedReader::jump_to(const char* chr, int32_t pos1) {
  //notice("[%s:%d %s] started", __FILE__, __LINE__, __PRETTY_FUNCTION__);
  std::pair<bool,bool> ret = chunk.jumpTo(chr,pos1);

  //notice("[%s:%d %s] jumped", __FILE__, __LINE__, __PRETTY_FUNCTION__);
  
  if ( !ret.first ) return false;

  if ( ret.second )
    open_current_file();

  //notice("[%s:%d %s] Finished open_current_file()", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    
  int32_t end = chunk.chromosome_is_chunked ? chunk.chunk_intervals.it->end0 : INT_MAX;
  if ( !target_intervals.empty() ) {
    if ( chunk.chunk_intervals.it->overlaps(*target_intervals.it) ) {
      if ( chunk.chunk_intervals.it->end0 > target_intervals.it->end0 )
	end = target_intervals.it->end0;
    }
  }

  //notice("[%s:%d %s] pos1 = %d, end = %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, pos1, end);

  if ( end < pos1 ) {
    return false;
  }
  
  //load_index();

  //notice("[%s:%d %s] Finished load_index()", __FILE__, __LINE__, __PRETTY_FUNCTION__);  

  char buf[65536];
  sprintf(buf,"%s:%d-%d",chr, pos1, end);
  load_index();
  if ( ftype.format == bcf ) {
    //notice("[%s:%d %s] Querying %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, buf);
    if ( itr ) hts_itr_destroy(itr);    
    itr = bcf_itr_querys(idx, hdr, buf);
    if (itr)
      return true;
  }
  else if (ftype.format==vcf && ftype.compression==bgzf) {
    if ( itr ) hts_itr_destroy(itr);    
    itr = tbx_itr_querys(tbx, buf);
    if (itr)
      return true;
  }
  return false;
}

bool BCFChunkedReader::load_index() {
  if (ftype.format==bcf) {
    if ( idx == NULL ) {
      if ( ( idx = bcf_index_load(chunk.current_file_name.c_str()) ) == NULL ) {
	warning("[E:%s] index cannot be loaded for %s for random access\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
	return false;
      }
    }
  }
  else if (ftype.format==vcf) {
    if (ftype.compression==bgzf) {
      if ( tbx == NULL ) {
	if ( (tbx = tbx_index_load(chunk.current_file_name.c_str())) == NULL ) {
	  warning("[E:%s] index cannot be loaded for %s for random access\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
	  return false;
	}
      }
    }
    else
      error("[E:%s] no random access support for non-BGZF VCF file: %s\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
  }
  return true;
}

bool BCFChunkedReader::open_current_file() {
  close(false);

  //notice("[%s], opening file %s",__PRETTY_FUNCTION__,chunk.current_file_name.c_str());
  // attempts to load the first file
  file = hts_open(chunk.current_file_name.c_str(), "r");
  while( file == NULL ) {
    notice("Skipping to read non-existent file %s", chunk.current_file_name.c_str());
    if ( chunk.setNextFileName() ) {
      file = hts_open(chunk.current_file_name.c_str(), "r");      
    }
    else {
      return false;
    }
  }
  ftype = file->format;

  if (ftype.format!=vcf && ftype.format!=bcf)
    error("[%s:%d %s] Not a VCF/BCF file: %s\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, chunk.current_file_name.c_str());

  s.l = s.m = 0; s.s = NULL;
  if (file==NULL) error("[%s:%d %s] Something went wrong. Open file hanle is lost", __FILE__, __LINE__, __PRETTY_FUNCTION__);
  if ( hdr ) bcf_hdr_destroy(hdr);
  hdr = bcf_alt_hdr_read(file);
  if (!hdr) error("[%s:%d %s]. Failed reading the BCF/VCF header from %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
  
  // do not load the index right now. Load only when needed.
  /*
  if (ftype.format==bcf) {
    if ( intervals_present ) {
      if ( ( idx = bcf_index_load(chunk.current_file_name.c_str()) ) != NULL )
	index_loaded = true;
      else
	error("[E:%s] index cannot be loaded for %s for random access\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
    }
  }
  else if (ftype.format==vcf) {
    if (ftype.compression==bgzf) {
      if ( intervals_present ) {
	if ((tbx = tbx_index_load(chunk.current_file_name.c_str())))
	  index_loaded = true;
	else
	  error("[E:%s] index cannot be loaded for %s for random access\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
      }
    }
    else {
      if (intervals_present)
	error("[E:%s] no random access support for non-BGZF VCF file: %s\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());
    }
  }
  */

  if ( !target_intervals.empty() ) {
    target_intervals.rewind();
    while( !initialize_current_interval() ) {
      if ( target_intervals.isend() ) {
	error("[E:%s] no overlapping interval found in : %s\n", __PRETTY_FUNCTION__, chunk.current_file_name.c_str());	
      }
      target_intervals.next();
    }
  }

  //notice("[%s], Finishing.. returning true",__PRETTY_FUNCTION__);
  return true;
  //random_access_enabled = intervals_present && index_loaded;
}

/**
 * Gets sequence name of a record.
 */
const char* BCFChunkedReader::get_seqname(bcf1_t *v)
{
    return bcf_get_chrom(hdr, v);
};

/**
 * Checks if index is loaded.
 */
//bool BCFChunkedReader::is_index_loaded()
//{
//    return index_loaded;
//}

/**
 * Gets bcf header.
 */
bcf_hdr_t* BCFChunkedReader::get_hdr()
{
    return hdr;
};

/**
 * Initialize next interval.
 * Returns false only if all intervals are accessed.
 */
bool BCFChunkedReader::initialize_current_interval() {
  if ( chunk.chunk_intervals.empty() || ( chunk.chunk_intervals.it->overlaps(*target_intervals.it) ) ) {
    char buf[65536];
    sprintf(buf,"%s:%d-%d",target_intervals.it->chrom.c_str(),
	    (chunk.chunk_intervals.empty() || chunk.chunk_intervals.it->beg1 < target_intervals.it->beg1 ) ? target_intervals.it->beg1 : chunk.chunk_intervals.it->beg1,
	    (chunk.chunk_intervals.empty() || chunk.chunk_intervals.it->end0 > target_intervals.it->end0 ) ? target_intervals.it->end0 : chunk.chunk_intervals.it->end0);
    load_index();
    if ( ftype.format == bcf ) {
      if ( itr ) hts_itr_destroy(itr);
      itr = bcf_itr_querys(idx, hdr, buf);
      if (itr)
	return true;
    }
    else if (ftype.format==vcf && ftype.compression==bgzf) {
      if ( itr ) hts_itr_destroy(itr);            
      itr = tbx_itr_querys(tbx, buf);
      if (itr)
	return true;
    }
  }
  return false;
};

/**
 * Reads next record, hides the random access of different regions from the user.
 */
bool BCFChunkedReader::read(bcf1_t *v) {
  //bcf_clear(v);
  if ( itr ) {
    while( true ) {
      if ( itr && ( ( ( ftype.format == bcf ) && ( bcf_itr_next(file,itr,v) >= 0 ) ) ||
		    ( ( ftype.format == vcf ) && ( tbx_itr_next(file,tbx,itr,&s) >= 0 ) ) ) ) {
	if ( ftype.format == vcf ) {
	  vcf_parse1(&s, hdr, v);
	}
	return true;
      }
      else {
	if ( target_intervals.next() ) {
	  initialize_current_interval();
	}
	else {
	  if ( chunk.setNextFileName() ) {
	    if ( open_current_file() ) {
	      // do nothing; current interval will be initialized
	    }
	    else return false;
	  }
	  else return false;
	}
      }
    }
  }
  else {
    if (bcf_read(file, hdr, v)==0)
      return true;
    else if ( chunk.setNextFileName() ) {
      if ( open_current_file() ) {
	return (read(v));
      }
      else
	return false;
    }
    else
      return false;
  }
}

/**
 * Closes the file.
 */
void BCFChunkedReader::close(bool destroy_hdr)
{
  if ( file ) {
    bcf_close(file);
    file = NULL;
  }
  if ( destroy_hdr && hdr ) {
    bcf_hdr_destroy(hdr);
    hdr = NULL;
  }
  if ( idx ) {
    hts_idx_destroy(idx);
    idx = NULL;
  }
  if ( tbx ) {
    tbx_destroy(tbx);
    tbx = NULL;
  }
  if ( itr ) {
    hts_itr_destroy(itr);
    itr = NULL;
  }  
}
