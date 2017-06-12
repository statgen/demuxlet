#include "sam_filtered_reader.h"

void SAMFilteredReader::init_params() {
  // check the sanity of input data
  if ( sam_file_name.empty() )
    error("[%s:%d %s] sam_file_name is empty", __FILE__, __LINE__, __PRETTY_FUNCTION__);

  if ( !target_interval_list.empty() ) {
    // stat processing the target region
    genomeLocus* tmp_target_locus = NULL;
    if ( !target_region.empty() ) {
      tmp_target_locus = new genomeLocus(target_region.c_str());
    }
    
    tsv_reader tsv_interval(target_interval_list.c_str());
    while( tsv_interval.read_line() ) {
      if ( tsv_interval.nfields < 3 )
	error("[E:%s:%d %s] Less than 3 columns observed in line %d of %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, tsv_interval.nlines, target_interval_list.c_str());
      const char* chrom = tsv_interval.str_field_at(0);
      int32_t beg1 = tsv_interval.int_field_at(1);
      int32_t end0 = tsv_interval.int_field_at(2);
      if ( target_region.empty() ) {
	target_loci.add(chrom, beg1, end0);
      }
      else if ( tmp_target_locus->overlaps(chrom, beg1, end0) ) {
	if ( tmp_target_locus->beg1 > beg1 )
	  beg1 = tmp_target_locus->beg1;
	if ( tmp_target_locus->end0 < end0 )
	  end0 = tmp_target_locus->end0;
	target_loci.add(chrom, beg1, end0);
      }
    }

    if ( tmp_target_locus )
      delete tmp_target_locus;
  }
  else if ( !target_region.empty() ) {
    target_loci.add(target_region.c_str());
  }

  // open the file handle
  file = hts_open(sam_file_name.c_str(), "r");
  if ( !file )
    error("[%s:%d %s] Cannot open %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sam_file_name.c_str());
  ftype = file->format.format;

  switch ( ftype ) {
  case sam: case bam: case cram:
    break;
  default:
    error("[%s:%d %s] File %s is not a SAM/BAM/CRAM file", __FILE__, __LINE__, __PRETTY_FUNCTION__, sam_file_name.c_str());
  }

  if ( !ref_file_name.empty() ) {
    char* fai = samfaipath(ref_file_name.c_str());
    hts_set_fai_filename(file, fai);
    free(fai);    
  }

  hdr = sam_hdr_read(file);
  if ( hdr == NULL )
    error("[%s:%d %s] Cannot read header from %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sam_file_name.c_str());    
  //s = bam_init1();

  if ( !target_loci.empty() ) {
    target_loci.rewind();
    while( !initialize_current_interval() ) {
      if ( target_loci.isend() ) {
	error("[E:%s] no informative interval found in\n", __PRETTY_FUNCTION__);	
      }
      target_loci.next();
    }    
  }
}

void SAMFilteredReader::set_buffer_size(int32_t new_buffer_size) {
  if ( new_buffer_size == 0 ) {
    if ( ridx + 1 != (int32_t)rbufs.size() )
      error("[E:%s:%d %s] Cannot set buffer size to be unlimited during limited buffer access access %d %u", __FILE__, __LINE__, __PRETTY_FUNCTION__, ridx, rbufs.size());
    unlimited_buffer = true;
  }
  else if ( new_buffer_size == (int32_t)rbufs.size() ) return;

  std::vector<bam1_t*> new_rbufs;  
  if ( new_buffer_size > (int32_t)rbufs.size() ) {
    // copy everything and add some empty ones at the beginning
    for(int32_t i=0; i < new_buffer_size-(int32_t)rbufs.size(); ++i)
      new_rbufs.push_back(bam_init1());
    for(int32_t i=0; i < (int32_t)rbufs.size(); ++i) {
      new_rbufs.push_back(rbufs[(ridx+1+i) % rbufs.size()]);      
    }
  }
  else {
    for(int32_t i=0; i < (int32_t)rbufs.size()-new_buffer_size; ++i)
      bam_destroy1(rbufs[(ridx+1+i) % rbufs.size()]);
    
    for(int32_t i=(int32_t)rbufs.size()-new_buffer_size; i < (int32_t)rbufs.size(); ++i)
      new_rbufs.push_back(rbufs[(ridx+1+i) % rbufs.size()]);
  }
  rbufs = new_rbufs;
  ridx = rbufs.size()-1;
}

bool SAMFilteredReader::load_index(bool continue_on_fail) {
  if ( idx == NULL ) {
    idx = sam_index_load(file, sam_file_name.c_str());
    if ( idx == NULL ) {
      if ( continue_on_fail ) {
	warning("[%s] Cannot load index", __PRETTY_FUNCTION__);
	return false;
      }
      else 
	error("[%s:%d %s] Cannot load index", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
    return idx ? true : false;
  }
  else return false;
}

bool SAMFilteredReader::jump_to(const char* chr, int32_t pos1) {
  int32_t cur_tid = cursor()->core.tid;
  int32_t cur_pos = cursor()->core.pos;
  int32_t new_tid = ( chr == NULL ) ? cur_tid : bam_name2id(hdr,chr);

  bool nojump = false;
  if ( ( cur_tid == new_tid ) && ( cur_pos + 1 < pos1 ) && ( pos1-cur_pos-1 < max_jumping_distance ) )
    nojump = true;

  if ( ignore_index || nojump ) {
    if ( ( cur_tid > new_tid ) || ( ( cur_tid == new_tid ) && ( cur_pos + 1 > pos1 ) ) )
      error("[E:%s Cannot go back to jump from %s:%d to %s:%d]", __PRETTY_FUNCTION__, bam_get_chrom(hdr,cursor()), cur_pos+1, chr, pos1);

    bam1_t* b = NULL;
    while( ( b = read() ) != NULL ) {
      if ( b->core.tid < new_tid ) continue;
      else if ( b->core.tid > new_tid ) break;
      else if ( b->core.pos + 1 < pos1 ) continue;
    }
    return ( b != NULL );
  }
  else {
    //int32_t end = INT_MAX;
    if ( !target_loci.empty() ) {
      if ( target_loci.moveTo(chr, pos1) ) { // move to specific position
	char buf[65536];
	sprintf(buf, "%s:%d-%d", chr, pos1, target_loci.it->end0);
	if ( itr ) hts_itr_destroy(itr);
	itr = sam_itr_querys(idx, hdr, buf);
	return itr ? true : false;
      }
      else {
	if ( target_loci.next() ) {
	  return initialize_current_interval(); // position is not within the target region, read the next target
	}
	else return false;
      }
    }
    else {
      if ( !load_index() )
	error("[E:%s] Cannot load index");
      if ( itr ) hts_itr_destroy(itr);
      char buf[65536];
      sprintf(buf, "%s:%d-%d", chr, pos1, INT_MAX);    
      itr = sam_itr_querys(idx, hdr, buf);
      return itr ? true : false;
    }
  }
}

bool SAMFilteredReader::initialize_current_interval() {
  if ( ignore_index ) return true;
  char buf[65536];
  load_index();
  sprintf(buf, "%s:%d-%d", target_loci.it->chrom.c_str(), target_loci.it->beg1, target_loci.it->end0);
  if ( itr ) hts_itr_destroy(itr);
  itr = sam_itr_querys(idx, hdr, buf);
  return itr ? true : false;
}

bam1_t* SAMFilteredReader::read() {
  // expand buffer if needed, and increase ridx by one
  if ( n_read % verbose == 0 ) {
    if ( cursor()->core.flag & 0x04 ) {
      notice("Reading %d reads (unmapped) and skipping %d",n_read,n_skip);      
    }
    else {
      notice("Reading %d reads at %s:%d and skipping %d",n_read,bam_get_chrom(hdr,cursor()),cursor()->core.pos+1,n_skip);
    }
  }
  
  if ( ( unlimited_buffer ) && ( nbuf == (int32_t)rbufs.size() ) ) {
    if ( ridx + 1 == nbuf ) {
      rbufs.push_back(bam_init1());
      ridx = rbufs.size() - 1;
    }
    else {
     std::vector<bam1_t*> new_rbufs;
     for(int32_t i=0; i < (int32_t)rbufs.size(); ++i) {
       new_rbufs.push_back(rbufs[(ridx+i) % rbufs.size()]);
     }
     new_rbufs.push_back(bam_init1());
     new_rbufs.swap(rbufs);
     ridx = rbufs.size() - 1;      
    }
  }
  else {
    ridx = (ridx + 1) % rbufs.size();    
  }
   
  if ( itr ) {
    while ( true ) {
      if ( itr && ( sam_itr_next(file,itr,rbufs[ridx]) >= 0 ) ) {
	++n_read;
	if ( passed_filter(hdr, rbufs[ridx]) ) {
	  ++nbuf;
	  return rbufs[ridx];
	}
	else {
	  ++n_skip;
	  continue;
	}
      }
      else if ( target_loci.next() ) { // if there are more target regions to read from
	initialize_current_interval();
      }
      else {
	ridx = (ridx + rbufs.size() - 1) % rbufs.size();
	eof = true;
	return NULL; }
    }
  }
  else {
    while( sam_read1(file, hdr, rbufs[ridx]) >= 0 ) {
      ++n_read;
      if ( passed_filter(hdr, rbufs[ridx]) ) {
	if ( target_loci.empty() ) {
	  ++nbuf;
	  return rbufs[ridx];
	}
	else {
	  if ( rbufs[ridx]->core.tid >= 0 ) {
	    const char* chr = bam_get_chromi(hdr, rbufs[ridx]->core.tid);
	    int32_t endpos = bam_endpos(rbufs[ridx]);
	    if ( target_loci.it->overlaps(chr, rbufs[ridx]->core.pos+1, endpos) ) {
	      ++nbuf;
	      return rbufs[ridx];
	    }
	    else if ( target_loci.overlaps(bam_get_chrom(hdr,rbufs[ridx]),rbufs[ridx]->core.pos+1, endpos) ) {
	      ++nbuf;
	      return rbufs[ridx];
	    }
	  }
	}
      }
      ++n_skip;
    }
    ridx = (ridx + rbufs.size() - 1) % rbufs.size();    
    eof = true;
    return NULL;
  }
}

int32_t SAMFilteredReader::clear_buffer_before(const char* chr, int32_t pos1) {
  int32_t tid;
  if ( chr == NULL ) {
    tid = rbufs[ridx]->core.tid;
    pos1 = rbufs[ridx]->core.pos+1;
  }
  else {
    tid = bam_name2id(hdr, chr);
    if ( tid < 0 ) error("[E:%s:%d %s] Cannot find chromosome %s from header", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr);    
  }
  int32_t n_rm = 0;
  for(int32_t i=0; i < nbuf; ++i) {
    bam1_t* b = rbufs[(ridx + nbuf - 1 -i) % rbufs.size()];
    if ( b->core.tid != tid ) ++n_rm;
    else if ( bam_endpos(b) < pos1 ) ++n_rm;
    else break;
  }
  nbuf -= n_rm;
  return n_rm;
}

bool SAMFilteredReader::passed_filter(bam_hdr_t* _hdr, bam1_t* b) {
  if ( filt.probThin < 1 )
    if ( rand()+0.5 > (RAND_MAX+1.0) * filt.probThin ) return false;
  
  if ( _hdr == NULL ) _hdr = hdr;
  if ( b == NULL ) b = cursor();

  if ( b->core.qual < filt.minMQ ) return false;
  //if ( filt.include_flag | b->core.flag ) {
  return ( filt.exclude_flag & b->core.flag ) ? false : true;
  //}
  //else return false;
}

void SAMFilteredReader::close() {
  if ( file ) { sam_close(file); file = NULL; }
  if ( hdr )  { bam_hdr_destroy(hdr); hdr = NULL; }
  if ( idx )  { hts_idx_destroy(idx); idx = NULL; }
  if ( itr )  { hts_itr_destroy(itr); itr = NULL; }
}
