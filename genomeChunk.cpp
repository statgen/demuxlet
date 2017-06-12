#include "genomeChunk.h"

void genomeChunk::init(const char* patternOrFileName,
		       const char* refFile,
		       const char* intervalFile,
		       int32_t unit,
		       genomeLoci* pTarget) {
  // parse the pattern info
  //notice("[%s:%d %s] started",__FILE__,__LINE__,__FUNCTION__);
  
  separated_by_chromosome = false;
  chromosome_is_chunked = false;
  custom_chunk_used = false;
  is_eof = false;
  chunk_unit = INT_MAX;
  bool beg_used = false;
  bool end_used = false;
  std::string sPatternOrFileName(patternOrFileName);
  size_t i, j, k, l, best;
  int32_t type;
  i=0;
  while(true) {
    j = sPatternOrFileName.find(GENOME_CHUNK_SEP_CHR,i);
    k = sPatternOrFileName.find(GENOME_CHUNK_SEP_BEG,i);
    l = sPatternOrFileName.find(GENOME_CHUNK_SEP_END,i);
    if ( ( j == k ) && ( k == l ) && ( l == std::string::npos ) ) break;
    //!= std::string::npos );
    //   ( (k = sPatternOrFileName.find("__BEG__",i)) != std::string::npos ) ||
    //	   ( (l = sPatternOrFileName.find("__END__",i)) != std::string::npos ) ; ) {    
    if ( ( j < k ) && ( j < l ) ) {
      best = j;
      type = GENOME_CHUNK_TYPE_CHR;
      separated_by_chromosome = true;
    }
    else if ( ( k < j ) && ( k < l ) ) {
      best = k;
      type = GENOME_CHUNK_TYPE_BEG;
      beg_used = true;
    }
    else if ( ( l < j ) && ( l < k ) ) {
      best = l;
      type = GENOME_CHUNK_TYPE_END;
      end_used = true;
    }
    else {
      error("[%s:%d %s] Cannot parse %s at position (i,j,k,l,std::string::npos)=(%u,%u,%u,%u,%u) \n", __FILE__, __LINE__, __FUNCTION__, patternOrFileName, i, j, k, l, std::string::npos);
      best = 0;
    }

    v_substrs.push_back( sPatternOrFileName.substr(i, best-i) );
    v_types.push_back( type );
    
    i = best + GENOME_CHUNK_SEP_LEN;

    //notice("[%s:%d %s] i=%u, best=%u, j=%u, k=%u, l=%u, v_substrs.size()=%u",__FILE__,__LINE__,__FUNCTION__,i,best,j,k,l,v_substrs.size());
  }
  v_substrs.push_back( sPatternOrFileName.substr(i) );

  if ( beg_used != end_used )
    error("[%s:%d %s] Cannot parse %s because both %s and %s do not exist\n", __FILE__, __LINE__, __FUNCTION__, patternOrFileName, GENOME_CHUNK_SEP_BEG, GENOME_CHUNK_SEP_END);

  if ( intervalFile != NULL ) {
    htsFile* hp = hts_open(intervalFile, "r");
    if ( hp == NULL )
      error("[E:%s:%d %s] Cannot open file %s for reading",__FILE__,__LINE__,__FUNCTION__, intervalFile);
    kstring_t str = {0,0,0};
    int32_t lstr = 0;
    int32_t nfields = 0;
    int32_t* fields = NULL;
    // model list is assumed to have [INFO_KEY] [MODEL_FILE] [INFO_DESCRIPTION = INFO_KEY if empty]
    for( int32_t i=0; ( lstr = hts_getline(hp, KS_SEP_LINE, &str) ) >= 0; ++i ) {
      fields = ksplit(&str, 0, &nfields);
      if ( nfields < 3 )
	error("[E:%s:%d %s] Less than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, intervalFile);
      else if ( nfields > 3 )
	warning("[E:%s:%d %s] More than three columns observed in line %d of %s",__FILE__,__LINE__,__FUNCTION__, i+1, intervalFile);

      // check whether the interval overlaps with target, if specified;
      const char* chrom = &str.s[fields[0]];
      int32_t beg1 = atoi(&str.s[fields[1]]);
      int32_t end0 = atoi(&str.s[fields[2]]);

      if ( ( pTarget == NULL ) || ( pTarget->overlaps(chrom,beg1,end0) ) ) 
	chunk_intervals.add( &str.s[fields[0]], atoi(&str.s[fields[1]]), atoi(&str.s[fields[2]]) );
    }
    custom_chunk_used = true;
  }
  else {
    //notice("[%s:%d %s]",__FILE__,__LINE__,__FUNCTION__);
    
    if ( separated_by_chromosome ) {
      if ( !beg_used ) {
	unit = INT_MAX;
      }

      // parse FASTA file info
      std::string sRefFile(refFile);
      ReferenceSequence refSeq(sRefFile);
      int32_t nseq = refSeq.fetch_nseq();
      std::vector<std::string> seqnames(nseq);
      std::vector<int32_t> seqlens(nseq);
      for(int32_t i=0; i < nseq; ++i) {
	seqnames[i] = refSeq.fetch_iseq_name(i);
	seqlens[i] = refSeq.fetch_seq_len(seqnames[i].c_str());
      }

      //notice("[%s:%d %s] nseq = %d",__FILE__,__LINE__,__FUNCTION__, nseq);      

      for(int32_t i=0; i < nseq; ++i) {
	for(int32_t j=0; j < seqlens[i]; j += unit) {
	
	  //notice("[%s:%d %s] i=%d, j=%d",__FILE__,__LINE__,__FUNCTION__, i, j);
	  
	  int32_t end = j + unit;
	  if ( end > seqlens[i] )
	    end = seqlens[i];
	  else
	    chromosome_is_chunked = true;
	
	  if ( ( pTarget == NULL ) || ( pTarget->overlaps(seqnames[i].c_str(),i+1,end) ) ) 
	    chunk_intervals.add( seqnames[i].c_str(), j+1, end );
	}
      }

      chunk_unit = unit;
    }
    else {
      if ( beg_used ) {
	error("[E:%s:%d %s] %s is expected but not observed in %s",__FILE__,__LINE__,__FUNCTION__, GENOME_CHUNK_SEP_CHR, patternOrFileName);      
      }
      // no chunking
    }
  }
  setFirstFileName();
}

void genomeChunk::setFirstFileName() {
  chunk_intervals.rewind();
  setFileName();
}

bool genomeChunk::setNextFileName() {
  if ( is_eof )
    return false;
  else if ( chunk_intervals.empty() ) {
    is_eof = true;
    return false;
  }
  else if ( chunk_intervals.next() ) {
    setFileName();
    return true;
  }
  else {
    is_eof = true;
    return false;
  }
}

// first bool - jump succceed
// second bool - file changed
std::pair<bool,bool> genomeChunk::jumpTo(const char* chr, int32_t pos1) {
  if ( chunk_intervals.empty() )
    return std::make_pair<bool,bool>(true,false); // file name has not been changed


  std::set<genomeLocus>::iterator old_it = chunk_intervals.it;
  
  bool ret = chunk_intervals.moveTo(chr, pos1);
  if ( ret == false ) return std::make_pair<bool,bool>(false,false);

  std::set<genomeLocus>::iterator new_it = chunk_intervals.it;
  if ( old_it != new_it ) {
    setFileName();
    return std::make_pair<bool,bool>(true,true);
  }
  else {
    return std::make_pair<bool,bool>(true,false);    
  }
}

void genomeChunk::setFileName() {
  current_file_name = v_substrs[0];
  char buf[256];
  for(int32_t i=0; i < (int32_t)v_types.size(); ++i) {
    switch ( v_types[i] ) {
    case GENOME_CHUNK_TYPE_CHR:
      current_file_name += chunk_intervals.it->chrom;
      break;
    case GENOME_CHUNK_TYPE_BEG:
      sprintf(buf, "%d", chunk_intervals.it->beg1);
      current_file_name += buf;
      break;
    case GENOME_CHUNK_TYPE_END:
      sprintf(buf, "%d", chunk_intervals.it->end0);
      current_file_name += buf;
      break;
    default:
      error("[%s:%d %s] Cannot recognize v_types[%d]=%d\n", __FILE__, __LINE__, __FUNCTION__, i, v_types[i]);
    }
    current_file_name += v_substrs[i+1];
  }
}
