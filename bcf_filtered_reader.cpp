#include "bcf_filtered_reader.h"

void BCFFilteredReader::init_params() {
  if ( bcf_file_name.empty() )
    error("[%s:%d %s] bcf_file_name is empty", __FILE__, __LINE__, __PRETTY_FUNCTION__);
    
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

    if ( tmp_target_locus)
      delete tmp_target_locus;
  }
  else if ( !target_region.empty() ) {
    target_loci.add(target_region.c_str());
  }
  // end processing the target region

  // read sex map
  if ( !sexMap.empty() ) {
    tsv_reader tsv_sex_map(sexMap.c_str());
    int32_t ncols = 0, icol = -1, scol;
    for(int32_t i=0; ; ++i) {
      if ( ncols == 0 ) {
	ncols = tsv_sex_map.read_line();
	if ( ncols == 2 ) { icol = 1; scol = 0; }
	else if ( ncols > 4 ) { icol = 4; scol = 1; }
	else error("[E:%s:%d %s] Cannot read sex map %s because of the column size = %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, sexMap.c_str(), ncols);
      }
      else if ( ncols != tsv_sex_map.read_line() ) {
	if ( tsv_sex_map.nfields <= 0 ) break;
	else error("[E:%s:%d %s] Sex map file has inconsistent number of columns at line %d.", __FILE__, __LINE__, __PRETTY_FUNCTION__, i);
      }
      int32_t sex = tsv_sex_map.int_field_at(icol);
      std::string id = tsv_sex_map.str_field_at(scol);
      if ( mSex.find(id) != mSex.end() )
	error("[E:%s:%d %s] Duplicate ID %s in %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, id.c_str(), sexMap.c_str());
      if ( sex == 0 ) {
	warning("Unknown sex for individual %s, assuming female", id.c_str());
	sex = 2;
      }
      else if ( ( sex < 0 ) || ( sex > 2 ) ) {
	error("[E:%s:%d %s] Unknown sex %d for %s in %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, sex, id.c_str(), sexMap.c_str());	
      }
      mSex[id] = (sex == 1 ? 1 : 0);
    }
  }


  ////////////////////////////////////////
  // Steps to initialize cdr object
  ////////////////////////////////////////
  if ( bcf_file_name.find(GENOME_CHUNK_SEP_CHR) == std::string::npos ) {
    // chunk is not necessary
    if ( ( !ref_file_name.empty() ) || ( !interval_file_name.empty() ) ) 
      warning("[W:%s:%d %s] ref_file_name or interval_file_name is used but will be ignored because chunking is not necessary. Use %s in the filename if you intended chunking the file", __FILE__, __LINE__, __PRETTY_FUNCTION__, GENOME_CHUNK_SEP_CHR);
    cdr.init(bcf_file_name.c_str(), NULL, NULL, INT_MAX, target_loci.empty() ? NULL : &target_loci);
  }
  else if ( !ref_file_name.empty() ) {
    if ( !interval_file_name.empty() ) 
      warning("[W:%s:%d %s] interval_file_name is used but will be ignored because ref_file_name is provided. Do not use ref_file_name if you intended chunking by the interval", __FILE__, __LINE__, __PRETTY_FUNCTION__, GENOME_CHUNK_SEP_CHR);
    
    cdr.init(bcf_file_name.c_str(), ref_file_name.c_str(), NULL, unit, target_loci.empty() ? NULL : &target_loci);
  }
  else if ( !interval_file_name.empty() ) {
    cdr.init(bcf_file_name.c_str(), NULL, interval_file_name.c_str(), INT_MAX, target_loci.empty() ? NULL : &target_loci);    
  }

  // default buffer size is 1
  set_buffer_size(1);

  // initialize filter
  vfilt.init(cdr.hdr);

  //notice("Processing sample info");
  if ( !sample_id_list.empty() ) {
    tsv_reader tsv_sample(sample_id_list.c_str());
    while( tsv_sample.read_line() > 0 ) {
      sm_ids.insert(tsv_sample.str_field_at(0));
    }
    notice("Finished loading %u IDs from %s",sm_ids.size(), sample_id_list.c_str());
  }

  // process sample info
  if ( !sm_ids.empty() ) {
    std::set<std::string>::iterator it = sm_ids.begin();
    while( it != sm_ids.end() ) {
      int32_t idx = bcf_hdr_sample_index(cdr.hdr, it->c_str());
      if ( idx < 0 )
	error("[E:%s:%d %s] Cannot find sample ID %s from the BCF file", __FILE__, __LINE__, __PRETTY_FUNCTION__, it->c_str());
      sm_icols.push_back(idx);
      v_sm_ids.push_back(bcf_hdr_sample_id(cdr.hdr,idx));
      if ( !mSex.empty() ) {
	if ( mSex.find(*it) == mSex.end() )
	  error("[E:%s:%d %s] Cannot find sex information of %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, it->c_str());
	sm_isexes.push_back(mSex[*it]);
      }
      else {
	sm_isexes.push_back(0);
      }
      ++it;
    }
  }
  else {
    int32_t nsamples = bcf_hdr_nsamples(cdr.hdr);
    for(int32_t i=0; i < nsamples; ++i) {
      sm_icols.push_back(i);
      v_sm_ids.push_back(bcf_hdr_sample_id(cdr.hdr,i));      
      if ( !mSex.empty() ) {
	std::string id(bcf_hdr_sample_id(cdr.hdr, i));
	if ( mSex.find(id) == mSex.end() )
	  error("[E:%s:%d %s] Cannot find sex information of %s", __FILE__, __LINE__, __PRETTY_FUNCTION__, id.c_str());
	sm_isexes.push_back(mSex[id]);
      }
      else {
	sm_isexes.push_back(0);	
      }
    }
  }

  notice("Finished identifying %u samples to load from VCF/BCF",sm_icols.size());
  if ( sm_icols.empty() ) {
    error("[E:%s] No sample to load from VCF/BCF",__PRETTY_FUNCTION__);
  }

  // set the rids for chrX, chrY, chrMT
  xRid  = bcf_hdr_name2id(cdr.hdr, xLabel.c_str());
  yRid  = bcf_hdr_name2id(cdr.hdr, yLabel.c_str());
  mtRid = bcf_hdr_name2id(cdr.hdr, mtLabel.c_str());
  ploidies = new int8_t[sm_isexes.size()];
  memset(ploidies, (int8_t)2, (int32_t)sm_isexes.size());

  sex_ploidies[0] = sex_ploidies[1] = 2;
}

void BCFFilteredReader::set_buffer_size(int32_t new_buffer_size) {
  if ( new_buffer_size == 0 ) { // keep the current buffer size
    if ( vidx + 1 != (int32_t)vbufs.size() )
      error("[E:%s:%d %s] Cannot set buffer size to be unlimited during limited buffer access access %d %u", __FILE__, __LINE__, __PRETTY_FUNCTION__, vidx, vbufs.size());
    unlimited_buffer = true;
  }
  else if ( new_buffer_size == (int32_t)vbufs.size() ) return;

  std::vector<bcf1_t*> new_vbufs;  
  if ( new_buffer_size > (int32_t)vbufs.size() ) {
    // copy everything and add some empty ones at the beginning
    for(int32_t i=0; i < new_buffer_size-(int32_t)vbufs.size(); ++i)
      new_vbufs.push_back(bcf_init());
    for(int32_t i=0; i < (int32_t)vbufs.size(); ++i) {
      new_vbufs.push_back(vbufs[(vidx+1+i) % vbufs.size()]);      
    }
  }
  else {
    for(int32_t i=0; i < (int32_t)vbufs.size()-new_buffer_size; ++i)
      bcf_destroy(vbufs[(vidx+1+i) % vbufs.size()]);
    
    for(int32_t i=(int32_t)vbufs.size()-new_buffer_size; i < (int32_t)vbufs.size(); ++i)
      new_vbufs.push_back(vbufs[(vidx+1+i) % vbufs.size()]);
  }
  vbufs.swap(new_vbufs);
  vidx = vbufs.size()-1;
}

bool BCFFilteredReader::parse_genotypes(bcf_hdr_t* hdr, bcf1_t* v) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  bcf_unpack(v, BCF_UN_ALL);

  //notice("Parsing genotypes at %s:%d:%s:%s", bcf_hdr_id2name(hdr,v->rid), v->pos+1, v->d.allele[0], v->d.allele[1]);

  if ( bcf_get_genotypes(hdr, v, &gts, &n_gts) < 0 ) {
    notice("Failed to parsing genotypes at %s:%d:%s:%s", bcf_hdr_id2name(hdr,v->rid), v->pos+1, v->d.allele[0], v->d.allele[1]);    
    return false;
  }
  int32_t nsamples = bcf_hdr_nsamples(cdr.hdr);
  acs.resize(v->n_allele);
  std::fill(acs.begin(), acs.end(), 0);
  an = 0;
  
  if ( gfilt.minDP > 0 ) {
    if ( !parse_int_fields("DP",hdr,v) ) {
	error("[E:%s:%d %s] Cannot find the field DP from the VCF file at position %s:%d",__FILE__,__LINE__,__PRETTY_FUNCTION__, bcf_hdr_id2name(hdr, v->rid), v->pos+1);      
    }
    const int32_t* iflds = (const int32_t*)flds;
    for(int32_t i=0; i < nsamples; ++i) {
      if ( iflds[i] < gfilt.minDP ) {
	gts[2*i] = bcf_gt_missing;
	gts[2*i+1] = bcf_gt_missing;	    
      }
    }      
  }

  if ( gfilt.minGQ > 0 ) {
    if ( !parse_int_fields("GQ",hdr,v) ) {
	error("[E:%s:%d %s] Cannot find the field DP from the VCF file at position %s:%d",__FILE__,__LINE__,__PRETTY_FUNCTION__, bcf_hdr_id2name(hdr, v->rid), v->pos+1);      
    }
    const int32_t* iflds = (const int32_t*)flds;
    for(int32_t i=0; i < nsamples; ++i) {
      if ( iflds[i] < gfilt.minGQ ) {
	gts[2*i] = bcf_gt_missing;
	gts[2*i+1] = bcf_gt_missing;	    
      }
    }      
  }

  // calculate allele count and allele balance for specified individuals
  int32_t tmp_gt;
  for(int32_t i=0; i < (int32_t)sm_icols.size(); ++i) {
    for(int32_t j=0; j < 2; ++j) {
      tmp_gt = gts[sm_icols[i]*2+j];
      if ( ( j < ploidies[i] ) && ( bcf_gt_allele(tmp_gt) >= 0 ) ) {
	//notice("%d %d %x %x", i, j, tmp_gt, bcf_gt_allele(tmp_gt));
	++an;
	++acs[bcf_gt_allele(tmp_gt)];
      }
    }
  }
  return true;
}

bool BCFFilteredReader::parse_likelihoods(bcf_hdr_t* hdr, bcf1_t* v, const char* name) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  bcf_unpack(v, BCF_UN_ALL);

  if ( bcf_get_format_int32(hdr, v, name, &pls, &n_pls) < 0 ) {
    return false;
  }

  // E-M frequency estimation
  int32_t niter = 10;
  int32_t nalleles = cursor()->n_allele;
  acs.resize(nalleles);
  for(int32_t i=0; i < nalleles; ++i)
    acs[i] = 1.0/nalleles; // acs represents initial AF
  int32_t ngenos = (nalleles+1)*nalleles/2;
  int32_t nsamples = bcf_hdr_nsamples(cdr.hdr);  
  
  gps = (float*) realloc(gps, sizeof(float)*ngenos*nsamples);
  n_gps = ngenos * nsamples; //sm_icols.size();
    
  double* gp = new double[ngenos];
  double sumgp;
  int32_t icol, i, j, k, l;

  for(int32_t it=0; it < niter; ++it) {
    std::vector<double> newacs(nalleles,0);
    an = 0;
    for(i=0; i < (int32_t)sm_icols.size(); ++i) {
      icol = sm_icols[i]*ngenos;
      if ( ploidies[i] == 2 ) {
	sumgp = 0;
	for(j=0, l=0; j < nalleles; ++j) {
	  for(k=0; k <= j; ++k, ++l) {
	    sumgp += (gp[l] = (j == k ? 1 : 2) * acs[j] * acs[k] * phredConv.toProb(pls[icol+l]));
	  }
	}
	for(j=0, l=0; j < nalleles; ++j) {
	  for(k=0; k <= j; ++k, ++l) {
	    gp[l] /= sumgp;
	    newacs[j] += gp[l];
	    newacs[k] += gp[l];
	  }
	}
	an += 2;
      }
      else if ( ploidies[i] == 1 ) {
	memset(gp, 0, ngenos*sizeof(double));
	sumgp = 0;
	for(j=0; j < nalleles; ++j) {
	  l = (j+1)*(j+2)/2-1;
	  sumgp += (gp[l] = acs[j] * phredConv.toProb(pls[icol+l]));
	}
	for(j=0, l=0; j < nalleles; ++j) {
	  l = (j+1)*(j+2)/2-1;
	  gp[l] /= sumgp;
	  newacs[j] += gp[l];
	}	
	++an;
      }
      if ( it + 1 == niter ) {
	for(l=0; l < ngenos; ++l) 
	  gps[icol+l] = (float)gp[l];
      }
    }
    for(i=0; i < nalleles; ++i)
      acs[i] = newacs[i]/an;
  }

  for(i=0; i < nalleles; ++i)
    acs[i] *= an;

  delete[] gp;
  
  return true;
}

bool BCFFilteredReader::parse_dosages(bcf_hdr_t* hdr, bcf1_t* v, const char* name) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  bcf_unpack(v, BCF_UN_ALL);

  if ( bcf_get_format_float(hdr, v, name, &dss, &n_dss) < 0 ) {
    return false;
  }

  // need to calculate genotype dosages and allele frequencies, but let's skip for now
  int32_t i, j;
  float f;
  bool missing;
  int32_t nalleles = cursor()->n_allele;
  for(int32_t i=0; i < nalleles; ++i)
    acs[i] = 0;
  an = 0;
  for(i=0; i < (int32_t)sm_icols.size(); ++i) {
    missing = false;
    for(j=0; j < nalleles-1; ++j) {
      f = dss[i*(nalleles-1)+j];
      if ( bcf_float_is_vector_end(f) || bcf_float_is_missing(f) ) {
	missing = true;
	// nothing
      }
      else {
	acs[j+1] += f;
      }
    }
    if ( !missing ) an += 2;
  }
  acs[0] = an;
  for(j=1; j < nalleles; ++j) acs[0] -= acs[j];
  
  return true;
}

bool BCFFilteredReader::parse_posteriors(bcf_hdr_t* hdr, bcf1_t* v, const char* name, double gt_error) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  //notice("Parsing posteriors at %s:%d:%s:%s, name=%s", bcf_hdr_id2name(hdr,v->rid), v->pos+1, v->d.allele[0], v->d.allele[1],name);  

  bcf_unpack(v, BCF_UN_ALL);

  if ( strcmp(name,"GT") == 0 ) {
    int32_t i, j, k, l, g, icol;
    int32_t nalleles = v->n_allele;
    int32_t ngenos = (nalleles+1)*nalleles/2;
    int32_t nsamples = bcf_hdr_nsamples(cdr.hdr);    
    //notice("Before parsing genotypes at %s:%d:%s:%s", bcf_hdr_id2name(hdr,v->rid), v->pos+1, v->d.allele[0], v->d.allele[1]);      
    if ( !parse_genotypes(hdr, v) ) return false;
    //notice("After parsing genotypes at %s:%d:%s:%s", bcf_hdr_id2name(hdr,v->rid), v->pos+1, v->d.allele[0], v->d.allele[1]);          
    gps = (float*) realloc(gps, sizeof(float)*ngenos*nsamples);
    n_gps = nsamples * 3;    
    for(i=0; i < (int32_t)sm_icols.size(); ++i) {
      g = get_genotype_at(i);
      icol = sm_icols[i]*ngenos;      
      if ( g < 0 ) { // missing genotype
	if ( ploidies[i] == 2 ) {
	  for(j=0, l=0; j < nalleles; ++j) {
	    for(k=0; k <= j; ++k, ++l) {
	      gps[icol+l] = (float)((j == k ? 1.0 : 2.0) * (acs[j]+1.0/nalleles)/(an+1.0) * (acs[k]+1.0/nalleles)/(an+1.0));
	    }
	  }
	}
	else if ( ploidies[i] == 1 ) {
	  memset(&gps[icol], 0, ngenos*sizeof(double));
	  for(j=0; j < nalleles; ++j) {
	    l = (j+1)*(j+2)/2-1;
	    gps[icol+l] = (float)((acs[j]+1.0/nalleles)/(an+1.0));
	  }
	}
      }
      else {
	for(j=0; j < ngenos; ++j) {
	  gps[icol+j] = (g == j) ? 1.0-gt_error : gt_error/(ngenos-1.0);
	}
      }
    }
    //if ( v->pos % 1000 == 0 )
    //notice("** %f %f %f %f %f %f",gps[3],gps[4],gps[5],gps[6],gps[7],gps[8]);
    return true;
  }
  else if ( strcmp(name,"PL") == 0 ) {
    return parse_likelihoods(hdr, v, name);
  }
  else { // GP as posterior
    if ( bcf_get_format_float(hdr, v, name, &gps, &n_gps) < 0 ) {
      return false;
    }

    int32_t i, j, icol;
    float sumgp;
    int32_t nalleles = v->n_allele;
    int32_t ngenos = (nalleles+1)*nalleles/2;
    std::vector<float> gpSums(ngenos,0);

    for(i=0; i < nalleles; ++i) {
      for(j=0; j <= i; ++j) {
	gpSums[(i+1)*i/2 + j] = ((i == j) ? 1.0 : 2.0)/(float)(nalleles*nalleles);
      }
    }

    // first, normalize by GP calculate the sum of genotype probabilities
    for(i=0; i < (int32_t)sm_icols.size(); ++i) {
      icol = sm_icols[i]*ngenos;
      sumgp = 0;
      for(j=0; j < ngenos; ++j) {
	//gps[icol+j] += gt_error;
	sumgp += gps[icol+j];
      }
      for(j=0; j < ngenos; ++j) {
	gps[icol+j] /= sumgp;
	gpSums[j] += gps[icol+j];
      }
    }

    for(j=0; j < ngenos; ++j)
      gpSums[j] /= (int32_t)(sm_icols.size()+1.0);

    // second, account for genotyping error as weighted sum of GP and gpSums
    for(i=0; i < (int32_t)sm_icols.size(); ++i) {
      icol = sm_icols[i]*ngenos;    
      for(j=0; j < ngenos; ++j) {
	gps[icol+j] = ((1.0-gt_error)*gps[icol+j] + gt_error*gpSums[j]);
      }
    }
    
    return true;
  }
}



bool BCFFilteredReader::parse_int_fields(const char* name, bcf_hdr_t* hdr, bcf1_t* v) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  bcf_unpack(v, BCF_UN_ALL);
  
  if ( bcf_get_format_int32(hdr, v, name, &flds, &n_flds) < 0 ) {  
    return false;
  }
  return true;  
}

bool BCFFilteredReader::parse_float_fields(const char* name, bcf_hdr_t* hdr, bcf1_t* v) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  bcf_unpack(v, BCF_UN_ALL);
  
  if ( bcf_get_format_float(hdr, v, name, &flds, &n_flds) < 0 ) {  
    return false;
  }
  return true;  
}

std::string& BCFFilteredReader::get_var_ID(bcf_hdr_t* hdr, bcf1_t* v) {
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();
  bcf_unpack(v, BCF_UN_STR);
  varID.assign(bcf_get_chrom(cdr.hdr,v));
  varID += ":";
  char buf[256];
  sprintf(buf,"%d",v->pos);
  varID += buf;
  for(int32_t i=0; i < v->n_allele; ++i) {
    varID += (i == 0 ? ":" : "_");
    varID += v->d.allele[i];
  }
  return varID;
}

bool BCFFilteredReader::passed_vfilter(bcf_hdr_t* hdr, bcf1_t* v) {
  if ( vfilt.probThin < 1 ) 
    if ( rand()+0.5 > (RAND_MAX+1.0) * vfilt.probThin ) return false;
  
  if ( hdr == NULL ) hdr = cdr.hdr;
  if ( v == NULL ) v = cursor();

  bool has_filter = vfilt.req_flt_ids.empty() ? true : false;

  if ( ( has_filter ) && ( vfilt.filt == NULL ) && ( !vfilt.require_GT) ) 
    return true;
  
  bcf_unpack(v, BCF_UN_FLT);
  if ( !has_filter ) {
    for(int32_t i=0; i < v->d.n_flt; ++i) {
      for(int32_t j=0; j < (int32_t)vfilt.req_flt_ids.size(); ++j) {
	if ( vfilt.req_flt_ids[j] == v->d.flt[i] )
	  has_filter = true;
      }
    }
  }

  if ( !has_filter ) return false;

  if ( vfilt.filt != NULL ) {
    int32_t ret = filter_test(vfilt.filt, v, NULL);
    if ( vfilt.filter_logic == FLT_INCLUDE ) {
      if (!ret) has_filter = false;
    }
    else if ( ret ) { has_filter = false; }
  }

  if ( !has_filter ) return false;

  if ( v->qual < vfilt.minQual ) return false;

  if ( v->n_allele > vfilt.maxAlleles ) return false;

  if ( v->n_allele < vfilt.minAlleles ) return false;

  if ( vfilt.snpOnly && (!bcf_is_snp(v)) ) return false;

  if ( vfilt.require_GT ) {
    //int32_t nsamples = bcf_hdr_nsamples(cdr.hdr);

    acs.resize(v->n_allele);
    std::fill(acs.begin(), acs.end(), 0);
    an = 0;

    //if ( rand() % 1000 == 0 ) notice("foo");        

    if ( !parse_genotypes(hdr,v) ) 
      error("[E:%s:%d %s] Cannot find the field GT from the VCF file at position %s:%d",__FILE__,__LINE__,__PRETTY_FUNCTION__, bcf_hdr_id2name(hdr, v->rid), v->pos+1);

    //if ( rand() % 1000 == 0 ) notice("foo");    

    if ( vfilt.minCallRate > (double)an/(2.0*(double)sm_icols.size()) ) return false;

    int32_t ac = an-acs[0];
    double af = (ac+5e-11)/(an+1e-10);

    //if ( rand() % 1000 ) notice("%d %d %lf",ac,an,af);

    //if ( rand() % 1000 == 0 ) notice("foo");
  
    if ( ac < vfilt.minAC ) return false;
    if ( ac > vfilt.maxAC ) return false;    
    if ( ( ac < vfilt.minMAC ) || ( an-ac < vfilt.minMAC ) ) return false;
    if ( ( ac > vfilt.maxMAC ) && ( an-ac > vfilt.maxMAC ) ) return false;    
    if ( af < vfilt.minAF ) return false;
    if ( ( af < vfilt.minMAF ) || ( 1.-af < vfilt.minMAF ) ) return false;
    if ( af > vfilt.maxAF ) return false;
    if ( ( af > vfilt.maxMAF ) && ( 1.-af > vfilt.maxMAF ) ) return false;
    return true;
  }
  else return true;
}

bool BCFFilteredReader::add_variant_to_extract(bcf_hdr_t* hdr, bcf1_t* v) {
  mode_extract = true;
  return variants2extract.insert(variantKeyS(hdr,v)).second;
}

bool BCFFilteredReader::jump_to(const char* chr, int32_t pos) {
  int32_t cur_rid, cur_pos, new_rid;
  cur_rid = cursor()->rid;
  cur_pos = cursor()->pos;
  new_rid = ( chr == NULL ) ? cur_rid : bcf_hdr_name2id(cdr.hdr, chr);
  if ( new_rid < 0 )
    error("[E:%s] Cannot move to %s", __PRETTY_FUNCTION__, chr);

  if ( ( cur_rid > new_rid ) || ( ( cur_rid == new_rid ) && ( cur_pos + 1 > pos ) ) )
    return cdr.jump_to(chr,pos);
  else if ( ( max_jumping_distance == INT_MAX) || ( ( cur_rid == new_rid ) && ( pos - 1 - cur_pos < max_jumping_distance ) ) ) {
    bcf1_t* v = NULL;
    while ( ( v = read() ) != NULL ) {
      if ( v->rid > new_rid ) return true;
      else if ( v->rid < new_rid ) continue;
      else if ( v->pos + v->rlen < pos ) continue;
      else return true;
    }
    return (v != NULL); 
  }
  else
    return cdr.jump_to(chr,pos);    
}

// set ploidies using sex. return true if ploidies are changed
bool BCFFilteredReader::set_ploidies_by_sex(bcf1_t* v) {
  if ( v == NULL ) v = cursor();
  if ( v->rid == xRid ) {
    if ( ( v->pos+1 >= xStart ) && ( v->pos+1 <= xStop) ) {
      if ( ( sex_ploidies[0] == 2 ) && ( sex_ploidies[1] == 1 ) ) return false;
      sex_ploidies[0] = 2;
      sex_ploidies[1] = 1;
      for(int32_t i=0; i < (int32_t)sm_isexes.size(); ++i)
	ploidies[i] = sex_ploidies[sm_isexes[i]];
      return true;
    }
    else {
      if ( ( sex_ploidies[0] == 2 ) && ( sex_ploidies[1] == 2 ) ) return false;
      sex_ploidies[0] = 2;
      sex_ploidies[1] = 2;      
      memset(ploidies, (int8_t)2, (int32_t)sm_isexes.size());
      return true;
    }
  }
  else if ( v->rid == yRid ) {
    if ( ( sex_ploidies[0] == 0 ) && ( sex_ploidies[1] == 1 ) ) return false;
    sex_ploidies[0] = 0;
    sex_ploidies[1] = 1;
    for(int32_t i=0; i < (int32_t)sm_isexes.size(); ++i)
      ploidies[i] = sex_ploidies[sm_isexes[i]];
    return true;    
  }
  else if ( v->rid == mtRid ) {
    if ( ( sex_ploidies[0] == 1 ) && ( sex_ploidies[1] == 1 ) ) return false;
    sex_ploidies[0] = 1;
    sex_ploidies[1] = 1;
    memset(ploidies, (int8_t)1, (int32_t)sm_isexes.size());    
    return true;
  }
  else {
    if ( ( sex_ploidies[0] == 2 ) && ( sex_ploidies[1] == 2 ) ) return false;    
    sex_ploidies[0] = 2;
    sex_ploidies[1] = 2;      
    memset(ploidies, (int8_t)2, (int32_t)sm_isexes.size());
    return true;    
  }
}

int32_t BCFFilteredReader::clear_buffer_before(const char* chr, int32_t pos1) {
  int32_t rid;
  if ( chr == NULL ) {
    rid = vbufs[vidx]->rid;
    pos1 = vbufs[vidx]->pos+1;
  }
  else {
    rid = bcf_hdr_name2id(cdr.hdr, chr);
    if ( rid < 0 ) error("[E:%s:%d %s] Cannot find chromosome %s from header", __FILE__, __LINE__, __PRETTY_FUNCTION__, chr);
  }

  int32_t n_rm = 0;
  for(int32_t i=0; i < nbuf; ++i) {
    bcf1_t* v = vbufs[(vidx + vbufs.size() - (nbuf-i-1) ) % vbufs.size()];
    if ( v->rid < rid ) ++n_rm;
    else if ( ( v->rid == rid ) && ( v->pos+v->rlen < pos1 ) ) ++n_rm;
    else break;
  }
  nbuf -= n_rm;
  return n_rm;
}
  
bcf1_t* BCFFilteredReader::read() {
  n_gts = 0; n_pls = 0; n_dss = 0; n_flds = 0;
  
  if ( ( unlimited_buffer ) && ( nbuf == (int32_t)vbufs.size() ) ) {
    if ( vidx + 1 == nbuf ) {
      vbufs.push_back(bcf_init());
      vidx = vbufs.size()-1;
    }
    else {
      std::vector<bcf1_t*> new_vbufs;
      for(int32_t i=0; i < nbuf; ++i) {
	new_vbufs.push_back(vbufs[(vidx+i+1) % vbufs.size()]);
      }
      new_vbufs.push_back(bcf_init());
      new_vbufs.swap(vbufs);
      vidx = vbufs.size()-1;
    }
  }
  else {
    vidx = (vidx + 1) % vbufs.size();
  }
  ++nbuf;

  if ( mode_extract ) {
    if ( variants2extract.empty() ) {
      vidx = (vidx + vbufs.size() - 1) % vbufs.size(); 
      eof = true;
      return NULL;
    }
    std::set<variantKeyS>::iterator it = variants2extract.begin();
    if ( cdr.jump_to(it->chrom.c_str(), it->pos) ) {
      while( cdr.read(vbufs[vidx]) ) {
	++nRead;

	if ( nRead % verbose == 0 )
	  notice("Reading %d variants, Skipping %d, Missing %d", nRead, nSkip, nMiss);
	
	if ( ( it->chrom.compare(bcf_get_chrom(cdr.hdr,vbufs[vidx])) != 0 ) || ( vbufs[vidx]->pos > it->pos ) ) {
	  // the marker was not found;
	  ++nMiss;
	  variants2extract.erase(it);
	  vidx = (vidx + vbufs.size() - 1) % vbufs.size();
	  --nbuf;
	  return read();
	}

	std::set<variantKeyS>::iterator it2 = variants2extract.find(*it);
	if ( it2 != variants2extract.end() ) { // found something
	  variants2extract.erase(it2);
	  if ( passed_vfilter() ) {
	    // set ploidies here
	    bcf1_t* v = cursor();
	    set_ploidies_by_sex(v);
	    return v;
	  }
	  else {
	    ++nSkip;
	    vidx = (vidx + vbufs.size() - 1) % vbufs.size();
	    --nbuf;
	    return read();
	  }
	}
      }
      ++nMiss;
      variants2extract.erase(it);
      vidx = (vidx + vbufs.size() - 1) % vbufs.size();
      --nbuf;      
      return read();
    }
    else {
      ++nMiss;
      variants2extract.erase(it);
      vidx = (vidx + vbufs.size() - 1) % vbufs.size();
      --nbuf;      
      return read();
    }
  }
  else {
    while( cdr.read(vbufs[vidx]) ) {
      ++nRead;
      if ( nRead % verbose == 0 )
	notice("Reading %d variants at %s:%d, Skipping %d, Missing %d.", nRead, bcf_hdr_id2name(cdr.hdr, vbufs[vidx]->rid), vbufs[vidx]->pos+1, nSkip, nMiss);      
      if ( passed_vfilter() ) return cursor();
      else {
	//vidx = (vidx + vbufs.size() - 1) % vbufs.size();
	//--nbuf;	
	++nSkip;
      }
    }
    vidx = (vidx + vbufs.size() - 1) % vbufs.size();
    --nbuf;    
    eof = true;
    return NULL;
  }
}

