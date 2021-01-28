#include "cramore.h"
#include "bcf_filtered_reader.h"
#include "sam_filtered_reader.h"
#include "sc_drop_seq.h"

int32_t main(int32_t argc, char** argv) {
  SAMFilteredReader sr;
  BCFFilteredReader vr;
  std::string field("GP");
  double genoError = 0.01;
  std::string outPrefix;
  std::string tagGroup("CB");
  std::string tagUMI("UB");
  int32_t capBQ = 40;
  int32_t minBQ = 13;
  int32_t minTD = 0;
  sr.filt.exclude_flag = 0x0f04;
  sr.filt.minMQ = 20;  
  std::vector<std::string> smIDs;
  std::vector<double> gridAlpha;
  //std::vector<double> gridASE;
  vr.verbose = 10000;
  sr.verbose = 1000000;  
  vr.vfilt.minMAC = 1;
  vr.vfilt.minCallRate = 0.5;
  vr.vfilt.maxAlleles = 2;  
  bool writePair = false;
  //bool fullPair = true;
  double doublet_prior = 0.5;
  std::string groupList;
  int32_t minTotalReads = 0;
  int32_t minUniqReads = 0;
  int32_t minCoveredSNPs = 0;

  paramList pl;

  BEGIN_LONG_PARAMS(longParameters)
    LONG_PARAM_GROUP("Options for input SAM/BAM/CRAM", NULL)
    LONG_STRING_PARAM("sam",&sr.sam_file_name, "Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed")
    LONG_STRING_PARAM("tag-group",&tagGroup, "Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB")
    LONG_STRING_PARAM("tag-UMI",&tagUMI, "Tag representing UMIs. For 10x genomiucs, use UB")

    LONG_PARAM_GROUP("Options for input VCF/BCF", NULL)
    LONG_STRING_PARAM("vcf",&vr.bcf_file_name, "Input VCF/BCF file, containing the individual genotypes (GT), posterior probability (GP), or genotype likelihood (PL)")
    LONG_STRING_PARAM("field",&field,"FORMAT field to extract the genotype, likelihood, or posterior from")
    LONG_DOUBLE_PARAM("geno-error",&genoError,"Genotype error rate (must be used with --field GT)")
    LONG_INT_PARAM("min-mac",&vr.vfilt.minMAC, "Minimum minor allele frequency")
    LONG_DOUBLE_PARAM("min-callrate",&vr.vfilt.minCallRate, "Minimum call rate")    
    LONG_MULTI_STRING_PARAM("sm",&smIDs, "List of sample IDs to compare to (default: use all)")
    LONG_STRING_PARAM("sm-list",&vr.sample_id_list, "File containing the list of sample IDs to compare")        

    LONG_PARAM_GROUP("Output Options", NULL)
    LONG_STRING_PARAM("out",&outPrefix,"Output file prefix")
    LONG_MULTI_DOUBLE_PARAM("alpha",&gridAlpha, "Grid of alpha to search for (default is 0, 0.5)")
    LONG_PARAM("write-pair",&writePair, "Writing the (HUGE) pair file")
    LONG_DOUBLE_PARAM("doublet-prior",&doublet_prior, "Prior of doublet")
    LONG_INT_PARAM("sam-verbose",&sr.verbose, "Verbose message frequency for SAM/BAM/CRAM")
    LONG_INT_PARAM("vcf-verbose",&vr.verbose, "Verbose message frequency for VCF/BCF")
    
    LONG_PARAM_GROUP("Read filtering Options", NULL)
    LONG_INT_PARAM("cap-BQ", &capBQ, "Maximum base quality (higher BQ will be capped)")
    LONG_INT_PARAM("min-BQ", &minBQ, "Minimum base quality to consider (lower BQ will be skipped)")
    LONG_INT_PARAM("min-MQ", &sr.filt.minMQ, "Minimum mapping quality to consider (lower MQ will be ignored)")
    LONG_INT_PARAM("min-TD", &minTD, "Minimum distance to the tail (lower will be ignored)")
    LONG_INT_PARAM("excl-flag", &sr.filt.exclude_flag, "SAM/BAM FLAGs to be excluded")

    LONG_PARAM_GROUP("Cell/droplet filtering options", NULL)
    LONG_STRING_PARAM("group-list",&groupList, "List of tag readgroup/cell barcode to consider in this run. All other barcodes will be ignored. This is useful for parallelized run")    
    LONG_INT_PARAM("min-total", &minTotalReads, "Minimum number of total reads for a droplet/cell to be considered")
    LONG_INT_PARAM("min-uniq", &minUniqReads, "Minimum number of unique reads (determined by UMI/SNP pair) for a droplet/cell to be considered")
    LONG_INT_PARAM("min-snp", &minCoveredSNPs, "Minimum number of SNPs with coverage for a droplet/cell to be considered")
  END_LONG_PARAMS();

  pl.Add(new longParams("Available Options", longParameters));
  pl.Read(argc, argv);
  pl.Status();

  if ( gridAlpha.empty() ) {
    gridAlpha.push_back(0);    
    //gridAlpha.push_back(0.05);
    //gridAlpha.push_back(0.1);
    //gridAlpha.push_back(0.15);
    //gridAlpha.push_back(0.2);
    //gridAlpha.push_back(0.25);        
    //gridAlpha.push_back(0.3);
    //gridAlpha.push_back(0.35);    
    //gridAlpha.push_back(0.4);
    //gridAlpha.push_back(0.45);    
    gridAlpha.push_back(0.5);    
  }

  std::set<std::string> bcdSet;
  if ( !groupList.empty() ) {
    tsv_reader tsv_group_list(groupList.c_str());
    while( tsv_group_list.read_line() > 0 ) {
      bcdSet.insert(tsv_group_list.str_field_at(0));
    }
    notice("Finished loading %u droplet/cell barcodes to consider", bcdSet.size());
  }

  for(int32_t i=0; i < (int32_t)smIDs.size(); ++i) {
    vr.add_specified_sample(smIDs[i].c_str());
  }

  vr.unlimited_buffer = true;
  vr.vfilt.maxAlleles = 2;
  sr.set_buffer_size(1);
  //sr.unlimited_buffer = true;
  
  vr.init_params();
  sr.init_params();

  int32_t n_warning_no_gtag = 0;
  int32_t n_warning_no_utag = 0;  
  
  if ( outPrefix.empty() )
    error("[E:%s:%d %s] --out parameter is missing",__FILE__,__LINE__,__PRETTY_FUNCTION__);

  char gtag[2] = {0,0};
  char utag[2] = {0,0};    

  if ( tagGroup.empty() ) { // do nothing
  }
  else if ( tagGroup.size() == 2 ) {
    gtag[0] = tagGroup.at(0);
    gtag[1] = tagGroup.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize group tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagGroup.c_str());
  }

  if ( tagUMI.empty() ) { // do nothing
  }
  else if ( tagUMI.size() == 2 ) {
    utag[0] = tagUMI.at(0);
    utag[1] = tagUMI.at(1);    
  }
  else {
    error("[E:%s:%d %s] Cannot recognize UMI tag %s. It is suppose to be a length 2 string",__FILE__,__LINE__,__FUNCTION__,tagUMI.c_str());
  }    

  // scan VCF and CRAM simultaneously
  // read a variant first

  sc_dropseq_lib_t scl;

  std::vector<int32_t> snpids;
  //std::vector<int32_t> cellids;  
  
  if ( !vr.read() )
    error("[E:%s Cannot read any single variant from %s]", __PRETTY_FUNCTION__, vr.bcf_file_name.c_str());
  
  if ( !vr.parse_posteriors(vr.cdr.hdr, vr.cursor(), field.c_str(), genoError) )
    error("[E:%s] Cannot parse posterior probability at %s:%d", __PRETTY_FUNCTION__, bcf_hdr_id2name(vr.cdr.hdr,vr.cursor()->rid), vr.cursor()->pos+1);


  // check if the chromosome names are in the same order between BCF and SAM
  std::map<int32_t,int32_t> rid2tids;
  std::map<int32_t,int32_t> tid2rids;
  int32_t ntids = bam_hdr_get_n_targets(sr.hdr);
  int32_t prevrid = -1;
  for(int32_t i=0; i < ntids; ++i) {
    const char* chrom = bam_get_chromi(sr.hdr, i);
    int32_t rid = bcf_hdr_name2id(vr.cdr.hdr, chrom);
    if ( rid >= 0 ) {
      if ( prevrid >= rid ) {
	const char* prevchrom = bcf_hdr_id2name(vr.cdr.hdr, prevrid);
	error("[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files have different ordering of chromosomes. SAM/BAM/CRAM file has %s before %s, but VCF/BCF file has %s after %s", __PRETTY_FUNCTION__, prevchrom, chrom, prevchrom, chrom);
      }
      rid2tids[rid] = i;
      tid2rids[i] = rid;
      prevrid = rid;
    }
  }

  if ( rid2tids.empty() || tid2rids.empty() || ( rid2tids.size() != tid2rids.size() ) ) {
    error("[E:%s] Your VCF/BCF files and SAM/BAM/CRAM files does not have any matching chromosomes, or some chromosome names are duplicated");
  }
  
  int32_t nv = vr.get_nsamples();
  double* gps = new double[nv*3];
  for(int32_t i=0; i < nv * 3; ++i) 
    gps[i] = vr.get_posterior_at(i);
  int32_t snpid = scl.add_snp( vr.cursor()->rid, vr.cursor()->pos, vr.cursor()->d.allele[0][0], vr.cursor()->d.allele[1][0], vr.get_af(1), gps);
  snpids.push_back(snpid);

  int32_t ibeg = 0;
  char base, qual;
  int32_t rpos;
  kstring_t readseq = {0,0,0};
  kstring_t readqual = {0,0,0};

  int32_t nReadsMultiSNPs = 0, nReadsSkipBCD = 0, nReadsPass = 0, nReadsRedundant = 0, nReadsN = 0, nReadsLQ = 0, nReadsTMP = 0, nNonBiallelicSNPs = 0;
    
  while( sr.read() ) { // read SAM file
    int32_t endpos = bam_endpos(sr.cursor());
    int32_t tid2rid = bcf_hdr_name2id(vr.cdr.hdr, bam_get_chrom(sr.hdr, sr.cursor()));
    if ( tid2rid < 0 ) { // no matching BCF entry in the chromosome, skip;
      continue;
    }

    int32_t n_cleared = vr.clear_buffer_before(  bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos );
    //for(int32_t i=ibeg; i < ibeg+n_cleared; ++i) {
    //  v_umis.clear();
    //}
    ibeg += n_cleared;

    // add new snps
    while( ( !vr.eof ) && ( ( vr.cursor()->rid < tid2rid ) || ( ( vr.cursor()->rid == tid2rid ) && ( vr.cursor()->pos < endpos ) ) ) ) {
      if ( vr.read() ) {
	if ( !vr.parse_posteriors(vr.cdr.hdr, vr.cursor(), field.c_str(), genoError) )
	  error("[E:%s] Cannot parse posterior probability at %s:%d", __PRETTY_FUNCTION__, bcf_hdr_id2name(vr.cdr.hdr,vr.cursor()->rid), vr.cursor()->pos+1);

	// check whether the variant is SNP and biallelic
	if ( ( vr.cursor()->rlen > 1 ) || ( vr.cursor()->n_allele != 2 ) || ( strlen(vr.cursor()->d.allele[0]) > 1) ) {
	  if ( nNonBiallelicSNPs < 10 ) {
	    warning("VCF record must be biallelic SNPs. Ignoring non-SNPs and/or multi-allelic variants at %d:%d", bcf_hdr_id2name(vr.cdr.hdr,vr.cursor()->rid), vr.cursor()->pos+1);
	  }

	  ++nNonBiallelicSNPs;
	  
	  if ( nNonBiallelicSNPs == 10 ) {
	    warning("Suppressing 10+ warnings of the same kind (non-SNP or multi-alleic variants)");
	  }
	}
	
	gps = new double[nv*3];
	for(int32_t i=0; i < nv * 3; ++i) {
	  gps[i] = vr.get_posterior_at(i);
	}
	snpid = scl.add_snp( vr.cursor()->rid, vr.cursor()->pos, vr.cursor()->d.allele[0][0], vr.cursor()->d.allele[1][0], vr.get_af(1), gps);
	snpids.push_back(snpid);
      }
      else {
	//error("Cannot read new SNP");
      }
    }

    // get barcode
    int32_t ibcd = 0;
    if ( tagGroup.empty() ) {
      ibcd = scl.add_cell(".");
    }
    else {
      uint8_t *bcd = (*gtag) ? (uint8_t*) bam_aux_get(sr.cursor(), gtag) : NULL;
      const char* sbcd = ".";
      if ( ( bcd != NULL ) && ( *bcd == 'Z' ) ) {
	sbcd = bam_aux2Z(bcd);
      }
      else {
	if ( n_warning_no_gtag < 10 ) {
	  notice("WARNING: Cannot find Droplet/Cell tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single group", tagGroup.c_str(), sr.n_read, bam_get_qname(sr.cursor()), bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos, bam_endpos(sr.cursor()));
	}
	else if ( n_warning_no_gtag == 10 ) {
	  notice("WARNING: Suppressing 10+ missing Droplet/Cell tag warnings...");
	}
	++n_warning_no_gtag;	
      }
      
      if ( bcdSet.empty() || ( bcdSet.find(sbcd) != bcdSet.end() ) ) {
	ibcd = scl.add_cell(sbcd);
	if ( ( ibcd + 1 == scl.nbcs ) && ( scl.nbcs % 1000 == 0 ) )
	  notice("Observed %d droplets with unique cell barcode", scl.nbcs);
      }
      else {
	++nReadsSkipBCD;
	continue;
      }
    }

    ++nReadsTMP;

    // get UMI
    std::string sumi(".");
    if ( tagUMI.empty() ) {
      catprintf(sumi,"%x",rand()); // give a random UMI
    }
    else {
      uint8_t *umi = (*utag) ? (uint8_t*) bam_aux_get(sr.cursor(), utag) : NULL;      
      if ( ( umi != NULL ) && ( *umi == 'Z' ) ) {
	sumi = bam_aux2Z(umi);
      }
      else {
	if ( n_warning_no_utag < 10 ) {
	  notice("WARNING: Cannot find UMI tag %s from %d-th read %s at %s:%d-%d. Treating all of them as a single UMI", tagUMI.c_str(), sr.n_read, bam_get_qname(sr.cursor()), bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos, bam_endpos(sr.cursor()));
	}
	else if ( n_warning_no_utag == 10 ) {
	  notice("WARNING: Suppressing 10+ UMI warnings...");
	}
	++n_warning_no_utag;
	//error("[E:%s] Cannot find UMI tag %d %d %x %s %s %x", __PRETTY_FUNCTION__, sr.nbuf, sr.ridx, sr.cursor(), bcd, utag, umi);
      }
    }

    ++scl.cell_totl_reads[ibcd];    
    
    // genotype all reads together
    int32_t nv_pass = 0;
    int32_t nv_redundant = 0;
    int32_t nv_valid = 0;
    int32_t allele, bq;

    //if ( rand() % 10000 == 0 )
    //notice("Reading between %s:%d-%d at %s:%d to %d i=beg=%d, nbuf=%d, vidx=%d, size=%u prevpos=%d", bam_get_chrom(sr.hdr, sr.cursor()), sr.cursor()->core.pos+1, bam_endpos(sr.cursor()), bcf_hdr_id2name(vr.cdr.hdr, scl.snps[ibeg].rid), scl.snps[ibeg].pos, scl.snps[ibeg+vr.nbuf-1].pos, ibeg, vr.nbuf, vr.vidx, vr.vbufs.size(), scl.snps[ibeg-1].pos);
    
    for(int32_t i=ibeg; i < ibeg+vr.nbuf; ++i) {
      bam1_t* b = sr.cursor();
      bam_get_base_and_qual_and_read_and_qual(b, (uint32_t)scl.snps[i].pos, base, qual, rpos, &readseq, &readqual);
      if ( rpos == BAM_READ_INDEX_NA ) {
	//if ( rand() % 1000 == 0 ) 
	//notice("Cannot find any informative read between %s:%d-%d at %s:%d", bam_get_chrom(sr.hdr, b), b->core.pos+1, bam_endpos(b), bcf_hdr_id2name(vr.cdr.hdr, scl.snps[i].rid), scl.snps[i].pos);
	continue;
      }
      if ( base == 'N' ) continue;

      ++nv_valid;

      if ( qual-33 < minBQ ) { continue; }
      if ( rpos < minTD-1 ) { continue; }
      if ( rpos + minTD > b->core.l_qseq ) { continue; }

      allele = ( base == scl.snps[i].ref ) ? 0 : ( ( base == scl.snps[i].alt ) ? 1 : 2 );
      bq = qual-33 > capBQ ? capBQ : qual-33;

      if ( scl.add_read(snpids[i], ibcd, sumi.c_str(), allele, bq) )
	++nv_pass;
      else
	++nv_redundant;
    }

    if ( nv_pass > 1 ) ++nReadsMultiSNPs;
    if ( nv_pass > 0 ) ++nReadsPass;
    else if ( nv_redundant > 0 ) ++nReadsRedundant;
    else if ( nv_valid > 0 ) ++nReadsLQ;
    else ++nReadsN;

    //if ( nv_valid > 0 ) ++scl.cell_totl_reads[ibcd];    
  }

  if ( n_warning_no_utag > 10 ) 
    notice("WARNING: Suppressed a total of %d UMI warnings...", n_warning_no_utag);
    
  if ( n_warning_no_gtag > 10 ) 
    notice("WARNING: Suppressed a total of %d droplet/cell barcode warnings...", n_warning_no_gtag);
  
  notice("Finished reading %d markers from the VCF file", (int32_t)snpids.size());

  int32_t nAlpha = (int32_t)gridAlpha.size();

  //notice("Finished processing %d reads across %d variants across %d barcodes", nReadsPass, (int32_t)v_poss.size(), (int32_t)bcMap.size(), (int32_t)bcMap.size());
  notice("Total number input reads : %d", sr.n_read);
  notice("Total number valid droplets observed : %d", scl.nbcs);
  notice("Total number valid SNPs observed     : %d", scl.nsnps);    
  notice("Total number of read-QC-passed reads : %d ", sr.n_read - sr.n_skip); //, nReadsN + nReadsUnique + nReadsLQ + nReadsPass);
  notice("Total number of skipped reads with ignored barcodes : %d", nReadsSkipBCD);
  notice("Total number of non-skipped reads with considered barcodes : %d", nReadsTMP);  
  notice("Total number of gapped/noninformative reads : %d", nReadsN);
  notice("Total number of base-QC-failed reads : %d", nReadsLQ);  
  notice("Total number of redundant reads : %d", nReadsRedundant);
  notice("Total number of pass-filtered reads : %d", nReadsPass);
  notice("Total number of pass-filtered reads overlapping with multiple SNPs : %d", nReadsMultiSNPs);

  //notice("Finished processing %d reads across %d variants across %d barcodes, filtering %d (%.2lf%%) reads, including %d (%.2lf%%) gapped reads, %d (%.2lf%%) low quality reads, and %d (%.2lf%%) redundant/qcfail reads from the BAM file %s", nReadsPass, (int32_t)v_poss.size(), (int32_t)bcMap.size(), nReadsLQ + nReadsUnique + nReadsN, 100.0 * (nReadsLQ + nReadsUnique + nReadsN) / (nReads, nReadsN, 100.0 * nReadsN / (nReadsPass + nReadsLQ + nReadsUnique + nReadsN), nReadsLQ, 100.0 * nReadsLQ / nR, nReadsRedundant, 100.0 * nReadsRedundant / nReadsAll,  inSam.c_str());

  sr.close();
  //vr.close();

  notice("Starting to prune out cells with too few reads...");
  int32_t nRemoved = 0;
  if ( minTotalReads + minUniqReads + minCoveredSNPs < 0 ) {
    for(int32_t i=0; i < scl.nbcs; ++i) {
      if ( ( scl.cell_totl_reads[i] < minTotalReads ) || ( scl.cell_uniq_reads[i] < minUniqReads) || ( (int32_t)scl.cell_umis[i].size() < minCoveredSNPs ) ) {
	for(std::map<int32_t,sc_snp_droplet_t*>::iterator it = scl.cell_umis[i].begin();
	    it != scl.cell_umis[i].end(); ++it) {
	  delete it->second;
	  scl.snp_umis[it->first].erase(i);
	}
	scl.cell_umis[i].clear();
	++nRemoved;
      }
    }
  }
  notice("Finishing pruning out %d cells with too few reads...", nRemoved);


  if ( (int32_t)snpids.size() != scl.nsnps )
    error("[E:%s snpids.size() = %u != scl.nsnps = %d",__PRETTY_FUNCTION__, snpids.size(), scl.nsnps);


  // Calculate average genotype probability
  double* gp0s = (double*) calloc(scl.nsnps * 3, sizeof(double)); 
  for(int32_t i=0; i < scl.nsnps; ++i) {
    for(int32_t j=0; j < nv; ++j) {
      gp0s[i*3] += scl.snps[i].gps[3*j];
      gp0s[i*3+1] += scl.snps[i].gps[3*j+1];
      gp0s[i*3+2] += scl.snps[i].gps[3*j+2];      
    }
    gp0s[i*3] /= nv;
    gp0s[i*3+1] /= nv;
    gp0s[i*3+2] /= nv;    
  }
  
  // start evaluating genotype concordances
  // calculate for (nBcd) x (nInds) to find the best matching genotypes first
  notice("Starting to identify best matching individual IDs");

  htsFile* wsingle = hts_open((outPrefix+".single").c_str(),"w");

  if ( wsingle == NULL )
    error("[E:%s:%d %s] Cannot create %s.single file",__FILE__,__LINE__,__FUNCTION__,outPrefix.c_str());
  
  std::vector<double> llks(scl.nbcs * nv, 0);
  std::vector<double> llk0s(scl.nbcs, 0);  
  double tmp;
  for(int32_t i=0; i < scl.nsnps; ++i) {
    if ( ( vr.verbose > 0 ) && ( (i+1) % vr.verbose == 0 ) )
      notice("Processing %d markers...",i+1);

    std::map<int32_t,sc_snp_droplet_t*>& cells = scl.snp_umis[i];
    if ( cells.empty() ) continue;
    std::map<int32_t,sc_snp_droplet_t*>::iterator it;
    //std::vector<double> GLs(scl.nbcs * 4, 0);

    double GLs[3];

    for(it = cells.begin(); it != cells.end(); ++it) {
      GLs[0] = GLs[1] = GLs[2] = 1.0;
      for(sc_snp_droplet_it_t it2=it->second->begin(); it2 != it->second->end(); ++it2) {
	uint8_t al = ( it2->second >> 24 ) & 0x00ff;
	uint8_t bq = ( it2->second >> 16 ) & 0x00ff;

	//if ( rand() % 1000 == 0 ) notice("bq = %d, al = %d", bq, al);
	//uint32_t ibcd = it->first;

	if ( al == 2 ) continue;

	GLs[0] *= ((al==0) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0);
	GLs[1] *= (0.5 - phredConv.phred2Err[bq]/3.0);
	GLs[2] *= ((al==1) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0);
	tmp = GLs[0] + GLs[1] + GLs[2];
	GLs[0] /= tmp;
	GLs[1] /= tmp;
	GLs[2] /= tmp;
      }

      GLs[0] += 1e-6;
      GLs[1] += 1e-6;
      GLs[2] += 1e-6;
      tmp = GLs[0] + GLs[1] + GLs[2];
      GLs[0] /= tmp;
      GLs[1] /= tmp;
      GLs[2] /= tmp;  
      
      gps = scl.snps[i].gps;
      for(int32_t k=0; k < nv; ++k) {
	llks[it->first * nv + k] += log(GLs[0]*gps[k*3] + GLs[1]*gps[k*3+1] + GLs[2]*gps[k*3+2]);
	//if ( rand() % 1000 == 0 ) notice("%lg %lg %lg",gps[k*3],gps[k*3+1],gps[k*3+2]);
      }
      llk0s[it->first] += log( GLs[0] * gp0s[i*3] + GLs[1] * gp0s[i*3+1] + GLs[2] * gp0s[i*3+2] );
    }
  }
 

  // find the best matching individual
  std::vector<int32_t> iBest(scl.nbcs,0);
  std::vector<int32_t> iNext(scl.nbcs,0);  
  std::vector<double> llkBest(scl.nbcs,0);
  notice("Identifying best-matching individual..");

  hprintf(wsingle, "BARCODE\tSM_ID\tRD.TOTL\tRD.PASS\tRD.UNIQ\tN.SNP\tLLK1\tLLK0\tPOSTPRB\n");
  int32_t i=0;
  for(std::map<std::string,int32_t>::iterator it = scl.bc_map.begin();
      it != scl.bc_map.end(); ++it) {
    int32_t imax = -1;
    int32_t inext = -1;
    double maxLLK = -1e300;
    double nextLLK = -1e300;
    double sumLLK = -1e300;

    if ( ( scl.cell_totl_reads[it->second] < minTotalReads ) || ( scl.cell_uniq_reads[it->second] < minUniqReads) || ( (int32_t)scl.cell_umis[it->second].size() < minCoveredSNPs ) ) continue;    

    for(int32_t j=0; j < nv; ++j) {
      double curLLK = llks[it->second * nv + j];
      if ( sumLLK > curLLK ) {
	sumLLK = sumLLK + log(1.0 + exp(curLLK - sumLLK));
      }
      else {
	sumLLK = curLLK + log(1.0 + exp(sumLLK - curLLK));
      }

      if ( curLLK > maxLLK ) {
	inext = imax;
	nextLLK = maxLLK;
	imax = j;
	maxLLK = curLLK;
      }
      else if ( curLLK > nextLLK ) {
	nextLLK = curLLK;
	inext = j;
      }
    }

    
    for(int32_t j=0; j < nv; ++j) {
      double curLLK = llks[it->second * nv + j];      
      hprintf(wsingle,"%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.3lg\n",
	      it->first.c_str(),
	      vr.get_sample_id_at(j),
	      scl.cell_totl_reads[it->second],
	      scl.cell_pass_reads[it->second],
	      scl.cell_uniq_reads[it->second],
	      (int32_t)scl.cell_umis[it->second].size(),
	      curLLK,
	      llk0s[it->second],
	      exp(curLLK-sumLLK)
	      );
    }
    iBest[it->second] = imax;
    iNext[it->second] = inext;
    llkBest[it->second] = maxLLK;

    ++i;
    if ( i % 1000 == 0 )
      notice("Processing %d droplets...", i);
  }
  notice("Finished processing %d droplets total", i);  
  hts_close(wsingle);

  htsFile* wsing2 = hts_open((outPrefix+".sing2").c_str(),"w");
  htsFile* wpair = (writePair ? hts_open((outPrefix+".pair").c_str(),"w") : NULL);
  htsFile* wbest = hts_open((outPrefix+".best").c_str(),"w");

  hprintf(wsing2, "BARCODE\tSM_ID\tRD.TOTL\tRD.PASS\tRD.UNIQ\tN.SNP\tLLK1\tLLK0\tPOSTPRB\n");  

  if ( ( writePair && wpair == NULL ) || ( wsingle == NULL ) )
    error("[E:%s:%d %s] Cannot create %s.single, %s.pair files",__FILE__,__LINE__,__FUNCTION__,outPrefix.c_str(), outPrefix.c_str());


  // start finding the next-best matching individual
  // here we iterate each cell separately.
  // pre-calculate nsnp*nv*nv*9, nv*1*9, 1*nv*9, 1*9
  double* gpAB = new double[scl.nsnps * nv * nv * 9];
  double* gpA0 = new double[scl.nsnps * nv * 9];
  double* gp00 = new double[scl.nsnps * 9];

  int32_t j, k, l, m, n;  
  for(i=0; i < scl.nsnps; ++i) {
    for(j=0; j < nv; ++j) {
      for(k=0; k < nv; ++k) {
	gps = scl.snps[i].gps;
	for(l=0; l < 3; ++l) {
	  for(m=0; m < 3; ++m) {
	    gpAB[i*nv*nv*9 + j*nv*9 + k*9 + l*3 + m] = gps[j*3+l] * gps[k*3+m];
	    gpA0[i*nv*9 + j*9 + l*3 + m] = gps[j*3+l] * gp0s[i*3+m];
	    gp00[i*9 + l*3 + m] = gp0s[i*3+l] * gp0s[i*3+m];	    
	  }
	}
      }
    }
  }

  // iterate each barcode
  int32_t n1 = nv;
  double* llksAB = new double[n1 * nv * nAlpha];
  double* llksA0 = new double[nv * nAlpha];
  double* llks00 = new double[nAlpha];
  //double* postAB = new double[n1 * nv * nAlpha];  

  if ( writePair )
    hprintf(wpair,"BARCODE\tSM1.ID\tSM2.ID\tLLK12\tPOSTPRB\n");
  hprintf(wbest,"BARCODE\tRD.TOTL\tRD.PASS\tRD.UNIQ\tN.SNP\tBEST\tSNG.1ST\tSNG.LLK1\tSNG.2ND\tSNG.LLK2\tSNG.LLK0\tDBL.1ST\tDBL.2ND\tALPHA\tLLK12\tLLK1\tLLK2\tLLK10\tLLK20\tLLK00\tPRB.DBL\tPRB.SNG1\n");
  //SINGLE.BEST.ID\tSINGLE.NEXT.ID\t
  //SM1.ID\tSM2.ID\tALPHA\tRD.TOTL\tRD.PASS\tRD.UNIQ\tN.SNP\tLLK12\tLLK1\tLLK0\tLLK10\tLLK00\tPOSTPRB\n");

  int ncells = 0;
  for(std::map<std::string,int32_t>::iterator it0 = scl.bc_map.begin(); it0 != scl.bc_map.end(); ++it0, ++ncells) {
    if ( ncells % 100 == 0 )
      notice("Processing %d cells..", ncells);
    
    i = it0->second;
    if ( ( scl.cell_totl_reads[i] < minTotalReads ) || ( scl.cell_uniq_reads[i] < minUniqReads) || ( (int32_t)scl.cell_umis[i].size() < minCoveredSNPs ) ) continue;
    
    memset(llksAB,0,sizeof(double)*n1*nv*nAlpha);
    //memset(postAB,0,sizeof(double)*n1*nv*nAlpha);    
    memset(llksA0,0,sizeof(double)*nv*nAlpha);
    memset(llks00,0,sizeof(double)*nAlpha);

    int32_t jbeg = 0; //fullPair ? 0 : iBest[i];
    int32_t jend = nv; //fullPair ? nv : iBest[i]+1;
    
    std::map<int32_t,sc_snp_droplet_t*>& snps = scl.cell_umis[i];
    if ( snps.empty() ) continue;
    std::map<int32_t,sc_snp_droplet_t*>::iterator it;
    std::vector<double> pGs(nAlpha*9,1.0);
    for(it = snps.begin(); it != snps.end(); ++it) {

      std::fill(pGs.begin(),pGs.end(),1.0);
      
      // calculate genotype likelihoods
      for(sc_snp_droplet_it_t it2=it->second->begin(); it2 != it->second->end(); ++it2) {
	uint8_t al = (it2->second >> 24) & 0x00ff;
	uint8_t bq = (it2->second >> 16) & 0x00ff;
	
	if ( al == 2 ) continue;

	double pR = (al == 0) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0;
	double pA = (al == 1) ? phredConv.phred2Mat[bq] : phredConv.phred2Err[bq]/3.0;
	double maxpG = 0;

	for(int32_t k=0; k < nAlpha; ++k) {
	  for(int32_t l=0; l < 3; ++l) {  // 1-Alpha
	    for(int32_t m=0; m < 3; ++m) { // Alpha
	      double p = 0.5*l + (m-l)*0.5*gridAlpha[k]; // %A (0, 0.5a, 1.0a, 0.5-0.5a, 0.5, 0.5+0.5a, 1.0-a, 1.0-0.5a, 1.0)
	      double& pG = pGs[k*9 + l*3 + m];
	      // l m p          pR   pA  
	      // 0 0 0          1-e  e/3  1-e    
	      // 0 1 a/2                  (1-e)(1-a/2) + e/3*a/2             
	      // 0 2 a
	      // 1 0 0.5-a/2
	      // 1 1 0.5        1-e  e/3  0.5-e/3  
	      // 1 2 0.5+a/2
	      // 2 0 1-a
	      // 2 1 1-a/2
	      // 2 2 1          1-e  e/3  e/3
	      pG *= (pR * (1.0-p) + pA * p);	      
	      if ( maxpG < pG )
		maxpG = pG;
	    }
	  }
	}

	for(int32_t k=0; k < nAlpha; ++k) {
	  // normalize
	  for(int32_t l=0; l < 3; ++l) {   // 1-Alpha
	    for(int32_t m=0; m < 3; ++m) { // Alpha
	      pGs[k*9 + l*3 + m] /= maxpG;
	    }
	  }
	}
      }

      // add marginal probability in genotype likelihood
      double maxpG = 0;
      for(int32_t k=0; k < nAlpha; ++k) {
	// normalize
	for(int32_t l=0; l < 3; ++l) {   // 1-Alpha
	  for(int32_t m=0; m < 3; ++m) { // Alpha
	    double& pG = pGs[k*9 + l*3 + m];
	    pG += 1e-6;
	    if ( maxpG < pG )
	      maxpG = pG;	    
	  }
	}
      }
      
      for(int32_t k=0; k < nAlpha; ++k) {
	// normalize
	for(int32_t l=0; l < 3; ++l) {   // 1-Alpha
	  for(int32_t m=0; m < 3; ++m) { // Alpha
	    pGs[k*9 + l*3 + m] /= maxpG;
	  }
	}
      }

      // calculate the sum of posterior probabilities
      int32_t isnp = it->first;
      std::vector<double> sumPs(nAlpha,0);
      double p;
      int32_t j, k, l, m, n;

      for(j=jbeg; j < jend; ++j) {
	// pairwise LLK
	for(k=0; k < nv; ++k) {
	  std::fill(sumPs.begin(), sumPs.end(), 0);
	  for(l=0; l < 3; ++l) {
	    for(m=0; m < 3; ++m) {
	      p = gpAB[isnp*nv*nv*9 + j*nv*9 + k*9 + l*3 + m];
	      for(n=0; n < nAlpha; ++n) 
		sumPs[n] += (p * pGs[n*9+l*3+m]);
	    }
	  }
	  for(n=0; n < nAlpha; ++n)
	    llksAB[(j-jbeg)*nv*nAlpha + k*nAlpha + n] += log(sumPs[n]);
	}

	// A0 LLK
	std::fill(sumPs.begin(), sumPs.end(), 0);
	for(l=0; l < 3; ++l) {
	  for(m=0; m < 3; ++m) {
	    p = gpA0[isnp*nv*9 + j*9 + l*3 + m];
	    for(n=0; n < nAlpha; ++n) 
	      sumPs[n] += (p * pGs[n*9+l*3+m]);
	  }
	}
	for(n=0; n < nAlpha; ++n)
	  llksA0[(j-jbeg)*nAlpha + n] += log(sumPs[n]);	
      }

      // 00 LLK
      std::fill(sumPs.begin(), sumPs.end(), 0);
      for(l=0; l < 3; ++l) {
	for(m=0; m < 3; ++m) {
	  p = gp00[isnp*9 + l*3 + m];
	  for(n=0; n < nAlpha; ++n) 
	    sumPs[n] += (p * pGs[n*9+l*3+m]);
	}
      }
      for(n=0; n < nAlpha; ++n)
	llks00[n] += log(sumPs[n]);
    }

    // normalize by max likelihood
    double maxLLK = -1e300;
    for(j=jbeg; j < jend; ++j) {
      for(k=0; k < nv; ++k) {
	for(n=0; n < nAlpha; ++n) {
	  if ( maxLLK < llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n] )
	    maxLLK = llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n];
	}
      }
    }

    // calculate posterior probability
    double sumSingle = 0, sumDouble = 0;
    for(j=jbeg; j < jend; ++j) {
      sumSingle += (exp(llksAB[(j-jbeg)*nv*nAlpha] - maxLLK)* (1.-doublet_prior) / (jend-jbeg));
      
      for(k=0; k < nv; ++k) {
	if ( j == k ) continue;
	for(n=1; n < nAlpha; ++n) {
	    sumDouble += ( exp(llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n] - maxLLK)* doublet_prior / (jend-jbeg) / (nv-1) / (nAlpha-1) / (gridAlpha[n] == 0.5 ? 2.0 : 1.0));
	}
      }
    }

    /*
    for(j=jbeg; j < jend; ++j) {
      for(k=0; k < nv; ++k) {
	for(n=0; n < nAlpha; ++n) {
	  postAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n] = (exp(llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n] - maxLLK) * (n == 0 ? 1.-doublet_prior : doublet_prior) / (jend-jbeg) / nv / (n == 0 ? 1 : (nAlpha-1)))/(sumSingle+sumDouble);
	}
      }
      }*/


    int32_t iSing1 = -1, iSing2 = -1;
    double maxSing1 = -1e300, maxSing2 = -1e300;
    for(j=jbeg; j < jend; ++j) {
      if ( maxSing1 < llksAB[(j-jbeg)*nv*nAlpha] ) {
	maxSing2 = maxSing1;
	iSing2 = iSing1;
	iSing1 = j;
	maxSing1 = llksAB[(j-jbeg)*nv*nAlpha];
      }
      else if ( maxSing2 < llksAB[(j-jbeg)*nv*nAlpha] ) {
	iSing2 = j;
	maxSing2 = llksAB[(j-jbeg)*nv*nAlpha];	
      }
      
      hprintf(wsing2,"%s\t%s\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.3lg\n",
	      it0->first.c_str(),
	      vr.get_sample_id_at(j),	    	      
	      scl.cell_totl_reads[i],
	      scl.cell_pass_reads[i],
	      scl.cell_uniq_reads[i],
	      (int32_t)scl.cell_umis[i].size(),
	      llksAB[(j-jbeg)*nv*nAlpha],
	      llks00[0],
	      exp(llksAB[(j-jbeg)*nv*nAlpha]-maxLLK) * (1.-doublet_prior) / (jend-jbeg) / sumSingle);
    }    

    if ( writePair )  {
      for(j=jbeg; j < jend; ++j) {
	hprintf(wpair,"%s\t%s\t%s\t%.3lf\t%.5lf\t%.5lg\n",
		it0->first.c_str(),
		vr.get_sample_id_at(j),
		vr.get_sample_id_at(j),
		gridAlpha[0],
		llksAB[(j-jbeg)*nv*nAlpha],
		exp(llksAB[(j-jbeg)*nv*nAlpha]-maxLLK)*(1.-doublet_prior)/(jend-jbeg)/(sumSingle+sumDouble));
	
	for(k=0; k < nv; ++k) {
	  for(n=0; n < nAlpha; ++n) {
	    if ( ( n > 0 ) && ( j != k ) ) {
	      if ( ( j > k ) && ( gridAlpha[n] == 0.5 ) ) continue;
	      hprintf(wpair,"%s\t%s\t%s\t%.3lf\t%.5lf\t%.5lg\n",
		      it0->first.c_str(),
		      vr.get_sample_id_at(j),
		      vr.get_sample_id_at(k),
		      gridAlpha[n],
		      llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n],
		      exp(llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n]-maxLLK)*doublet_prior/(jend-jbeg)/(nv-1)/(nAlpha-1)/(sumSingle+sumDouble));
	    }
	  }
	}
      }
    }

    int jBest = -1, kBest = -1, alphaBest = -1;
    double maxAB = -1e300;

    for(j=jbeg; j < jend; ++j) {
      for(k=0; k < nv; ++k) {
	if ( j == k ) continue;
	for(n=1; n < nAlpha; ++n) {
	  if ( maxAB < llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n] ) {
	    jBest = j;
	    kBest = k;
	    alphaBest = n;
	    maxAB = llksAB[(j-jbeg)*nv*nAlpha+k*nAlpha+n];
	  }
	}
      }
    }

    double singLLK1  = llksAB[(iSing1-jbeg)*nv*nAlpha];
    double singLLK2  = llksAB[(iSing2-jbeg)*nv*nAlpha];    
    //double singLLK1  = llksAB[(iBest[i]-jbeg)*nv*nAlpha];
    //double singLLK2  = llksAB[(iNext[i]-jbeg)*nv*nAlpha];
    double singLLK0  = llks00[0];
    double pairLLK12 = llksAB[(jBest-jbeg)*nv*nAlpha+kBest*nAlpha+alphaBest];
    double pairLLK1  = llksAB[(jBest-jbeg)*nv*nAlpha];
    double pairLLK2  = llksAB[(kBest-jbeg)*nv*nAlpha];
    double pairLLK10 = llksAB[(jBest-jbeg)*nv*nAlpha+alphaBest];
    double pairLLK20 = llksAB[(kBest-jbeg)*nv*nAlpha+alphaBest];        
    double pairLLK00 = llks00[alphaBest];
    double postDoublet = sumDouble/(sumSingle+sumDouble);
    double postSinglet = exp(singLLK1 - maxLLK) * (1.-doublet_prior) / (jend-jbeg) / sumSingle;

    hprintf(wbest,"%s\t%d\t%d\t%d\t%d\t",
	    it0->first.c_str(),
	    scl.cell_totl_reads[i],
	    scl.cell_pass_reads[i],
	    scl.cell_uniq_reads[i],
	    (int32_t)scl.cell_umis[i].size());

    if ( ( pairLLK12 > pairLLK1 ) && ( pairLLK12 > pairLLK2 ) && ( pairLLK12 > singLLK1 + 2 ) )  {
      // best interpretation is doublet
      hprintf(wbest,"DBL-%s-%s-%.3lf",
	    vr.get_sample_id_at(jBest),
	    vr.get_sample_id_at(kBest),
	      gridAlpha[alphaBest]);	      
    }
    else if ( singLLK1 > singLLK2 + 2 ) {
      hprintf(wbest,"SNG-%s",
	      vr.get_sample_id_at(iSing1));
      //vr.get_sample_id_at(iBest[i]));
    }
    else {
      hprintf(wbest,"AMB-%s-%s-%s/%s",
	      vr.get_sample_id_at(iSing1),
	      vr.get_sample_id_at(iSing2),	      	      
	      //vr.get_sample_id_at(iBest[i]),
	      //vr.get_sample_id_at(iNext[i]),	      
	      vr.get_sample_id_at(jBest),
	      vr.get_sample_id_at(kBest));	      
    }

      //hprintf(wbest,"\t%s\t%.4lf",vr.get_sample_id_at(iBest[i]), singLLK1);
    hprintf(wbest,"\t%s\t%.4lf",vr.get_sample_id_at(iSing1), singLLK1);
    //hprintf(wbest,"\t%s\t%.4lf\t%.4lf",vr.get_sample_id_at(iNext[i]), singLLK2, singLLK0);        
    hprintf(wbest,"\t%s\t%.4lf\t%.4lf",vr.get_sample_id_at(iSing2), singLLK2, singLLK0);    
    hprintf(wbest,"\t%s\t%s\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.3lg\t%.3lg\n",
	    vr.get_sample_id_at(jBest),
	    vr.get_sample_id_at(kBest),
	    gridAlpha[alphaBest],
	    pairLLK12,
	    pairLLK1,
	    pairLLK2,
	    pairLLK10,
	    pairLLK20,
	    pairLLK00,
	    postDoublet,
	    postSinglet);
  }
  notice("Finished writing output files");

  if ( writePair ) hts_close(wpair);
  hts_close(wbest);
  hts_close(wsing2);  
  
  return 0;
}
