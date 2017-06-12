#ifndef __BCF_FILTERED_READER_H
#define __BCF_FILTERED_READER_H

#include <cstdlib>
#include <cstring>
#include <vector>

#include "htslib/vcf.h"
#include "htslib/hts.h"
#include "htslib/tbx.h"
#include "filter.h"
#include "hts_utils.h"
#include "cramore.h"
#include "bcf_chunked_reader.h"
#include "genomeChunk.h"
#include "tsv_reader.h"
#include "bcf_filter_arg.h"
#include "bcf_variant_key.h"
#include "Error.h"

class BCFFilteredReader {
public:
  // arguments to be exposed to parameters
  std::string bcf_file_name;
  std::string ref_file_name;
  std::string interval_file_name;
  std::string target_region;
  std::string target_interval_list;
  std::string sample_id_list;

  // argument related to gender map;
  int32_t xStart;
  int32_t xStop;
  std::string xLabel;
  std::string yLabel;
  std::string mtLabel;
  std::string sexMap;
  std::map<std::string,int8_t> mSex;
  bool isX;
  int32_t xRid;
  int32_t yRid;
  int32_t mtRid;
  int8_t sex_ploidies[2];
  
  int32_t unit;

  //bool jumping_extract; // not used yet
  int32_t max_jumping_distance; // not used yet
  
  bcf_vfilter_arg vfilt;
  bcf_gfilter_arg gfilt;
  
  int32_t verbose;

  // internal arguments irrelevant to initial parameters
  bool mode_extract;
  std::set<variantKeyS> variants2extract;
  genomeLoci target_loci;
  
  int32_t nRead;
  int32_t nSkip;
  int32_t nMiss;

  int32_t vidx;
  std::vector<bcf1_t*> vbufs;
  bool unlimited_buffer;
  int32_t nbuf;
  bool eof;

  BCFChunkedReader cdr;
  
  int32_t* gts;
  int32_t n_gts;
  int32_t* pls;
  int32_t n_pls;
  float* dss;
  int32_t n_dss;
  float* gps;
  int32_t n_gps;
  void* flds;
  int32_t n_flds;
  int8_t* ploidies;
  
  std::vector<double> acs;
  int32_t an;

  std::string varID;

  std::set<std::string> sm_ids;  // sample ids to focus on
  std::vector<std::string> v_sm_ids; //
  std::vector<int32_t> sm_icols;  // columns of samples to extract
  std::vector<int8_t> sm_isexes; // sex information for samples to extract

 BCFFilteredReader() : xStart(2699520), xStop(154931044), xLabel("X"), yLabel("Y"), mtLabel("MT"), isX(false), xRid(-1), yRid(-1), mtRid(-1), unit(INT_MAX), max_jumping_distance(0), verbose(10000), mode_extract(false), nRead(0), nSkip(0), nMiss(0), vidx(-1), unlimited_buffer(false), nbuf(0), eof(false), gts(NULL), n_gts(0), pls(NULL), n_pls(0), dss(NULL), n_dss(0), gps(NULL), n_gps(0), flds(NULL), n_flds(0), an(0) {}

  ~BCFFilteredReader() {
    for(int32_t i=0; i < (int32_t)vbufs.size(); ++i) {
      bcf_destroy(vbufs[i]);
    }
    if ( pls ) free(pls);
    if ( gts ) free(gts);
    if ( dss ) free(dss);
    if ( flds ) free(flds);
    if ( ploidies ) delete [] ploidies;
  }

  // initialization function
  void init_params();

  // buffer management 
  void set_buffer_size(int32_t buffer_size);
  inline bcf1_t* cursor() { return vbufs[vidx]; }
  inline bcf1_t* cursor(int32_t idx) {
    if ( idx < nbuf ) {
      return vbufs[(vidx + vbufs.size() - idx) % vbufs.size()];
    }
    else {
      error("[E:%s:%d %s] Cannot move cursor to %d-th variant, which is out of bound from %d", __FILE__, __LINE__, __PRETTY_FUNCTION__, idx, nbuf);
      return 0;      
    }
  }  
  bcf1_t* read();
  void pop() { --nbuf; }
  int32_t clear_buffer_before(const char* chr = NULL, int32_t pos1 = INT_MAX);

  // filtering
  bool passed_vfilter(bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL);
  bool jump_to(const char* chr = NULL, int32_t pos = INT_MAX);

  // variant extraction
  bool add_variant_to_extract(bcf_hdr_t* hdr, bcf1_t* v);

  // parse contents
  bool parse_genotypes(bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL);
  bool parse_likelihoods(bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL, const char* name = "PL");
  bool parse_posteriors(bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL, const char* name = "GP", double gt_error = 1e-4);  
  bool parse_dosages(bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL, const char* name = "DS");
  bool parse_int_fields(const char* name, bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL);
  bool parse_float_fields(const char* name, bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL);
  std::string& get_var_ID(bcf_hdr_t* hdr=NULL, bcf1_t* v=NULL);
  bool set_ploidies_by_sex(bcf1_t* v=NULL);

  // inline functions
  inline int32_t get_genotype_at(int32_t i) {
    int32_t a1 = bcf_gt_allele(gts[sm_icols[i]*2]);
    int32_t a2 = bcf_gt_allele(gts[sm_icols[i]*2+1]);
    if ( ( a1 < 0 ) || ( a2 < 0 ) ) return -1;
    else return bcf_alleles2gt(a1,a2);
  }

  inline int32_t get_allele_at(int32_t i) {
    return bcf_gt_allele(gts[sm_icols[i/2]*2+(i%2)]);
  }

  inline int32_t get_likelihood_at(int32_t i, int32_t ngenotypes = 3) {
    return pls[sm_icols[i/ngenotypes]*ngenotypes + i % ngenotypes];
  }

  inline float get_posterior_at(int32_t i, int32_t ngenotypes = 3) {
    return gps[sm_icols[i/ngenotypes]*ngenotypes + (i % ngenotypes)];
  }
  
  inline float get_dosage_at(int32_t i, int32_t nalleles = 2) {
    return dss[sm_icols[i/(nalleles-1)]+i%(nalleles-1)];
  }

  inline const char* get_sample_id_at(int32_t i) {
    return bcf_hdr_sample_id(cdr.hdr, sm_icols[i]);
  }

  inline int32_t get_nsamples() { return (int32_t)sm_icols.size(); }

  inline const std::string& get_sample_id(int32_t idx) { return v_sm_ids[sm_icols[idx]]; }

  inline bool add_specified_sample(const char* id) { return sm_ids.insert(id).second; }

  inline double get_af(int32_t allele, double pseudocount = 1.0) {
    return (acs[allele]+pseudocount/acs.size())/(an+pseudocount);
  }
};

#endif
