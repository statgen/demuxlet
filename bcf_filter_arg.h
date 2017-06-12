#ifndef __BCF_FILTER_ARG_H
#define __BCF_FILTER_ARG_H

//extern "C" {
#include "filter.h"
//}
#include "Error.h"
#include "hts_utils.h"

// others
#define MASK_GT_MISS   0x01
#define MASK_GT_HOMREF 0x02
#define MASK_GT_HET    0x04
#define MASK_GT_HOMALT 0x08
#define MASK_GT_NONREF (MASK_GT_HET|MASK_GT_HOMALT)
#define MASK_GT_NOMISS (MASK_GT_HOMREF|MASK_GT_HET|MASK_GT_HOMALT)
#define MASK_GT_ALL    (MASK_GT_MISS|MASK_GT_HOMREF|MASK_GT_HET|MASK_GT_HOMALT)

#define FLT_INCLUDE 1
#define FLT_EXCLUDE 2


struct _bcf_vfilter_arg {
  std::vector<std::string> required_filters; // require at least one of the
  std::string include_expr;
  std::string exclude_expr;
  std::vector<int32_t> req_flt_ids;
  filter_t* filt;
  int32_t filter_logic;

  // allele frequency filters
  int32_t minAC;
  int32_t minMAC;
  int32_t maxAC;
  int32_t maxMAC;
  double minAF;
  double minMAF;
  double maxAF;
  double maxMAF;
  int32_t minAN;
  double minCallRate;
  double minQual;
  int32_t maxAlleles;
  int32_t minAlleles;
  bool snpOnly;

  bool require_GT;
  bool require_PL;

  double probThin;

  _bcf_vfilter_arg() {
    minAC = 0;
    minMAC = 0;
    maxAC = INT_MAX;
    maxMAC = INT_MAX;
    minAF = 0;
    minMAF = 0;
    maxAF = 1.;
    maxMAF = 1.;
    
    minAN = 0;
    minCallRate = 0;

    minQual = 0;
    maxAlleles = INT_MAX;
    minAlleles = 0;
    snpOnly = false;

    filt = NULL;
    filter_logic = 0;
    
    require_GT = false;
    require_PL = false;

    probThin = 1.0;
  }

  void init(bcf_hdr_t* hdr) {
    // initialize variant filter
    std::string filter_str;
    if ( include_expr.empty() ) {
      if ( exclude_expr.empty() ) {
	// do nothing
      }
      else {
	filter_str = exclude_expr;
	filter_logic |= FLT_EXCLUDE;
      }
    }
    else {
      if ( exclude_expr.empty() ) {
	filter_str = include_expr;
	filter_logic |= FLT_INCLUDE;      
      }
      else {
	error("[E:%s:%d %s] Cannot use both --include-expr and --exclude-expr options",__FILE__,__LINE__,__FUNCTION__);
      }    
    }

    if ( filter_logic != 0 )
      filter_init(hdr, filter_str.c_str());

    if ( !required_filters.empty() ) {
      for(int32_t i=0; i < (int32_t)required_filters.size(); ++i) {
	req_flt_ids.push_back(bcf_hdr_id2int(hdr, BCF_DT_ID, required_filters[i].c_str()));
      }
    }

    if ( ( minAC > 0 ) || ( minMAC > 0 ) || ( minAF > 0 ) || ( minMAF > 0 ) ||
	 ( maxAC < INT_MAX ) || ( maxMAC < INT_MAX ) || ( maxAF < 1 ) || ( maxMAF < 0 ) ||
	 ( minAN > 0 ) || ( minCallRate > 0 ) )
      require_GT = true;
  }
};

struct _bcf_gfilter_arg {
  int32_t minDP;
  int32_t minGQ;
  //int32_t minAD;

  _bcf_gfilter_arg() {
    minDP = 0;
    minGQ = 0;
  }
};

typedef struct _bcf_vfilter_arg bcf_vfilter_arg;
typedef struct _bcf_gfilter_arg bcf_gfilter_arg;


#endif
