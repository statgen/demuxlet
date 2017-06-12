#ifndef __BCF_VARIANT_KEY_H
#define __BCF_VARIANT_KEY_H

#include <string>
#include <vector>

#include "hts_utils.h"

class variantKeyS {
 public:
  std::string chrom;
  int32_t pos;
  int32_t rlen;
  std::vector<std::string> alleles;

  variantKeyS(bcf_hdr_t* hdr, bcf1_t* v) : chrom(bcf_get_chrom(hdr,v)), pos(v->pos), rlen(v->rlen) {
    alleles.resize(v->n_allele);
    bcf_unpack(v,BCF_UN_STR);    
    for(int32_t i=0; i < v->n_allele; ++i) {
      alleles[i].assign(v->d.allele[i]);
    }

    //notice("Parsed variant key is %s : %d : %s : %u", chrom.c_str(), pos, rlen, alleles[0].c_str(), alleles.size());    
  }

  variantKeyS(const char* s) {
    const char* p = s;
    const char* pp = NULL;    
    int32_t ncolons = 0;
    alleles.resize(2);
    while( *p != '\0' ) {
      if ( *p == ':' ) {
	if ( ncolons == 0 ) {
	  chrom.assign(s, p-s);
	  pp = p;
	  ++ncolons;
	}
	else if ( ncolons == 1 ) {
	  pos = atoi(pp+1)-1;
	  pp = p;
	  ++ncolons;
	}
	else if ( ncolons == 2 ) {
	  rlen = p-pp-1;
	  alleles[0].assign(pp+1,rlen);
	  pp = p;				
	  ++ncolons;
	}
	else if ( ncolons == 3 ) {
	  break;
	}
      }
      ++p;
    }

    if ( ncolons == 3 ) {
      const char* p1 = pp+1;
      const char* pp1 = pp+1;
      //alleles.resize(2);
      while( *p1 != *p ) {
	if ( *p1 == ',' ) {
	  alleles.back().assign(pp1,p1-pp1);
	  pp1 = p1+1;
	  alleles.resize(alleles.size()+1);
	}
	++p1;
      }
      alleles.back().assign(pp1, p-pp1);
    }
    else
      error("[E:%s] Cannot parse variant key %s",__PRETTY_FUNCTION__, s);

    //notice("Parsed variant key for %s is %s : %d : %d : %s : %u", s, chrom.c_str(), pos, rlen, alleles[0].c_str(), alleles.size());
  }

  bool operator<(const variantKeyS& rhs) const {
    if ( chrom == rhs.chrom ) {
      if ( pos == rhs.pos ) {
	if ( rlen == rhs.rlen ) {
	  for(int32_t i=0; i < (int32_t)alleles.size(); ++i) {
	    if ( i >= (int)rhs.alleles.size() ) return false;
	    else if ( alleles[i] != rhs.alleles[i] ) return alleles[i] < rhs.alleles[i];
	  }
	  return false;
	}
	else return ( rlen < rhs.rlen );
      }
      else return ( pos < rhs.pos );
    }
    else {
      int32_t ichr1 = atoi(chrom.c_str());
      int32_t ichr2 = atoi(rhs.chrom.c_str());
      if ( ichr1 == ichr2 ) { return chrom < rhs.chrom; }
      else { return ichr1 < ichr2; }
    }
  }

  bool operator==(const variantKeyS& rhs) const {
    if ( chrom == rhs.chrom ) {
      if ( pos == rhs.pos ) {
	if ( rlen == rhs.rlen ) {
	  if ( alleles.size() != rhs.alleles.size() ) return false;
	  for(int32_t i=0; i < (int32_t)alleles.size(); ++i) {
	    if ( alleles[i] != rhs.alleles[i] ) return false;
	  }
	  return true;
	}
	else return false;
      }
      else return false;
    }
    else return false;
  }  
};

#endif
