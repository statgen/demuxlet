#ifndef __SC_DROP_SEQ_H
#define __SC_DROP_SEQ_H

#include <map>
#include <string>
#include <vector>
#include <stdint.h>

#include "Error.h"

// read : read ID
// unique read : (barcode, snp, umi)
// snp  : snp ID

// read -> cell -> UMI


// 8bit  - BQ
// 8bit  - allele
// 16bit - count
typedef std::map<std::string,uint32_t> sc_snp_droplet_t;
typedef std::map<std::string,uint32_t>::iterator sc_snp_droplet_it_t;

class sc_snp_t {
 public:
  int32_t rid;
  int32_t pos;
  char ref;
  char alt;
  double af;
  double* gps;
};

class sc_dropseq_lib_t {
 public:
  // vector containing SNP & genotype info, index is snp_id
  std::vector<sc_snp_t> snps;

  // mapper between barcode -> bcd_id  
  std::map<std::string,int32_t> bc_map;
  
  // cell_umis[i]->[j] contains the map of UMIs overlapping with snp j in cell i
  std::vector< std::map<int32_t,sc_snp_droplet_t*> > cell_umis;

  // Number of pass-filtered reads and unique reads
  std::vector<int32_t> cell_totl_reads;  
  std::vector<int32_t> cell_pass_reads;
  std::vector<int32_t> cell_uniq_reads;  
  
  std::vector< std::map<int32_t,sc_snp_droplet_t*> > snp_umis;
  int32_t nbcs;
  int32_t nsnps;
  int32_t add_snp(int32_t _rid, int32_t _pos, char _ref, char _alt, double _af, double* _gps);
  int32_t add_cell(const char* barcode);
  bool add_read(int32_t snpid, int32_t cellid, const char* umi, char allele, char qual);

 sc_dropseq_lib_t() : nbcs(0), nsnps(0) {}
};

#endif
