#ifndef __GENOME_LOCI_H
#define __GENOME_LOCI_H

#include <vector>
#include <set>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <climits>
#include "Error.h"

// A single genomic int32_terval
class genomeLocus {
 public:
  std::string chrom; // chromosome name
  int32_t beg1; // includes 1-based, excludes 0-based
  int32_t end0; // excludes 1-based, includes 0-based
  char buf[255];

  genomeLocus(const char* c, int32_t b, int32_t e) : chrom(c), beg1(b), end0(e) {
    sprintf(buf,"%s:%d-%d",c,b,e);
  }

  // convert [chr]:[beg1]-[end0] string int32_to int32_terval
  // 20:100-110 means [100,110] in 1-based [100,111) in 1-based [99,110) in 0-based
  genomeLocus(const char* region) {
    strcpy(buf,region);
    const char* pcolon = strchr(region,':');
    const char* pminus = strchr(pcolon+1,'-');
    //if ( ( pcolon == NULL ) || ( pminus == NULL ) )
    if ( pcolon == NULL )
      error("Cannot parse %s in genomeLocus::genomeLocus()");
    chrom = std::string(region,0,pcolon-region);
    beg1 = atoi(pcolon+1);
    if ( pminus == NULL ) end0 = INT_MAX;
    else {
      end0 = atoi(pminus+1);
      if ( end0 == 0 ) end0 = INT_MAX;
    }
  }

  const char* toString() const { 
    return buf;
  }

  // compare between genomeLocus
  bool operator< (const genomeLocus& l) const {
    if ( chrom == l.chrom ) {
      if ( beg1 == l.beg1 ) {
	return ( end0 < l.end0 );
      }
      else {
	return ( beg1 < l.beg1 );
      }
    }
    else {
      int32_t n1 = atoi(chrom.c_str());
      int32_t n2 = atoi(l.chrom.c_str());
      if ( ( n1 == 0 ) && ( n2 == 0 ) ) {
	return chrom < l.chrom;
      }
      else if ( ( n1 > 0 ) && ( n2 > 0 ) ) {
	return n1 < n2;
      }
      else { // treat n1 == 0 as infinite
	return ( n1 > 0 ) ? true : false;
      }
    }
  }

  // length
  unsigned long length() const { return end0-beg1+1; }

  bool overlaps(const char* _chrom, int32_t _beg1, int32_t _end0) const {
    if ( chrom == _chrom ) {
      if ( ( beg1 <= _end0 )  && ( _beg1 <= end0 ) ) {
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }    
  }

  // check overlap with other locus
  bool overlaps (const genomeLocus& l) const {
    if ( chrom == l.chrom ) {
      if ( ( beg1 <= l.end0 )  && ( l.beg1 <= end0 ) ) {
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }
  }

  // merge two locus if possible
  bool merge (const genomeLocus& l) {
    if ( chrom == l.chrom ) {
      if ( ( beg1-1 <= l.end0 )  && ( l.beg1-1 <= end0 ) ) {
	if ( l.beg1 < beg1 ) beg1 = l.beg1;
	if ( l.end0 > end0 ) end0 = l.end0;
	return true;
      }
      else {
	return false;
      }
    }
    else {
      return false;
    }
  }

  // check if it contains pos
  bool contains0(const char* chr, int32_t pos0) const { return contains1(chr,pos0+1); }

  bool contains1(const char* chr = NULL, int32_t pos1 = INT_MAX) const {
    if ( ( chr == NULL ) || ( chrom == chr ) ) {
      return ( ( pos1 >= beg1 ) && ( pos1 <= end0 ) );
    }
    else {
      return false;
    }
  }
};

// Collection of genomic locus
class genomeLoci {
 public:
  std::set<genomeLocus> loci;
  std::set<genomeLocus>::iterator it;
  bool overlapResolved;
  int32_t maxLength;

  genomeLoci() : overlapResolved(false), maxLength(0) {}
  genomeLoci(const char* reg) : overlapResolved(false), maxLength(0) {
    add(reg);
    resolveOverlaps();
  }

  // functions for iterating each locus
  void rewind() { it = loci.begin(); }
  bool next() { ++it; return ( it != loci.end() ); }
  bool isend() { return ( it == loci.end() ); }
  const genomeLocus& currentLocus() { return (*it); }

  // check the size 
  bool empty() { return loci.empty(); }
  int32_t numLocus() const { return (int32_t)loci.size(); }

  // add a locus
  bool add(const char* chr, int32_t beg1, int32_t end0) {
    overlapResolved = false;
    if ( end0-beg1+1 > maxLength ) maxLength = end0-beg1+1;
    return loci.insert(genomeLocus(chr,beg1,end0)).second;
  }

  // add a locus
  bool add(const char* region) {
    overlapResolved = false;
    std::pair<std::set<genomeLocus>::iterator, bool> ret = loci.insert(genomeLocus(region));
    int32_t l = ret.first->end0 - ret.first->beg1 + 1;
    if ( l > maxLength ) maxLength = l;
    return ret.second;
  }

  // Resolve overlapping int32_tervals
  int32_t resolveOverlaps() {
    if ( !overlapResolved ) {
      std::set<genomeLocus>::iterator it;
      std::set<genomeLocus>::iterator prev;
      int32_t numMerged = 0;
      for(it = loci.begin(); it != loci.end(); ++it) {
	if ( it != loci.begin() ) {
	  if ( prev->overlaps(*it) ) {
	    // if overlaps, erase both and insert merged one
	    genomeLocus locus = *prev;
	    locus.merge(*it);
	    if ( (int32_t)locus.length() > maxLength ) maxLength = locus.length();
	    loci.erase(it);
	    loci.erase(prev);
	    prev = it = loci.insert(locus).first;
	    ++numMerged;
	  }
	  else {
	    prev = it;
	  }
	}
	else {
	  prev = it;
	}
      }
      overlapResolved = true;
      return numMerged;
    }
    else {
      return 0;
    }
    return 0;
  }

  unsigned long totalLength() const {
    //resolveOverlaps();
    unsigned long sz = 0;
    std::set<genomeLocus>::iterator it2;
    for(it2 = loci.begin(); it2 != loci.end(); ++it2) {
      sz += it2->length();
    }
    return sz;
  }

  bool moveTo(const char* chr = NULL, int32_t pos1 = INT_MAX) {
    if ( it->contains1(chr, pos1) ) return true;
				     
    chr = it->chrom.c_str();
    genomeLocus locus(chr, pos1, pos1);
    it = loci.lower_bound(locus);
    if ( it == loci.begin() ) { // do nothing
      return (it->contains1(chr,pos1));
    }
    else if ( it == loci.end() ) {
      std::set<genomeLocus>::iterator i = it;
      --i;
      if ( i->contains1(chr,pos1) ) { it = i; return true; }
      else { return false; }
    }
    else {
      if ( it->contains1(chr,pos1) ) return true;
      else {
	std::set<genomeLocus>::iterator i = it;
	--i;
	if ( i->contains1(chr,pos1) ) { it = i; return true; }
	else { return false; }
      }
    }
  }

  bool contains1(const char* chr, int32_t pos1) {
    genomeLocus locus(chr, pos1, pos1);
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= pos1 ) ) {
      if ( it2->end0 >= pos1 ) return true;
      ++it2;
    }
    return false;
  }

  bool overlaps(const char* chr, int32_t beg1, int32_t end0) {
    genomeLocus locus(chr, overlapResolved ? beg1 : beg1-maxLength, overlapResolved ? beg1 : beg1-maxLength);
    if ( loci.empty() ) return false;
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= end0 ) ) {
      if ( ( it2->beg1 <= end0 ) && ( beg1 <= it2->end0 ) )
	return true;
      ++it2;
    }
    //notice("%s:%d-%d",it2->chrom.c_str(),it2->beg1,it2->end0);
    return false;
  }

  bool contains(const char* chr, int32_t beg1, int32_t end0) {
    if ( loci.empty() ) return false;

    resolveOverlaps();
    genomeLocus locus(chr, beg1-maxLength, beg1-maxLength);
    std::set<genomeLocus>::iterator it2 = loci.lower_bound(locus);
    if ( it2 != loci.begin() ) --it2;
    if ( it2->chrom != chr ) ++it2;    
    while( it2 != loci.end() && ( it2->chrom == chr ) && ( it2->beg1 <= end0 ) ) {
      if ( ( it2->beg1 <= beg1 ) && ( end0 <= it2->end0 ) )
	return true;
      ++it2;
    }
    return false;    
  }
};

#endif
