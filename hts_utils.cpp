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

#include "hts_utils.h"
extern "C" {
#include "htslib/hfile.h"
}
#include "Error.h"
#include <cassert>

/********
 *General
 ********/

KHASH_MAP_INIT_STR(vdict, bcf_idinfo_t)
typedef khash_t(vdict) vdict_t;

//#ifdef __cplusplus
//extern "C" {
//#endif
//  int ks_resize2(kstring_t*, unsigned long);
//#ifdef __cplusplus
//}
//#endif

//struct faidx_t {
//    BGZF *bgzf;
//    int n, m;
//    char **name;
//    khash_t(s) *hash;
//    enum fai_format_options format;
//};


/**********
 *FAI UTILS
 **********/

/**
 * An alternate sequence fetcher for upper case sequence.
 */
char *faidx_fetch_uc_seq(const faidx_t *fai, const char *c_name, int p_beg_i, int p_end_i, int *len)
{
  char* seq = faidx_fetch_seq(fai, c_name, p_beg_i, p_end_i, len);
  if ( *len > 0 ) {
    for(int i=0; i < *len; ++i)
      if ( isgraph(seq[i]) ) seq[i] = toupper(seq[i]);
  }
  return seq;
    // int l;
    // char c;
    // khiter_t iter;
    // faidx1_t val;
    // char *seq=NULL;

    // // Adjust position
    // iter = kh_get(s, fai->hash, c_name);
    // if(iter == kh_end(fai->hash)) return 0;
    // val = kh_value(fai->hash, iter);
    // if(p_end_i < p_beg_i) p_beg_i = p_end_i;
    // if(p_beg_i < 0) p_beg_i = 0;
    // else if(val.len <= p_beg_i) p_beg_i = val.len - 1;
    // if(p_end_i < 0) p_end_i = 0;
    // else if(val.len <= p_end_i) p_end_i = val.len - 1;

    // // Now retrieve the sequence
    // int ret = bgzf_useek(fai->bgzf, val.offset + p_beg_i / val.line_blen * val.line_len + p_beg_i % val.line_blen, SEEK_SET);
    // if ( ret<0 )
    // {
    //     *len = -1;
    //     fprintf(stderr,"[fai_fetch_seq] Error: fai_fetch failed. (Seeking in a compressed, .gzi unindexed, file?)\n");
    //     return NULL;
    // }
    // l = 0;
    // seq = (char*)malloc(p_end_i - p_beg_i + 2);
    // while ( (c=bgzf_getc(fai->bgzf))>=0 && l < p_end_i - p_beg_i + 1)
    //     if (isgraph(c)) seq[l++] = toupper(c);
    // seq[l] = '\0';
    // *len = l;
    // return seq;
}

/**********
 *HTS UTILS
 **********/

/**
 * Checks file extension for use in writing files.
 */
bool str_ends_with(std::string& file_name, const char* ext)
{
    size_t len = file_name.size();
    const char* suffix = file_name.c_str();
    size_t ext_len = strlen(ext);
    suffix = (len>ext_len) ? suffix + len - ext_len : suffix;
    if (!strcmp(ext, suffix))
    {
        return true;
    }

    return false;
}

/**************
 *BAM HDR UTILS
 **************/

/**
 * Copies contigs found in bam header to bcf header.
 */
void bam_hdr_transfer_contigs_to_bcf_hdr(const bam_hdr_t *sh, bcf_hdr_t *vh)
{
    kstring_t s = {0,0,0};
    for (size_t i=0; i<(size_t)bam_hdr_get_n_targets(sh); ++i)
    {
        s.l = 0;
        ksprintf(&s, "##contig=<ID=%s,length=%d>", bam_hdr_get_target_name(sh)[i], bam_hdr_get_target_len(sh)[i]);
        bcf_hdr_append(vh, s.s);
    }
    if (s.m) free(s.s);
}

/**********
 *BAM UTILS
 **********/

/**
 * Gets the end position of the last mapped base in the read.
 */
int32_t bam_get_end_pos1(bam1_t *s)
{
    int32_t end_pos1 = bam_get_pos1(s);
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        for (int32_t i = 0; i < (int32_t)n_cigar_op; ++i)
        {
            int32_t opchr = bam_cigar_opchr(cigar[i]);
            
            if (opchr=='M' || opchr=='D' || opchr=='N' || opchr=='=' || opchr=='X')
            {
                end_pos1 += bam_cigar_oplen(cigar[i]);
            }
        }
    }
    
    return end_pos1-1;
}

/**
 * Gets the read sequence from a bam record
 */
void bam_get_seq_string(bam1_t *s, kstring_t *seq)
{
    seq->l=0;
    uint8_t* sq = bam_get_seq(s);
    for (uint16_t i = 0; i < bam_get_l_qseq(s); ++i)
    {
        kputc("=ACMGRSVTWYHKDBN"[bam_seqi(sq, i)], seq);
    }
};

/**
 * Gets the base qualities from a bam record, when N is observed, a placeholder value of 0(!, 33 adjusted) is entered
 */
void bam_get_qual_string(bam1_t *s, kstring_t *qual)
{
    qual->l=0;
    uint32_t offset = 0;
    uint8_t* q = bam_get_qual(s);
    for (int32_t i = 0; i < bam_get_l_qseq(s); ++i)
    {
        kputc(q[i-offset] + 33, qual);
    }
};

/**
 * Gets the cigar from a BAM record
 */
void bam_get_cigar_string(bam1_t *s, kstring_t *cigar_string)
{
    cigar_string->l=0;
    int32_t n_cigar_op = bam_get_n_cigar_op(s);
    if (n_cigar_op)
    {
        uint32_t *cigar = bam_get_cigar(s);
        for (int32_t i = 0; i < (int32_t)n_cigar_op; ++i)
        {
            kputw(bam_cigar_oplen(cigar[i]), cigar_string);
            kputc(bam_cigar_opchr(cigar[i]), cigar_string);
        }
    }
}

/**
 * Gets the cigar string from a bam record
 */
void bam_get_cigar_expanded_string(bam1_t *s, kstring_t *cigar_expanded_string)
{
    kstring_t cigar_string = {0,0,0};
    bam_get_cigar_string(s, &cigar_string);

    cigar_expanded_string->l = 0;
    int32_t lastIndex = cigar_string.l;
    int32_t i = 0;
    kstring_t token = {0,0,0};

    if (lastIndex<0)
    {
        return;
    }
    char c;
    bool seenM = false;

    while (i<=lastIndex)
    {
        c = cigar_string.s[i];

        //captures the numeric count
        if (c<'A')
        {
            kputc(c, &token);
        }

        if (c>'A' || i==lastIndex)
        {
            //it is possible for I's to be observed before the first M's in the cigar string
            //in this case, we treat them as 'S'
            if (!seenM)
            {
                if (c=='I')
                {
                    c = 'S';
                }
                else if (c=='M')
                {
                    seenM = true;
                }
            }

            int32_t count = atoi(token.s);
            for (int32_t j=0; j<count; ++j)
                kputc(c, cigar_expanded_string);
            token.l = 0;;
        }

        ++i;
    }

    if (cigar_string.m) free(cigar_string.s);
    if (token.m) free(token.s);
}

/**
 * Gets the base in the read that is mapped to a genomic position.
 * Extracts the read sequence and qualities too.
 */
void bam_get_base_and_qual_and_read_and_qual(bam1_t *srec, uint32_t pos, char& base, char& qual, int32_t& rpos, kstring_t* readseq, kstring_t* readqual)
{
    bam1_core_t *c = &srec->core;
    int32_t rlen = c->l_qseq;
    uint32_t cpos = c->pos; //reference coordinates of the first mapped base
    rpos = 0; //read coordinates

    kstring_t str = {0,0,0};
    base = 'N';
    qual = 0;

    if (c->n_cigar)
    {
        uint32_t *cigar = bam_get_cigar(srec);
        for (uint32_t i = 0; i < c->n_cigar; ++i)
        {
            char op = bam_cigar_opchr(cigar[i]);
            str.l = 0;
            kputw(bam_cigar_oplen(cigar[i]), &str);
            char* stop;
            uint32_t len = strtol(str.s, &stop, 10);
            assert(stop);

	    //fprintf(stderr,"%d%c",bam_cigar_oplen(cigar[i]),op);

            if (op=='M')
            {
                if (pos>=cpos && pos<=cpos+len-1)
                {
                    rpos += pos-cpos;
                    break;
                }

                cpos += len;
                rpos += len;
            }
            else if ( ( op=='D' ) || ( op=='N' ) )
            {
                if (pos>=cpos && pos<=cpos+len-1)
                {
                    rpos = -1;
                    break;
                }

                cpos += len;
            }
            else if (op=='S' || op=='I')
            {
                rpos += len;
            }
        }

        //std::cout << "bpos " << bpos << "\n";
        if (rpos>=0 && rpos<=rlen)
        {
            //sequence
            bam_get_seq_string(srec, readseq);
            base = readseq->s[rpos];

            //qual
            bam_get_qual_string(srec, readqual);
            qual = readqual->s[rpos];
        }
        else
        {
            rpos = BAM_READ_INDEX_NA;
        }
    }
    
    if ( str.s ) free(str.s);

    if ( rpos >= rlen ) {
      rpos = BAM_READ_INDEX_NA;
      base = '.';
    }

    //if ( rand() % 1000 == 0 )
    //fprintf(stderr,", pos = %d, cpos = %d, b = %c, q = %c, rpos=%d, seq=%s, qual = %s\n", pos, cpos, base, qual, rpos, readseq->s, readqual->s);
    //    std::cout << "q: " << s[bpos-1] << " " << q << "\n";
//    for (uint32_t i = 0; i < c->l_qseq; ++i) std::cerr << ((char)(s[i] + 33));
};

/**
 * Prints a bam.
 */
void bam_print(bam_hdr_t *h, bam1_t *s)
{
    const char* chrom = bam_get_chrom(h, s);
    uint32_t pos1 = bam_get_pos1(s);
    kstring_t seq = {0,0,0};
    bam_get_seq_string(s, &seq);
    uint32_t len = bam_get_l_qseq(s);
    kstring_t qual = {0,0,0};
    bam_get_qual_string(s, &qual);
    kstring_t cigar_string = {0,0,0};
    bam_get_cigar_string(s, &cigar_string);
    kstring_t cigar_expanded_string = {0,0,0};
    bam_get_cigar_expanded_string(s, &cigar_expanded_string);
    //uint16_t flag = bam_get_flag(s);
    uint32_t mapq = bam_get_mapq(s);

    std::cerr << "##################" << "\n";
    std::cerr << "chrom-pos: " << chrom << "-" << pos1 << "\n";
    std::cerr << "read     : " << seq.s << "\n";
    std::cerr << "qual     : " << qual.s << "\n";
    std::cerr << "cigar_str: " << cigar_string.s << "\n";
    std::cerr << "cigar    : " << cigar_expanded_string.s << "\n";
    std::cerr << "len      : " << len << "\n";
    std::cerr << "mapq     : " << mapq << "\n";
    std::cerr << "mpos1    : " << bam_get_mpos1(s) << "\n";
    std::cerr << "mtid     : " << bam_get_mtid(s) << "\n";

    if (seq.m) free(seq.s);
    if (qual.m) free(qual.s);
    if (cigar_string.m) free(cigar_string.s);
    if (cigar_expanded_string.m) free(cigar_expanded_string.s);
}

/**************
 *BCF HDR UTILS
 **************/

/**
 * Copies contigs found in bcf header to another bcf header.
 */
void bcf_hdr_transfer_contigs(const bcf_hdr_t *hsrc, bcf_hdr_t *hdest)
{
    vdict_t *d = (vdict_t*)hsrc->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    const char **names = (const char**) calloc(m,sizeof(const char*));
    int len[m];
    khint_t k;

    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        len[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if ( sscanf(kh_val(d, k).hrec[0]->vals[len[tid]],"%d",&j) )
            len[tid] = j;
        names[tid] = kh_key(d,k);
    }

    kstring_t s = {0,0,0};
    for (tid=0; tid<m; tid++)
    {
        s.l = 0;
        ksprintf(&s, "##contig=<ID=%s,length=%d>", names[tid], len[tid]);
        bcf_hdr_append(hdest, s.s);
    }
    if (s.m) free(s.s);
}

/**
 * Extracts sequence length by rid.
 */
int32_t* bcf_hdr_seqlen(const bcf_hdr_t *hdr, int32_t *nseq)
{
    vdict_t *d = (vdict_t*)hdr->dict[BCF_DT_CTG];
    int tid, m = kh_size(d);
    int32_t *len = (int32_t*) malloc(m*sizeof(int32_t));
    khint_t k;

    for (k=kh_begin(d); k<kh_end(d); k++)
    {
        if ( !kh_exist(d,k) ) continue;
        tid = kh_val(d,k).id;
        len[tid] = bcf_hrec_find_key(kh_val(d, k).hrec[0],"length");
        int j;
        if ( sscanf(kh_val(d, k).hrec[0]->vals[len[tid]],"%d",&j) )
            len[tid] = j;
    }

    return len;
}

/**
 * Prints a VCF record to STDERR.
 */
void bcf_print(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    vcf_format(h, v, &s);
    std::cerr << s.s;
    if (s.m) free(s.s);
};

/**
 * Prints a VCF record in compact string representation to STDERR.
 */
void bcf_print_liten(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    bcf_variant2string(h, v, &s);
    std::cerr << s.s << "\n";
    if (s.m) free(s.s);
};

/**
 * Prints a VCF record in compact string representation to STDERR.
 */
void bcf_print_lite(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    bcf_variant2string(h, v, &s);
    std::cerr << s.s;
    if (s.m) free(s.s);
};

/**
 * Prints a VCF record in compact string representation to STDERR with alleles sorted.
 */
void bcf_print_lite_sorted(bcf_hdr_t *h, bcf1_t *v)
{
    kstring_t s = {0,0,0,};
    bcf_variant2string_sorted(h, v, &s);
    std::cerr << s.s;
    if (s.m) free(s.s);
};

/**
 * Reads header of a VCF file and returns the bcf header object.
 * This wraps around vcf_hdr_read from the original htslib to
 * allow for an alternative header file to be read in.
 *
 * this searches for the alternative header saved as <filename>.hdr
 * If the VCF files is BCF, any alternative header is ignored.
 */
bcf_hdr_t *bcf_alt_hdr_read(htsFile *fp)
{
    bcf_hdr_t *h = NULL;

    //check for existence of alternative header
    kstring_t alt_hdr_fn = {0, 0, 0};
    kputs(fp->fn, &alt_hdr_fn);
    kputs(".hdr", &alt_hdr_fn);

    FILE *file = NULL;
    struct stat st;
    if (stat(alt_hdr_fn.s, &st)==0 && st.st_size)
    {
        file = fopen(alt_hdr_fn.s, "r");
    }

    if (fp->format.format ==bcf || !file)
    {
        h = bcf_hdr_read(fp);
    }
    else
    {
        fprintf(stderr,"[I:%s:%d %s] read alternative header for %s\n", __FILE__, __LINE__, __FUNCTION__, fp->fn);
        fclose(file);
        htsFile *alt_hdr = hts_open(alt_hdr_fn.s, "r");
        h = bcf_hdr_read(alt_hdr);
        hts_close(alt_hdr);

        //helps move the pointer to the right place
        bcf_hdr_t *temp_h = bcf_hdr_read(fp);
        bcf_hdr_destroy(temp_h);
    }

    if (alt_hdr_fn.m) free(alt_hdr_fn.s);
    return h;
}

int32_t bcf_hdr_sample_index(bcf_hdr_t* h, const char* id) {
  return bcf_hdr_id2int(h, BCF_DT_SAMPLE, id);
}

const char* bcf_hdr_sample_id(bcf_hdr_t* h, int32_t idx) {
  //return (const char*)(h->id[BCF_DT_SAMPLE][idx]);
  return bcf_hdr_int2id(h, BCF_DT_SAMPLE, idx);
}

/**
 * Get number of samples in bcf header
 */
int32_t bcf_hdr_get_n_sample(bcf_hdr_t *h)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_SAMPLE];
    return kh_size(d);
}

/**
 * Help function for adding a header with a backup tag name.
 * If the <tag> is already present, a new tag is attempted
 * in the format <tag>_1 to <tag>_9.  If <tag>_9 failed,
 * the function will not add any new tag and will return
 * an empty string.
 *
 * Returns the tag that was inserted or updated.
 */
std::string bcf_hdr_append_info_with_backup_naming(bcf_hdr_t *h, std::string tag, std::string number, std::string type, std::string description, bool rename)
{
    if (bcf_hdr_id2int(h,  BCF_DT_ID, tag.c_str())==-1)
    {
        std::string meta_hdr = "##INFO=<ID=" + tag +
                                 ",Number=" + number +
                                 ",Type=" + type +
                                 ",Description=\"" + description + "\">";

        bcf_hdr_append(h, meta_hdr.c_str());
    }
    else
    {
        if (rename)
        {
            std::string new_tag = "";
            for (uint32_t i=0; i<=9; ++i)
            {
                char c = 49+i;
                new_tag = tag +  "_" + c;
                if (bcf_hdr_id2int(h,  BCF_DT_ID, new_tag.c_str())==-1)
                {
                    std::string meta_hdr = "##INFO=<ID=" + new_tag +
                                 ",Number=" + number +
                                 ",Type=" + type +
                                 ",Description=\"" + description + "\">";

                    bcf_hdr_append(h, meta_hdr.c_str());

                    break;
                }

                new_tag = "";
            }

            return new_tag;
        }
        else
        {
            return tag;
        }
    }

    return tag;
}

/**
 *  bcf_get_format_*() - same as bcf_get_info*() above
 *
 *  The function bcf_get_format_string() is a higher-level (slower) variant of bcf_get_format_char().
 *  see the description of bcf_update_format_string() and bcf_update_format_char() above.
 *  Unlike other bcf_get_format__*() functions, bcf_get_format_string() allocates two arrays:
 *  a single block of \0-terminated strings collapsed into a single array and an array of pointers
 *  to these strings. Both arrays must be cleaned by the user.
 *
 *  Returns negative value on error or the number of written values on success.
 *
 *  Example:
 *      int ndst = 0; char **dst = NULL;
 *      if ( bcf_get_format_string(hdr, line, "XX", &dst, &ndst) > 0 )
 *          for (i=0; i<bcf_hdr_nsamples(hdr); i++) printf("%s\n", dst[i]);
 *      free(dst[0]); free(dst);
 *
 *  Example:
 *      int ngt, *gt_arr = NULL, ngt_arr = 0;
 *      ngt = bcf_get_genotypes(hdr, line, &gt_arr, &ngt_arr);
 *
 *  todo: modify to allow direct reading, instead of copying to char*
 */
int32_t bcf_get_format_string_ro(const bcf_hdr_t *hdr, bcf1_t *line, const char *tag, char ***dst, int *ndst)
{
    int i,tag_id = bcf_hdr_id2int(hdr, BCF_DT_ID, tag);
    if ( !bcf_hdr_idinfo_exists(hdr,BCF_HL_FMT,tag_id) ) return -1;    // no such FORMAT field in the header
    if ( bcf_hdr_id2type(hdr,BCF_HL_FMT,tag_id)!=BCF_HT_STR ) return -2;     // expected different type

    if ( !(line->unpacked & BCF_UN_FMT) ) bcf_unpack(line, BCF_UN_FMT);

    for (i=0; i<line->n_fmt; i++)
        if ( line->d.fmt[i].id==tag_id ) break;
    if ( i==line->n_fmt ) return -3;                               // the tag is not present in this record
    bcf_fmt_t *fmt = &line->d.fmt[i];

    int nsmpl = bcf_hdr_nsamples(hdr);
    if ( !*dst )
    {
        *dst = (char**) malloc(sizeof(char*)*nsmpl);
        if ( !*dst ) return -4;     // could not alloc
        (*dst)[0] = NULL;
    }
    int n = (fmt->n+1)*nsmpl;
    if ( *ndst < n )
    {
        (*dst)[0] = (char*) realloc((*dst)[0], n);
        if ( !(*dst)[0] ) return -4;    // could not alloc
        *ndst = n;
    }
    for (i=0; i<nsmpl; i++)
    {
        uint8_t *src = fmt->p + i*fmt->n;
        uint8_t *tmp = (uint8_t*)(*dst)[0] + i*(fmt->n+1);
        memcpy(tmp,src,fmt->n);
        tmp[fmt->n] = 0;
        (*dst)[i] = (char*) tmp;
    }
    return n;
}

/**********
 *BCF UTILS
 **********/

/**
 * n choose r.
 */
uint32_t choose(uint32_t n, uint32_t r)
{
    if (r>n)
    {
        return 0;
    }
    else if (r==n)
    {
        return 1;
    }
    else if (r==0)
    {
        return 1;
    }
    else
    {
        if (r>(n>>1))
        {
            r = n-r;
        }

        uint32_t num = n;
        uint32_t denum = 1;

        for (uint32_t i=1; i<r; ++i)
        {
            num *= n-i;
            denum *= i+1;
        }

        return num/denum;
    }
}

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_ap2g(uint32_t no_allele, uint32_t no_ploidy)
{
    if (no_ploidy==1)
    {
        return no_allele;
    }
    else if (no_ploidy==2)
    {
        return (((no_allele+1)*(no_allele))>>1); ;
    }
    else
    {
        return choose(no_ploidy+no_allele-1, no_allele-1);
    }
}

/**
 * Gets number of genotypes from number of alleles and genotypes.
 */
uint32_t bcf_ag2p(uint32_t no_alleles, uint32_t no_genotypes)
{
    if (no_alleles==2 && no_genotypes==3)
    {
        return 2;
    }
    else if (no_alleles==3 && no_genotypes==6)
    {
        return 2;
    }
    else if (no_alleles==4 && no_genotypes==10)
    {
        return 2;
    }
    else if (no_alleles == no_genotypes)
    {
        return 1;
    }

    uint32_t no_ploidy = 1;
    while (true)
    {
        uint32_t k = choose(no_ploidy+no_alleles-1, no_alleles-1);

        if (k==no_genotypes)
        {
            return no_ploidy;
        }
        else if (k>no_genotypes)
        {
            return 0;
        }

        ++no_ploidy;
    }
}

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_g2i(int32_t* g, uint32_t n)
{
    if (n==1)
    {
        return g[0];
    }
    if (n==2)
    {
        return g[0] + (((g[1]+1)*(g[1]))>>1);
    }
    else
    {
        uint32_t index = 0;
        for (uint32_t i=0; i<n; ++i)
        {
            index += bcf_ap2g(g[i], i+1);
        }
        return index;
    }
}

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_g2i(int32_t g0, int32_t g1)
{
    return g0 + (((g1+1)*(g1))>>1);
}

/**
 * Gets number of genotypes from number of alleles and ploidy.
 */
uint32_t bcf_g2i(std::string genotype)
{
    uint32_t index = 0;
    for (uint32_t i=0; i<genotype.size(); ++i)
    {
        uint32_t allele = genotype.at(i)-65;
        index += bcf_ap2g(allele, i+1);
    }
    return index;
}

/**
 * Gets a string representation of a variant.
 */
void bcf_variant2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;
    kputs(bcf_get_chrom(h, v), var);
    kputc(':', var);
    kputw(bcf_get_pos1(v), var);
    kputc(':', var);
    for (size_t i=0; i<bcf_get_n_allele(v); ++i)
    {
        if (i) kputc('/', var);
        kputs(bcf_get_alt(v, i), var);
    }
}

/**
 * strcmp wrapper for qsort.
 */
int32_t cmpstr(void const *a, void const *b)
{
    char const *aa = *(char const **)a;
    char const *bb = *(char const **)b;

    return strcmp(aa, bb);
}

/**
 * Returns true if a is before b, false otherwise.
 */
bool bcf_is_in_order(bcf1_t *a, bcf1_t *b)
{
    if (bcf_get_rid(a)==bcf_get_rid(b))
    {
        return bcf_get_pos0(a)<=bcf_get_pos0(b);
    }

    return bcf_get_rid(a)<bcf_get_rid(b);
}

/**
 * Gets a sorted string representation of a variant.
 */
void bcf_variant2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;
    kputs(bcf_get_chrom(h, v), var);
    kputc(':', var);
    kputw(bcf_get_pos1(v), var);
    kputc(':', var);

    if (v->n_allele==2)
    {
        kputs(bcf_get_alt(v, 0), var);
        kputc(',', var);
        kputs(bcf_get_alt(v, 1), var);
    }
    else
    {
        char** allele = bcf_get_allele(v);
        char** temp = (char**) malloc((bcf_get_n_allele(v)-1)*sizeof(char*));
        for (size_t i=1; i<v->n_allele; ++i)
        {
            temp[i] = allele[i];
        }
        std::qsort(temp, bcf_get_n_allele(v), sizeof(char*), cmpstr);
        kputs(bcf_get_alt(v, 0), var);
        for (int32_t i=0; i<v->n_allele-1; ++i)
        {
            kputc(',', var);
            kputs(temp[i], var);
        }
        free(temp);
    }
}

/**
 * Gets a string representation of the alleles of a variant.
 */
void bcf_alleles2string(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;

    if (v->n_allele==2)
    {
        kputs(bcf_get_alt(v, 0), var);
        kputc(',', var);
        kputs(bcf_get_alt(v, 1), var);
    }
    else
    {
        char** allele = bcf_get_allele(v);
        for (int32_t i=0; i<v->n_allele; ++i)
        {
            if (i) kputc(',', var);
            kputs(allele[i], var);
        }
    }
}

/**
 * Gets a sorted string representation of the alleles of a variant.
 */
void bcf_alleles2string_sorted(bcf_hdr_t *h, bcf1_t *v, kstring_t *var)
{
    bcf_unpack(v, BCF_UN_STR);
    var->l = 0;

    if (v->n_allele==2)
    {
        kputs(bcf_get_alt(v, 0), var);
        kputc(',', var);
        kputs(bcf_get_alt(v, 1), var);
    }
    else
    {
        char** allele = bcf_get_allele(v);
        char** temp = (char**) malloc((bcf_get_n_allele(v)-1)*sizeof(char*));
        for (int32_t i=1; i<v->n_allele; ++i)
        {
            temp[i-1] = allele[i];
        }

        std::qsort(temp, bcf_get_n_allele(v)-1, sizeof(char*), cmpstr);

        kputs(bcf_get_alt(v, 0), var);
        for (int32_t i=0; i<v->n_allele-1; ++i)
        {
            kputc(',', var);
            kputs(temp[i], var);
        }

        free(temp);
    }
}

/**
 * Get chromosome name
 */
const char* bcf_get_chrom(bcf_hdr_t *h, bcf1_t *v)
{
    if (v->rid >= h->n[BCF_DT_CTG])
    {
      error("[E:%s:%d %s] [E:%s:%d %s] rid '%d' does not have an associated contig defined in the header.  Try tabix workaround or just add the header.\n",__FILE__,__LINE__,__FUNCTION__, __FILE__, __LINE__, __FUNCTION__, v->rid);
      //exit(1);
      //return NULL;
    }
    else if ( v->rid < 0 ) return NULL;
    
    return h->id[BCF_DT_CTG][v->rid].key;
}

/**
 *Set chromosome name.
 */
void bcf_set_chrom(bcf_hdr_t *h, bcf1_t *v, const char* chrom)
{
    vdict_t *d = (vdict_t*)h->dict[BCF_DT_CTG];
    khint_t k = kh_get(vdict, d, chrom);
    if (k == kh_end(d))
    {
      error("[E:%s:%d %s] [E:%s:%d %s] contig '%s' is not defined in the header\n",__FILE__,__LINE__,__FUNCTION__, __FILE__, __LINE__, __FUNCTION__, chrom);
        kstring_t contig = {0,0,0};
        ksprintf(&contig, "##contig=<ID=%s,length=2147483647>", chrom);
        bcf_hdr_append(h, contig.s);
        if (contig.m) free(contig.s);
        k = kh_get(vdict, d, chrom);
    }
    v->rid = kh_val(d, k).id;
};

/**
 * Set id.
 */
void bcf_set_id(bcf1_t *v, char* id)
{
    if (v->d.id)
    {
        free(v->d.id);
    }

    v->d.id = strdup(id);
};

void hprintf(htsFile* fp, const char * msg, ...) {
  va_list ap;

  va_start(ap, msg);

  kstring_t tmp = {0,0,0};
  kvsprintf(&tmp, msg, ap);

  int ret;
  if ( fp->format.compression != no_compression )
    ret = bgzf_write(fp->fp.bgzf, tmp.s, tmp.l);
  else
    ret = hwrite(fp->fp.hfile, tmp.s, tmp.l);

  free(tmp.s);

  if ( ret < 0 ) {
    error("[E:%s:%d %s] [E:%s:%d %s] hprintf failed. Aborting..",__FILE__,__LINE__,__FUNCTION__, __FILE__, __LINE__, __FUNCTION__);
  }

  va_end(ap);
}

void parse_intervals(std::vector<GenomeInterval>& intervals, std::string interval_list, std::string interval_string)
{
    intervals.clear();
    std::map<std::string, uint32_t> m;

    if (interval_list!="")
    {
        htsFile *file = hts_open(interval_list.c_str(), "r");
        if (file)
        {
            kstring_t *s = &file->line;
            while (hts_getline(file, '\n', s)>=0)
            {
                std::string ss = std::string(s->s);
                if (m.find(ss)==m.end())
                {
                    m[ss] = 1;
                    GenomeInterval interval(ss);
                    intervals.push_back(interval);
                }
            }
            hts_close(file);
        }
    }

    std::vector<std::string> v;
    if (interval_string!="")
        split(v, ",", interval_string);

    for (size_t i=0; i<v.size(); ++i)
    {
        if (m.find(v[i])==m.end())
        {
            m[v[i]] = 1;
            GenomeInterval interval(v[i]);
            intervals.push_back(interval);
        }
    }
}

std::string bam_hdr_get_sample_name(bam_hdr_t* hdr) {
  if ( !hdr )
    error("[E:%s:%d %s] [E:%s:%d %s] Failed to read the BAM header",__FILE__,__LINE__,__FUNCTION__,__FILE__, __LINE__, __FUNCTION__);

  char *ptext = strdup(hdr->text);  
  const char *p = ptext; 
  const char *q, *r;
  int32_t n = 0;
  std::string sm;
  while( ( q = strstr(p, "@RG" ) ) != 0 ) {
    p = q + 3;
    r = q = 0;
    if ( ( q = strstr(p, "\tID:" ) ) != 0 ) q += 4;
    if ( ( r = strstr(p, "\tSM:" ) ) != 0 ) r += 4;
    if ( r && q ) {
      char *u, *v;
      for (u = (char*)q; *u && *u != '\t' && *u != '\n'; ++u);
      for (v = (char*)r; *v && *v != '\t' && *v != '\n'; ++v);
      *u = *v = '\0';
      if ( sm.empty() )
	sm = r;
      else if ( sm.compare(r) != 0 ) {
	error("[E:%s:%d %s] [E:%s:%d %s] Multiple sample IDs are included in one BAM file - %s, %s",__FILE__,__LINE__,__FUNCTION__, __FILE__, __LINE__, __FUNCTION__, sm.c_str(), r);
	//abort();
      }
    }
    else break;
    p = q > r ? q : r;
    ++n;
  }
  if ( sm.empty() ) {
    warning("[W:%s:%d %s] Sample ID information cannot be found",__FILE__,__LINE__,__FUNCTION__);
  }
  free(ptext);
  return sm;
}

int32_t bam_get_unclipped_start(bam1_t* b) {
  uint32_t* cigar = bam_get_cigar(b);
  bam1_core_t* c = &b->core;
  int32_t i, y;
  for(i = y =0; i < (int32_t)c->n_cigar; ++i) {
    int l = cigar[i]>>4, op = cigar[i]&0xf;
    if ( ( op == BAM_CSOFT_CLIP ) || ( op == BAM_CHARD_CLIP ) )
      y += l;
    else
      return ( c->pos - y );
  }
  return ( c->pos - y );
}

int32_t bam_get_unclipped_end(bam1_t* b) {
  uint32_t* cigar = bam_get_cigar(b);
  bam1_core_t* c = &b->core;
  int32_t i, y;
  for(i = y = 0; i < (int32_t)c->n_cigar; ++i) {
    int l = cigar[i]>>4, op = cigar[i]&0xf;
    switch( op ) {
    case BAM_CMATCH:
    case BAM_CEQUAL:
    case BAM_CDIFF:
    case BAM_CDEL:
    case BAM_CREF_SKIP:      
    case BAM_CSOFT_CLIP:
    case BAM_CHARD_CLIP:      
      y += l;
      //case BAM_CINS:
      //case BAM_CPAD:
      //case BAM_CBACK:
    }
  }
  return ( c->pos + y );  
}

int32_t bam_get_clipped_end(bam1_t* b) {
  uint32_t* cigar = bam_get_cigar(b);
  bam1_core_t* c = &b->core;
  int32_t i, y;
  for(i = y = 0; i < (int32_t)c->n_cigar; ++i) {
    int l = cigar[i]>>4, op = cigar[i]&0xf;
    switch( op ) {
    case BAM_CMATCH:
    case BAM_CEQUAL:
    case BAM_CDIFF:
    case BAM_CDEL:
    case BAM_CREF_SKIP:      
      //case BAM_CSOFT_CLIP:
      //case BAM_CHARD_CLIP:      
      y += l;
      //case BAM_CINS:
      //case BAM_CPAD:
      //case BAM_CBACK:
    }
  }
  return ( c->pos + y );  
}


bool same_hrecs(bcf_hdr_t* dst_hdr, bcf_hrec_t* dst, bcf_hdr_t* src_hdr, bcf_hrec_t* src) {
  // Check that both records are of the same type. The bcf_hdr_id2length
  // macro cannot be used here because dst header is not synced yet.
  vdict_t *d_src = (vdict_t*)src_hdr->dict[BCF_DT_ID];
  vdict_t *d_dst = (vdict_t*)dst_hdr->dict[BCF_DT_ID];
  khint_t k_src  = kh_get(vdict, d_src, src->vals[0]);
  khint_t k_dst  = kh_get(vdict, d_dst, src->vals[0]);
  if ( (kh_val(d_src,k_src).info[src->type]>>8 & 0xf) != (kh_val(d_dst,k_dst).info[dst->type]>>8 & 0xf) ) {
    warning("Warning: trying to combine \"%s\" tag definitions of different lengths\n", src->vals[0]);
    return false;
  }
  if ( (kh_val(d_src,k_src).info[src->type]>>4 & 0xf) != (kh_val(d_dst,k_dst).info[dst->type]>>4 & 0xf) ) {
    warning("Warning: trying to combine \"%s\" tag definitions of different types\n", src->vals[0]);
    return false;
  }
  return true;
}

char *samfaipath(const char *fn_ref)
{
    char *fn_list = 0;
    if (fn_ref == 0) return 0;
    fn_list = (char*) calloc(strlen(fn_ref) + 5, 1);
    strcat(strcpy(fn_list, fn_ref), ".fai");
    if (access(fn_list, R_OK) == -1) { // fn_list is unreadable
        if (access(fn_ref, R_OK) == -1) {
            fprintf(stderr, "[samfaipath] fail to read file %s.\n", fn_ref);
        } else {
            if (fai_build(fn_ref) == -1) {
                fprintf(stderr, "[samfaipath] fail to build FASTA index.\n");
                free(fn_list); fn_list = 0;
            }
        }
    }
    return fn_list;
};

/*
// Minimal sanitisation of a header to ensure.
// - null terminated string.
// - all lines start with @ (also implies no blank lines).
//
// Much more could be done, but currently is not, including:
// - checking header types are known (HD, SQ, etc).
// - syntax (eg checking tab separated fields).
// - validating n_targets matches @SQ records.
// - validating target lengths against @SQ records.
bam_hdr_t *sam_hdr_sanitise(bam_hdr_t *h) {
    if (!h)
        return NULL;

    // Special case for empty headers.
    if (h->l_text == 0)
        return h;

    uint32_t i, lnum = 0;
    char *cp = h->text, last = '\n';
    for (i = 0; i < h->l_text; i++) {
        // NB: l_text excludes terminating nul.  This finds early ones.
        if (cp[i] == 0)
            break;

        // Error on \n[^@], including duplicate newlines
        if (last == '\n') {
            lnum++;
            if (cp[i] != '@') {
                error("Malformed SAM header at line %u", lnum);
                bam_hdr_destroy(h);
                return NULL;
            }
        }

        last = cp[i];
    }

    if (i < h->l_text) { // Early nul found.  Complain if not just padding.
        uint32_t j = i;
        while (j < h->l_text && cp[j] == '\0') j++;
        if (j < h->l_text)
            warning("Unexpected NUL character in header. Possibly truncated");
    }

    // Add trailing newline and/or trailing nul if required.
    if (last != '\n') {
        warning("Missing trailing newline on SAM header. Possibly truncated");

        if (h->l_text == UINT32_MAX) {
            error("No room for extra newline");
            bam_hdr_destroy(h);
            return NULL;
        }

        if (i >= h->l_text - 1) {
	  cp = (char*)realloc(h->text, (size_t) h->l_text+2);
            if (!cp) {
                bam_hdr_destroy(h);
                return NULL;
            }
            h->text = cp;
        }
        cp[i++] = '\n';

        // l_text may be larger already due to multiple nul padding
        if (h->l_text < i)
            h->l_text = i;
        cp[h->l_text] = '\0';
    }

    return h;
}
*/

/*
bam_hdr_t* bam_hdr_merge(std::vector<bam_hdr_t*> hdrs) {
  bam_hdr_t* merged_hdr* merged_hdr =

merged_header_t *merged_hdr;

    merged_hdr = calloc(1, sizeof(*merged_hdr));
    if (merged_hdr == NULL) return NULL;

    merged_hdr->targets_sz   = 16;
    merged_hdr->target_name = malloc(merged_hdr->targets_sz
                                     * sizeof(*merged_hdr->target_name));
    if (NULL == merged_hdr->target_name) goto fail;

    merged_hdr->target_len = malloc(merged_hdr->targets_sz
                                    * sizeof(*merged_hdr->target_len));
    if (NULL == merged_hdr->target_len) goto fail;

    merged_hdr->sq_tids = kh_init(c2i);
    if (merged_hdr->sq_tids == NULL) goto fail;

    merged_hdr->rg_ids = kh_init(cset);
    if (merged_hdr->rg_ids == NULL) goto fail;

    merged_hdr->pg_ids = kh_init(cset);
    if (merged_hdr->pg_ids == NULL) goto fail;

    return merged_hdr;
    
}
*/
