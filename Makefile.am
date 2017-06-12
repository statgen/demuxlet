## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

AM_CPPFLAGS = -fopenmp -I ../../htslib/ -I ../../htslib/htslib -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 

noinst_HEADERS = \
	Constant.h
	bcf_variant_key.h bcf_filter_arg.h

bin_PROGRAMS = demuxlet

## vcfast files
demulxet_SOURCES = \
	Error.cpp PhredHelper.cpp params.cpp \
	hts_utils.cpp bcf_ordered_reader.cpp interval.cpp interval_tree.cpp utils.cpp genome_interval.cpp \
	bcf_ordered_writer.cpp reference_sequence.cpp \
	bcf_chunked_reader.cpp genomeChunk.cpp \
	bcf_filtered_reader.cpp \
	bam_ordered_reader.cpp \
	filter.cpp \
	tsv_reader.cpp \
	cmd_cram_demuxlet.cpp

demuxlet_LDADD = ../../htslib/libhts.a -lpthread -llzma -lz -lbz2 -lgomp -lcurl -lcrypto
