DEBUG = -O2

FRAGREP_TARGET = fragrep
FASPLIT_TARGET = fasplit
ALIMATCH_TARGET = alimatch
SNGLNALN_TARGET = snglnaln
ALN2PATTERN_TARGET = aln2pattern
PATTERN2EPS_TARGET = pattern2eps

SEQUTIL_PATH = ..

MAKE = make
CC = g++
OPTIMIZE = -O2
VERSION = 2

ifdef use-gsl
GSL_LIB = -lgsl
GSL_FLAG = -DUSE_GSL
else
GSL_LIB = 
GSL_FLAG =
endif

CXXFLAGS = -I$(SEQUTIL_PATH)

OBJECTS = fragrep.o omit_pattern_matching.o 

SEQUTIL_OBJECTS = fasta.o alignment_reader.o bbq_util.o bbq_tokenizer.o position_weight_matrix.o aln_util.o utils.o

ALIMATCH_OBJECTS = alimatch.o 

SNGLNALN_OBJECTS = snglnaln.o 

ALN2PATTERN_OBJECTS = aln2pattern.o pattern_matching.o eps_output.o

PATTERN2EPS_OBJECTS = pattern2eps.o pattern_matching.o eps_output.o

FRAGEXPAND_OBJECTS = fragexpand.o bbq_util.o

FASPLIT_OBJECTS = fasta.o bbq_util.o

all: 
	$(MAKE) fragrep DEBUG=-O2  
	$(MAKE) alimatch DEBUG=-O2  
	$(MAKE) snglnaln DEBUG=-O2 
	$(MAKE) aln2pattern DEBUG=-O2  
	$(MAKE) pattern2eps DEBUG=-O2  	
	$(MAKE) fragexpand DEBUG=-O2  
	$(MAKE) fasplit DEBUG=-O2  

debug: 
	$(MAKE) fragrep DEBUG="-g -D_STLP_DEBUG"
	$(MAKE) alimatch DEBUG="-g -D_STLP_DEBUG"
	$(MAKE) snglnaln DEBUG="-g -D_STLP_DEBUG"
	$(MAKE) aln2pattern DEBUG="-g -D_STLP_DEBUG"
	$(MAKE) pattern2eps DEBUG="-g -D_STLP_DEBUG"
	$(MAKE) fragexpand DEBUG="-g -D_STLP_DEBUG"
	$(MAKE) fasplit DEBUG="-g -D_STLP_DEBUG"

release: 
	$(MAKE) clean
	$(MAKE) all SEQUTIL_PATH=.


fragrep:  $(OBJECTS) $(SEQUTIL_OBJECTS)
	$(CC) -o $(FRAGREP_TARGET) $(OBJECTS) $(SEQUTIL_OBJECTS) $(LINKOPTS) $(GSL_LIB)

snglnaln:  $(SNGLNALN_OBJECTS) $(SEQUTIL_OBJECTS)
	$(CC) -o $(SNGLNALN_TARGET) $(SNGLNALN_OBJECTS) $(SEQUTIL_OBJECTS) $(LINKOPTS) $(GSL_LIB)

template.tps: template.eps
	./make_ps_template.sh

aln2pattern:  $(ALN2PATTERN_OBJECTS) $(SEQUTIL_OBJECTS)
	$(CC) -o $(ALN2PATTERN_TARGET) $(ALN2PATTERN_OBJECTS) $(SEQUTIL_OBJECTS) $(LINKOPTS)  $(GSL_LIB)

pattern2eps:  $(PATTERN2EPS_OBJECTS) $(SEQUTIL_OBJECTS)
	$(CC) -o $(PATTERN2EPS_TARGET) $(PATTERN2EPS_OBJECTS) $(SEQUTIL_OBJECTS) $(LINKOPTS)  $(GSL_LIB)

alimatch:  $(ALIMATCH_OBJECTS) $(SEQUTIL_OBJECTS)
	$(CC) -o $(ALIMATCH_TARGET) $(ALIMATCH_OBJECTS) $(SEQUTIL_OBJECTS) $(LINKOPTS) $(GSL_LIB)

fragexpand:  fragexpand.cpp $(FRAGEXPAND_OBJECTS)
	$(CC) -I $(SEQUTIL_PATH) -o fragexpand -O2 fragexpand.cpp $(SEQUTIL_OBJECTS)  $(GSL_LIB)

fasplit:  fasplit.cpp $(FASPLIT_OBJECTS)
	$(CC) -I $(SEQUTIL_PATH) -o fasplit -O2 fasplit.cpp $(FASPLIT_OBJECTS) $(GSL_LIB)

aln2pattern.o: aln2pattern.cpp template.tps
	$(CC) -c $(CXXFLAGS) $(DEBUG) aln2pattern.cpp

pattern2eps.o: pattern2eps.cpp template.tps
	$(CC) -c $(CXXFLAGS) $(DEBUG) pattern2eps.cpp

bbq_tokenizer.o: $(SEQUTIL_PATH)/sequtil/bbq_tokenizer.h $(SEQUTIL_PATH)/sequtil/bbq_tokenizer.cpp
	$(CC) -o bbq_tokenizer.o -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/bbq_tokenizer.cpp

bbq_util.o: $(SEQUTIL_PATH)/sequtil/bbq_util.h $(SEQUTIL_PATH)/sequtil/bbq_util.cpp
	$(CC) -o bbq_util.o -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/bbq_util.cpp

position_weight_matrix.o: $(SEQUTIL_PATH)/sequtil/position_weight_matrix.h $(SEQUTIL_PATH)/sequtil/position_weight_matrix.cpp
	$(CC) -o position_weight_matrix.o -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/position_weight_matrix.cpp $(GSL_FLAG)

fasta.o: $(SEQUTIL_PATH)/sequtil/fasta.h $(SEQUTIL_PATH)/sequtil/fasta.cpp
	$(CC) -o fasta.o -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/fasta.cpp

alignment_reader.o: $(SEQUTIL_PATH)/sequtil/alignment_reader.hpp $(SEQUTIL_PATH)/sequtil/alignment_reader.cpp
	$(CC) -o alignment_reader.o -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/alignment_reader.cpp

aln_util.o: $(SEQUTIL_PATH)/sequtil/aln_util.h $(SEQUTIL_PATH)/sequtil/aln_util.cpp
	$(CC) -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/aln_util.cpp

alimatch.o: alimatch.cpp 
	$(CC) -c $(CXXFLAGS) $(DEBUG) alimatch.cpp

utils.o: $(SEQUTIL_PATH)/sequtil/utils.h $(SEQUTIL_PATH)/sequtil/utils.cpp
	$(CC) -c $(CXXFLAGS) $(DEBUG) $(SEQUTIL_PATH)/sequtil/utils.cpp

snglnaln.o: snglnaln.cpp
	$(CC) -c $(CXXFLAGS) $(DEBUG) snglnaln.cpp

fragexpand.o: fragexpand.cpp
	$(CC) -c $(CXXFLAGS) $(DEBUG) fragexpand.cpp

fasplit.o: fasplit.cpp
	$(CC) -c $(CXXFLAGS) $(DEBUG) fasplit.cpp

fragrep.o: fragrep.cpp
	$(CC) -c $(CXXFLAGS) $(CXXFLAGS) $(DEBUG) fragrep.cpp 

pattern_matching.o: pattern_matching.cpp pattern_matching.h
	$(CC) -c $(CXXFLAGS) $(DEBUG) pattern_matching.cpp

omit_pattern_matching.o: omit_pattern_matching.cpp omit_pattern_matching.h
	$(CC) -c  $(CXXFLAGS) $(DEBUG) omit_pattern_matching.cpp

eps_output.o: template.tps eps_output.cpp eps_output.h
	$(CC) -c  $(CXXFLAGS) $(DEBUG) eps_output.cpp

clean :
	rm -f fragrep
	rm -f aln2pattern
	rm -f alimatch
	rm -f snglnaln
	rm -f fragexpand
	rm -f fasplit
	rm -f template.tps
	rm -f fragrep-$(VERSION).tar.gz
	rm -f fragrep-$(VERSION).tar
	rm -rf fragrep-$(VERSION)
	rm -rf *.o

package: 
	rm -rf fragrep-$(VERSION)
	mkdir fragrep-$(VERSION)
	cp *.cpp fragrep-$(VERSION)
	cp *.h fragrep-$(VERSION)
	cp template.eps fragrep-$(VERSION)
	cp make_ps_template.sh fragrep-$(VERSION)	
	mkdir fragrep-$(VERSION)/man
	cp man/aln2pattern.1 fragrep-$(VERSION)/man
	cp man/pattern2eps.1 fragrep-$(VERSION)/man
	cp man/fragrep.1 fragrep-$(VERSION)/man
	mkdir fragrep-$(VERSION)/sequtil
	cp $(SEQUTIL_PATH)/sequtil/*.h fragrep-$(VERSION)/sequtil
	cp $(SEQUTIL_PATH)/sequtil/*.hpp fragrep-$(VERSION)/sequtil
	cp $(SEQUTIL_PATH)/sequtil/*.cpp fragrep-$(VERSION)/sequtil
	cp Makefile fragrep-$(VERSION)
	cp -r doc fragrep-$(VERSION)
	cp AUTHORS fragrep-$(VERSION)
	cp COPYING fragrep-$(VERSION)
	cp INSTALL fragrep-$(VERSION)
	cp README fragrep-$(VERSION)
	tar cvvf fragrep-$(VERSION).tar fragrep-$(VERSION)/
	gzip fragrep-$(VERSION).tar
	rm -rf fragrep-$(VERSION)/
