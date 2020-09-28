TARGET 			= 	PALMER
CPP_FILES 		= 	$(shell ls *.cpp)
BASE 			= 	$(basename $(CPP_FILES))

CPP_FLAGS 		= 	-I. -fpermissive -lpthread -std=c++11
# C_FLAGS   		= 	-g -Wall -O2
C_FLAGS   		= 	-g -w -O2

CC				= 	g++

SAMTOOLS_DIR	= 	samtools
HTSLIB_DIR 		= 	htslib

HTSLIB 			= 	$(HTSLIB_DIR)/libhts.a
HTSLIB_LIB 		= 	$(HTSLIB) $(HTSLIB_static_LIBS)
HTSLIB_LDFLAGS 	= 	$(HTSLIB_static_LDFLAGS)
HTSLIB_CPPFLAGS =	-I$(HTSLIB_DIR)

SAMTOOLS_CPPFLAGS 		=	-I$(SAMTOOLS_DIR) -L$(HTSLIB_DIR)
SAMTOOLS_OBJS			=	$(SAMTOOLS_DIR)/sam.o \
							$(SAMTOOLS_DIR)/bam.o \
							$(SAMTOOLS_DIR)/bam_plbuf.o \
							$(SAMTOOLS_DIR)/bedidx.o \
							$(SAMTOOLS_DIR)/libst.a 

ALL_CPPFLAGS 	= 	$(SAMTOOLS_CPPFLAGS) $(HTSLIB_CPPFLAGS) $(CPP_FLAGS) $(C_FLAGS) 

all: $(TARGET)

HTSDIR=$(HTSLIB_DIR)
include $(HTSLIB_DIR)/htslib_vars.mk
include $(HTSLIB_DIR)/htslib_static.mk
include $(HTSLIB_DIR)/htslib.mk

OBJS = samview.o samfaidx.o $(HTSLIB_LIB)

samtools_h = samtools/samtools.h $(htslib_hts_defs_h) $(htslib_sam_h)

$(TARGET): $(OBJS)
	$(CC) -o $(TARGET) $(BASE).cpp $(OBJS) $(SAMTOOLS_OBJS) $(ALL_CPPFLAGS) 

samview.o: extension/samview.cpp $(HTSLIB_PUBLIC_HEADERS) 
	$(CC) $(ALL_CPPFLAGS) -c $^ 
samfaidx.o: extension/samfaidx.cpp $(HTSLIB_PUBLIC_HEADERS) $(samtools_h)
	$(CC) $(ALL_CPPFLAGS) -c $^ 

clean :
	rm -rf $(TARGET) *.o 0000/*
