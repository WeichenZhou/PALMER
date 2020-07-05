TARGET 			= 	PALMER
CPP_FILES 		= 	$(shell ls *.cpp)
BASE 			= 	$(basename $(CPP_FILES))

CPP_FLAGS 		= 	-I. -fpermissive -lpthread
C_FLAGS   		= 	-g -w -O3

SAMTOOLS_DIR		= 	samtools
HTSLIB_DIR 			= 	htslib

HTSLIB 			= 	$(HTSLIB_DIR)/libhts.a
HTSLIB_LIB 		= 	$(HTSLIB) $(HTSLIB_static_LIBS)
HTSLIB_LDFLAGS 	= 	$(HTSLIB_static_LDFLAGS)
HTSLIB_CPPFLAGS =	-I$(HTSLIB_DIR)

SAMTOOLS_CPPFLAGS 		=	-I$(SAMTOOLS_DIR) -L$(HTSLIB_DIR)
SAMTOOLS_OBJS			=	$(SAMTOOLS_DIR)/sample.o \
							$(SAMTOOLS_DIR)/sam.o \
							$(SAMTOOLS_DIR)/sam_utils.o \
							$(SAMTOOLS_DIR)/bam.o \
							$(SAMTOOLS_DIR)/bam_plbuf.o \
							$(SAMTOOLS_DIR)/bedidx.o \
							$(SAMTOOLS_DIR)/sam_opts.o

							# $(SAMTOOLS_DIR)/sam_view.o \

ALL_CPPFLAGS 	= 	$(SAMTOOLS_CPPFLAGS) $(HTSLIB_CPPFLAGS) $(CPP_FLAGS) $(C_FLAGS)

all: $(TARGET) samview.o

HTSDIR=$(HTSLIB_DIR)
include $(HTSLIB_DIR)/htslib_static.mk
include $(HTSLIB_DIR)/htslib.mk

OBJS = samview.o $(HTSLIB_LIB) 

$(TARGET): $(OBJS)
	g++ -o $(TARGET) $(BASE).cpp $(OBJS) $(SAMTOOLS_OBJS) -std=c++11 $(ALL_CPPFLAGS) 

samview.o: extension/samview.cpp $(HTSLIB_PUBLIC_HEADERS) 
	g++ $(ALL_CPPFLAGS) -c $^ 

clean :
	rm -rf $(TARGET) *.o workdir/*
