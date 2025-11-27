TARGET = PALMER
CPP_FILES = $(shell ls *.cpp)
BASE = $(basename $(CPP_FILES))
HTS_CFLAGS = $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LIBS = $(shell pkg-config --libs htslib 2>/dev/null)
HTS_LIBS ?= -lhts

ifeq ($(strip $(HTS_CFLAGS)),)
$(error "htslib headers not found. Please install htslib and ensure pkg-config can locate it (see README.md).")
endif

$(TARGET): $(OBJS)
	g++ -o $(TARGET) $(HTS_CFLAGS) $(BASE).cpp -O3 -w -std=c++17 -pthread -lstdc++fs $(HTS_LIBS)

clean:
	rm -f PALMER
