TARGET = PALMER
CPP_FILES = $(shell ls *.cpp)
BASE = $(basename $(CPP_FILES))
HTS_CFLAGS = $(shell pkg-config --cflags htslib 2>/dev/null)
HTS_LIBS = $(shell pkg-config --libs htslib 2>/dev/null)
HTS_LIBS ?= -lhts

$(TARGET): $(OBJS)
	g++ -o $(TARGET) $(HTS_CFLAGS) $(BASE).cpp -O3 -w -std=c++11 -pthread $(HTS_LIBS)

clean:
	rm -f PALMER
