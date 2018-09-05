TARGET = PALMER
CPP_FILES = $(shell ls *.cpp)
BASE = $(basename $(CPP_FILES))

$(TARGET):$(OBJS)
	-rm -f $@
	g++ -o $(TARGET) $(BASE).cpp
