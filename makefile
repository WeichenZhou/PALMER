TARGET = PALMER
CPP_FILES = $(shell ls *.cpp)
BASE = $(basename $(CPP_FILES))

$(TARGET):$(OBJS)
	g++ -o $(TARGET) $(BASE).cpp -O3 -w -std=c++11

clean :
	rm PALMER
