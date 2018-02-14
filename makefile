TARGET = PALMER
CPP_FILES = $(shell ls *.cpp)
BASE = $(basename $(CPP_FILES))
OBJS = $(addsuffix .o, $(addprefix obj/,$(BASE)))

$(TARGET):$(OBJS)
	-rm -f $@
	g++ -o $(TARGET) $(OBJS)
 
obj/%.o:%.cpp
	@if test ! -d "obj"; then\
		mkdir -p obj;\
	fi;
	g++ -c -o $@ $<

clean:
	-rm -rf obj/
