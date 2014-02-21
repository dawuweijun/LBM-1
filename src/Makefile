TARGET=lbm
SRC=main.cpp 
HEADERS=$(wildcard %.h,%.hpp)
OBJS=$(patsubst %.cpp,%.o,$(SRC))
DEPENDS=$(patsubst %.cpp,vstream%.d,$(SRC))
LIBS= 
INCLUDE= 
FLAGS= -Wall

CC=g++ -std=c++0x
CC+= $(FLAGS) $(INCLUDE) $(LIBS)

.PHONY : release debug clean depend
release : CC+= -fopenmp -O2
release : all
debug : CC+= -g -O0
debug : all
clean:
	-$(RM) $(TARGET) $(OBJS) $(DEPENDS)

%.d: %.cpp %.hpp
	@set -e; $(CC) -MM $(CPPFLAGS) $< \
		| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@
	-include $(DEPENDS)


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS)

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) -c $<
