PROG=lbm
SRC=main.cpp 
HEADERS=$(wildcard %.h,%.hpp)
OBJS=$(patsubst %.cpp,%.o,$(SRC))
DEPENDS=$(patsubst %.cpp,vstream%.d,$(SRC))
LIBS= 
INCLUDE= 
FLAGS= -fopenmp -O2


CC=g++ -std=c++0x

CC+= $(FLAGS) $(INCLUDE) $(LIBS)

$(PROG): $(OBJS)
	$(CC) -o $@ $(OBJS)

.cpp.o: $(HEADERS)
	$(CC) $(CFLAGS) -c $<

.PHONY : clean depend
clean:
	-$(RM) $(PROG) $(OBJS) $(DEPENDS)

%.d: %.cpp %.hpp
	@set -e; $(CC) -MM $(CPPFLAGS) $< \
		| sed 's/\($*\)\.o[ :]*/\1.o $@ : /g' > $@; \
		[ -s $@ ] || rm -f $@
	-include $(DEPENDS)
