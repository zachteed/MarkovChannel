# SRC := src/*.cpp

MKLROOT = /opt/intel/mkl

INCLUDE_DIRS = -I include -I proto -I$(MKLROOT)/include

CC = gcc

PROTOC = protoc

CFLAGS = $(INCLUDE_DIRS) -m64 -D USE_MKL

LFLAGS = -L$(MKLROOT)/lib/intel64 -L/opt/intel/lib/intel64

LIBS = -lmkl_intel_lp64 -lmkl_core -lmkl_intel_thread -liomp5 -lpthread -lm -ldl -lprotobuf -lstdc++

SRCS = proto/MarkovChannel.pb.cc src/math_functions.cpp src/graph_functions.cpp src/Model.cpp src/tests/model.cpp

OBJS = $(SRCS:.c=.o)

MAIN = MarkovChannel

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean

all:    $(MAIN)

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)


# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)

proto/MarkovChannel.pb.cc:
	$(PROTOC) -I=proto --cpp_out=proto proto/MarkovChannel.proto

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it
