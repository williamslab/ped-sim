CPPSRCS= main.cc cmdlineopts.cc cointerfere.cc
CSRCS= memcpy.c
OBJS= $(patsubst %.cc,%.o,$(CPPSRCS)) $(patsubst %.c,%.o,$(CSRCS))
EXEC= ped-sim

GPP = g++
GCC = gcc
DEFINES= 
CFLAGS = -Wall $(DEFINES)
CPPFLAGS = -std=c++11 $(CFLAGS)
ifdef DEBUG           # to use run `make DEBUG=1`
  CFLAGS += -g
else
  CFLAGS += -O2
endif

# profiling: run program normally then do:
# (note: I haven't read about the options below, I just found them in a
#  Dr. Dobb's article.)
#     `gprof -b -z hapi gmon.out`
ifdef PROFILE       # to use run `make PROFILE=1
  CFLAGS += -pg
endif

LIBS = 

# dependency variables / commands
DEPDIR = .deps
df = $(DEPDIR)/$(*F)

all: $(EXEC)

$(EXEC): $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS)

# for minimal dependencies on libraries:
distribute: $(OBJS) $(HEADERS)
	$(GPP) -o $(EXEC) $(OBJS) $(CFLAGS) $(LIBS) -static-libstdc++ -static-libgcc -Wl,--wrap=memcpy


# This way of building dependencies (per-file) described at
# http://make.paulandlesley.org/autodep.html

.c.o:
	@mkdir -p $(DEPDIR)
	$(GCC) -MMD $(CFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d

.cc.o:
	@mkdir -p $(DEPDIR)
	$(GPP) -MMD $(CPPFLAGS) -o $@ -c $<
	@cp $*.d $(df).P; \
	  sed -e 's/#.*//' -e 's/^[^:]*: *//' -e 's/ *\\$$//' \
	  -e '/^$$/ d' -e 's/$$/ :/' < $*.d >> $(df).P; \
	  rm -f $*.d


# include the .P dependency files, but don't warn if they don't exist (the -)
-include $(CPPSRCS:%.cc=$(DEPDIR)/%.P)
-include $(CSRCS:%.c=$(DEPDIR)/%.P)
# The following applies if we don't use a dependency directory:
#-include $(SRCS:.cc=.P)

tags: $(SRCS) *.h
	ctags --language-force=c++ --extra=+q --fields=+i --excmd=n *.c *.cc *.h

clean:
	rm -f $(EXEC) $(OBJS)

clean-deps:
	rm -f $(DEPDIR)/*.P
