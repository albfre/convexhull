CC = g++
CFLAGS = -g -funroll-loops -std=gnu++0x -O2 -Wall -pedantic
FILE = ConvexHull
SRCS = $(FILE).cpp $(FILE).h
OBJS = $(FILE).o
TEST = Test$(FILE)
TESTSRC = $(TEST).cpp
PROFILECFLAGS = -pg -funroll-loops -std=gnu++0x -O2 -Wall -pedantic

all: product test

product: $(SRCS) \
; $(CC) $(CFLAGS) -c $(SRCS)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

profile: productprofile testprofile

productprofile: $(SRCS) \
; $(CC) $(PROFILECFLAGS) -c $(SRCS)

testprofile: $(TESTSRC) \
; $(CC) $(PROFILECFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS)

write: \
; echo $(OBJS)
