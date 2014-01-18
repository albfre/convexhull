CC = g++
CFLAGS = -g -O2 -Wall -pedantic
FILE = ConvexHull
SRCS = $(FILE).cpp $(FILE).h
OBJS = $(FILE).o
TEST = Test$(FILE)
TESTSRC = $(TEST).cpp
PROFILECFLAGS = -p -O2 -Wall -pedantic

all: $(SRCS) \
; $(CC) $(CFLAGS) -c $(SRCS)

profile: $(SRCS) \
; $(CC) $(PROFILECFLAGS) -c $(SRCS)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

testprofile: $(TESTSRC) \
; $(CC) $(PROFILECFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS)

write: \
; echo $(OBJS)
