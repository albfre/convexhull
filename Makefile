CC = g++
CFLAGS = -g -Wall -pedantic
FILE = ConvexHull
SRCS = $(FILE).cpp $(FILE).h
OBJS = $(FILE).o
TEST = Test$(FILE)
TESTSRC = $(TEST).cpp
PROFILECFLAGS = -pg -Wall -pedantic

all: $(SRCS) \
; $(CC) $(CFLAGS) -c $(SRCS)

profile: $(SRCS) \
; $(CC) $(PROFILECFLAGS) -c $(SRCS)

test: $(TESTSRC) \
; $(CC) $(CFLAGS) -o $(TEST) $(TESTSRC) $(OBJS)

clean: \
; $(RM) $(OBJS)

write: \
; echo $(OBJS)
