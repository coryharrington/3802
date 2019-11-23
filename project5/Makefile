TARGETS=ladder
OBJS=ladder.o matrix.o madd.o mmul.o
SRCS=$(patsubst %,%.c,$(TARGETS))
HEADERS=matrix.h
CC=gcc
CFLAGS= -Wall -pthread -g

all : $(TARGETS)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS)
	
clean: 
	@rm -f $(OBJS) $(TARGET)

testing: matrix.o madd.o mmul.o
	gcc -Wall -pthread -g ladder.c -o ladder matrix.o -lm

%.o : %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $<
