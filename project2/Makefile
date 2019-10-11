CC=g++
CFLAGS=-std=c++11 -O2 -larmadillo
TARGETS=rf-antena

all: $(TARGETS)

$(TARGETS): % : %.cpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm -rf *.o *~ $(TARGETS)