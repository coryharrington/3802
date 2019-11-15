SHELL := /bin/bash
C_FLAGS += -Wall -g -pthread
TARGETS=test_server
all : $(TARGETS)

$(TARGETS): %: %.c
	$(CC) $(CFLAGS) -o $@ $<

clean: 
	@rm -f $(TARGETS) *.o a.out 