SHELL := /bin/bash
CMAKE_CPPFLAGS += -Wno-c++11-extensions

all : sommerfeld
sommerfeld: sommerfeld.cpp
	g++ $(CMAKE_CPPFLAGS) -o summerfeld summerfeld.cpp -O 3
clean:
	rm -rf sommerfeld
