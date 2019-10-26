SHELL := /bin/bash
CMAKE_CPPFLAGS += -Wno-c++11-extensions

all : sommerfeld
sommerfeld: sommerfeld.cpp
	g++ $(CMAKE_CPPFLAGS) -o sommerfeld sommerfeld.cpp
clean:
	rm -rf sommerfeld
