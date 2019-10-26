SHELL := /bin/bash
CMAKE_CPPFLAGS += -Wno-c++11-extensions

all : sommerfeld
bashtesting:
	./test CMAKE_CPPFLAGS
rf-mode: rf-mode.cpp
	g++ $(CMAKE_CPPFLAGS) -o rf-mode rf-mode.cpp
clean:
	rm -rf .bessel rf-mode