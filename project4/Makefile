SHELL := /bin/bash
CMAKE_CPPFLAGS += -Wno-c++11-extensions

all : ivp p_table
ivp:
	g++ $(CMAKE_CPPFLAGS) -o ivp ivp.cpp
p_table: p_table.cpp
	g++ $(CMAKE_CPPFLAGS) -o p_table p_table.cpp
clean:
	rm -rf p_table ivp e*_out.csv m*_out.csv r*out.csv