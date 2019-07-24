cc = g++
prom = Dancing_link
CXXFLAGS=-g -std=c++11 -Wall -pedantic -Wextra -lstdc++ -g
deps = $(shell find ./ -name "*.h")
src = $(shell find ./ -name "*.cpp")
obj = $(src:%.cpp=%.o) 

$(prom): $(obj)
	$(cc) -std=c11 -o  $(prom) $(obj) 

%.o: %.c $(deps)
	$(cc) -std=c11 -c $< -o  $@

clean:
	rm -rf $(obj) $(prom)