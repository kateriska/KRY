CXX=g++
#CXX=clang
XLOGINXX=xforto00

OBJ=kry.cpp#[PUT OBJECT FILES HERE]
BIN=kry

CXXFLAGS:=-Wall -Wextra -Wsuggest-override -Wnull-dereference -Wshadow -Wold-style-cast -pedantic -std=c++11
GMPFLAGS:=-L/gmp_install/lib -lgmp

LINK.o = $(LINK.cpp)

all: CXXFLAGS += -Ofast -march=native -flto
all: kry

debug: CXXFLAGS += -g3 -fsanitize=address,undefined -fno-omit-frame-pointer
debug: kry

kry: $(OBJ)
		$(CXX) $(CXXFLAGS) $(LDFLAGS) $(OBJ) -o $(BIN) $(GMPFLAGS)

pack:
		zip $(XLOGINXX).zip *.cpp *.hpp  Makefile doc.pdf

dep:
		g++ *.cpp -MM >> Makefile
