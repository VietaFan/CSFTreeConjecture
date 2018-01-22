all: matchfinder.exe

matchfinder.exe: matchfinder.o printers.o DegSeqTree.o 
	g++ -Wall -std=c++11 -O3 -o matchfinder.exe matchfinder.o DegSeqTree.o printers.o -pthread
	
matchfinder.o: matchfinder.cpp
	g++ -Wall -O3 -c matchfinder.cpp -std=c++0x -pthread

DegSeqTree.o: DegSeqTree.cpp
	g++ -O3 -c DegSeqTree.cpp -std=c++0x

printers.o: printers.cpp
	g++ -O3 -c printers.cpp -std=c++0x
	
matchfinder2.exe: matchfinder2.o TreeGenerator.o Tree.o
	g++ -Wall -std=c++11 -O3 -o matchfinder2.exe matchfinder2.o TreeGenerator.o Tree.o

matchfinder2.o: matchfinder2.cpp
	g++ -Wall -O3 -c matchfinder2.cpp -std=c++11 
	
TreeGenerator.o: TreeGenerator.cpp
	g++ -Wall -O3 -c TreeGenerator.cpp -std=c++11 

Tree.o: Tree.cpp
	g++ -Wall -O3 -c Tree.cpp -std=c++11
