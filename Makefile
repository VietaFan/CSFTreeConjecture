all: matchfinder.exe

matchfinder.exe: matchfinder.o printers.o DegSeqTree.o 
	g++ -Wall -std=c++11 -O3 -o matchfinder3.exe matchfinder3.o DegSeqTree.o printers.o -pthread
	
matchfinder.o: matchfinder.cpp
	g++ -Wall -O3 -c matchfinder.cpp -std=c++0x -pthread

DegSeqTree.o: DegSeqTree.cpp
	g++ -O3 -c DegSeqTree.cpp -std=c++0x

printers.o: printers.cpp
	g++ -O3 -c printers.cpp -std=c++0x

