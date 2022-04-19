buildsa:	buildsa.cpp
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib buildsa.cpp -o buildsa -lsdsl -ldivsufsort -ldivsufsort64

querysa:	querysa.cpp
	g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib querysa.cpp -o querysa -lsdsl -ldivsufsort -ldivsufsort64 