all:
	clang++ main.cpp -o main -Wall -std=c++11 -O3 -DWUP_NO_OPENCV -DWUP_NO_MPICH -lpthread -I ../wup/cpp/include

run:
	time ./main -maxGen 50 -popSize 5

debug:
	clang++ main.cpp -o main -Wall -std=c++11 -g -DWUP_NO_OPENCV -DWUP_NO_MPICH -lpthread -I ../wup/cpp/include
	gdb main

val:
	clang++ main.cpp -o main -Wall -std=c++11 -O1 -g -DWUP_NO_OPENCV -DWUP_NO_MPICH -lpthread -I ../wup/cpp/include
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes --verbose ./main	-maxFitEval 100000 -maxGen 100 -popSize 50 -dims 10 -ntupleDims 20 -maxRun 1 -threads 1
	
