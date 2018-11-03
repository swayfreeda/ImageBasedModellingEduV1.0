CC=g++
INC=-I ../../src
LFLAGS=-L../../src -lcoldet

unit: main.o ../../src/libcoldet.a
	$(CC) -O2 main.o -o unit $(LFLAGS)

../../src/libcoldet.a:
	make -C ../../src
	
main.o: main.cpp
	$(CC) $(INC) -DGCC -c main.cpp
	
	
