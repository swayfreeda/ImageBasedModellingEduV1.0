OBJS=\
box.o\
box_bld.o\
cdmath3d.o\
coldet.o\
coldet_bld.o\
multiobject.o\
mytritri.o\
sysdep.o\
tritri.o

libcoldet.a: $(OBJS)
	ar cr libcoldet.a $(OBJS)
	
box.o: box.cpp
	g++ -c -DGCC -O2 box.cpp

box_bld.o: box_bld.cpp
	g++ -c -DGCC -O2 box_bld.cpp

cdmath3d.o: cdmath3d.cpp
	g++ -c -DGCC -O2 cdmath3d.cpp

coldet.o: coldet.cpp
	g++ -c -DGCC -O2 coldet.cpp

coldet_bld.o: coldet_bld.cpp
	g++ -c -DGCC -O2 coldet_bld.cpp

multiobject.o: multiobject.cpp
	g++ -c -DGCC -O2 multiobject.cpp

mytritri.o: mytritri.cpp
	g++ -c -DGCC -O2 mytritri.cpp

sysdep.o: sysdep.cpp
	g++ -c -DGCC -O2 sysdep.cpp

tritri.o: tritri.c
	gcc -c -DGCC -O2 tritri.c
	
	
