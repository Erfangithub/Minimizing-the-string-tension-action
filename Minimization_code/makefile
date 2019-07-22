
All: main

main: main.o set_init.o maction1ij.o maction1.o maction2ij.o fresnel.o maction2.o mfintz2sin.o mfintzsin.o fsinint.o fcosint.o partderij2.o
	g++ main.o set_init.o maction1ij.o maction1.o fresnel.o maction2ij.o maction2.o mfintz2sin.o mfintzsin.o fsinint.o fcosint.o partderij2.o -o main -lgsl -lgslcblas -lgmp -lmpfr

main.o: main.cpp set_init.h maction1ij.h maction1.h fresnel.h maction2ij.h maction2.h mfintz2sin.h mfintzsin.h fsinint.h fcosint.h partderij2.h
	g++ -c main.cpp

maction1ij.o: maction1ij.cpp maction1ij.h
	g++ -c maction1ij.cpp

set_init.o: set_init.cpp set_init.h
	g++ -c set_init.cpp

maction1.o: maction1.cpp maction1.h
	g++ -c maction1.cpp

maction2ij.o: maction2ij.cpp maction2ij.h fresnel.h
	g++ -c maction2ij.cpp

fresnel.o: fresnel.cpp fresnel.h
	g++ -c fresnel.cpp

maction2.o: maction2.cpp maction2.h
	g++ -c maction2.cpp

mfintz2sin.o: mfintz2sin.cpp mfintz2sin.h fresnel.h
	g++ -c mfintz2sin.cpp

mfintzsin.o: mfintzsin.cpp mfintzsin.h fresnel.h
	g++ -c mfintzsin.cpp

fsinint.o: fsinint.cpp fsinint.h fresnel.h
	g++ -c fsinint.cpp

fcosint.o: fcosint.cpp fcosint.h fresnel.h
	g++ -c fcosint.cpp

partderij2.o: partderij2.cpp partderij2.h fsinint.h fresnel.h mfintzsin.h mfintz2sin.h
	g++ -c partderij2.cpp

clean:
	rm -f *o main


# main.cpp maction1ij.cpp set_init.cpp maction1.cpp maction2ij.cpp fresnel.cpp maction2.cpp -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lgmp -lmpfr