all: Correlate

Correlate: Correlate.o
	g++ $(CPPFLAGS) -o Correlate Correlate.o


Correlate.o: Correlate.cpp image.cpp image.h ReadImage.cpp ReadImageHeader.cpp WriteImage.cpp

clean:
	rm -f Correlate; rm *.o
