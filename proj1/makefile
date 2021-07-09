all: Driver

Driver: Driver.o
	g++ $(CPPFLAGS) -o Driver Driver.o


Driver.o: Driver.cpp image.cpp image.h ReadImage.cpp ReadImageHeader.cpp WriteImage.cpp

clean:
	rm -f Driver; rm *.o
