CC=g++
VTKLIBDIR=$(shell dirname $$( ldd $$( which vtk ) | grep vtkIO.so | sed -r 's/.* => (.*) \(.*\)/\1/g' ))
VTKVER=$(shell basename $$( dirname $$( ldd $$( which vtk ) | grep vtkIO.so | sed -r 's/.* => (.*) \(.*\)/\1/g' )) | cut -f2 -d-)
CXXFLAGS=-Wno-deprecated -std=c++11 -I/usr/include/vtk-$(VTKVER) -I/usr/include/gmsh
LDFLAGS=-L$(VTKLIBDIR) -lgmsh -lvtkCommon -lvtkIO -lvtkFiltering

default: all

clean:
	rm -f mesher

all: mesher

mesher: mesher.cpp
	echo $(VTKLIBDIR)
	$(CC) $(CXXFLAGS) $(LDFLAGS) mesher.cpp -o mesher

.PHONY: all default clean
