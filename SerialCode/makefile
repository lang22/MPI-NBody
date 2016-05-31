#make file for my serial version of the project
CC=gcc
CFLAGS=-pg -O2 -c -Wall
LDFLAGS=-lm -lGL -lGLU -lglfw3 -lX11 -lXxf86vm -lpthread -ldl -lXrandr -lXi -lXinerama -lXcursor
SOURCES=nbfast.c
OBJECTS=nbfast.o
EXECUTABLE=NBODY.x

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) -pg $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm *.o $(EXECUTABLE); rm res/*
