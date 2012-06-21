
CC = gcc
CFLAGS = -O3 -Wno-unused-result -Wno-format

INCLUDES = snapshot_io.h data_ops.h
SOURCES =  snapshot_io.c data_ops.c
OBJECTS = $(SOURCES:.c=.o)
LIBS = -lm

EXECUTABLES = cutout info reduce getBBox cutoutDm

all: $(EXECUTABLES)

$(EXECUTABLES): $(OBJECTS)
	$(CC) $@.c $(OBJECTS) -o $@ $(LIBS)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@ $(LIBS)


.PHONY:clean
clean:
	rm -rf *.o $(EXECUTABLES)
