
CC = gcc
CFLAGS = -O3

reduce: reduce.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o reduce reduce.c snapshot_io.c data_ops.c -lm -O3

convert2ifirt: convert2ifrit.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o convert2ifrit convert2ifrit.c snapshot_io.c data_ops.c -lm -O3

getBBox: getBBox.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o getBBox getBBox.c snapshot_io.c data_ops.c -lm -O3

santaBarbara: santaBarbara.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o santaBarbara santaBarbara.c snapshot_io.c data_ops.c -lm -O3

all: reduce convert2ifirt getBBox santaBarbara

clean:
	rm -rf reduce convert2ifrit getBBox
