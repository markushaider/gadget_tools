
CC = gcc
CFLAGS = -O3 -Wno-unused-result

reduce: reduce.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o reduce reduce.c snapshot_io.c data_ops.c -lm -O3

convert2ifirt: convert2ifrit.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o convert2ifrit convert2ifrit.c snapshot_io.c data_ops.c -lm -O3

getBBox: getBBox.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o getBBox getBBox.c snapshot_io.c data_ops.c -lm -O3

cutout: cutout.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o cutout cutout.c snapshot_io.c data_ops.c -lm -O3

santaBarbara: santaBarbara.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o santaBarbara santaBarbara.c snapshot_io.c data_ops.c -lm -O3

info: info.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o info info.c snapshot_io.c data_ops.c -lm -O3

santaBarbara50: santaBarbara50.c snapshot_io.c data_ops.c snapshot_io.h data_ops.h
	$(CC) $(CFLAGS) -o santaBarbara50 santaBarbara50.c snapshot_io.c data_ops.c -lm -O3

all: reduce convert2ifirt getBBox santaBarbara info cutout

clean:
	rm -rf reduce convert2ifrit getBBox info santaBarbara cutout
