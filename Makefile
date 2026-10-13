# pi4decoder - plain C99, libm only.
CC      ?= gcc
CFLAGS  ?= -O3 -std=c99 -Wall -Wextra -pedantic
LDLIBS   = -lm

LIB_SRC  = src/pi4_spec.c src/pi4_fano.c src/pi4_fft.c src/pi4_decode.c
LIB_OBJ  = $(LIB_SRC:.c=.o)

all: pi4decoder pi4gen libpi4.a

libpi4.a: $(LIB_OBJ)
	$(AR) rcs $@ $^

pi4decoder: src/main.o src/wav.o libpi4.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

pi4gen: src/pi4gen.o src/wav.o libpi4.a
	$(CC) $(CFLAGS) -o $@ $^ $(LDLIBS)

tests/test_spec: tests/test_spec.c libpi4.a
	$(CC) $(CFLAGS) -Isrc -o $@ $^ $(LDLIBS)

tests/test_decode: tests/test_decode.c libpi4.a
	$(CC) $(CFLAGS) -Isrc -o $@ $^ $(LDLIBS)

tests/test_noise: tests/test_noise.c libpi4.a
	$(CC) $(CFLAGS) -Isrc -o $@ $^ $(LDLIBS)

test: tests/test_spec tests/test_decode
	./tests/test_spec
	./tests/test_decode

# small static build for size measurement
small: $(LIB_SRC) src/main.c src/wav.c
	$(CC) -Os -std=c99 -s -ffunction-sections -fdata-sections -Wl,--gc-sections -o pi4decoder-small $^ $(LDLIBS)

clean:
	rm -f src/*.o libpi4.a pi4decoder pi4gen pi4decoder-small tests/test_spec tests/test_decode tests/test_noise

.PHONY: all test clean small
