CC=gcc
CFLAGS=-I.
DEPS = openptv.h

openptv: openptv.c
	$(CC) -o openptv openptv.c -I. -I/usr/local/lib -loptv