CC = gcc

default: libicsi2022.a

libicsi2022.a: icsi2022.o
	ar rcs $@ $^

icsi2022.o: icsi2022.c icsi2022.h
	$(CC) -c $< -fPIC

clean:
	rm *.o *.a
