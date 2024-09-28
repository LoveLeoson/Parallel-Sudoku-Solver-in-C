CFLAGS=-W -O3 -g -march=native
INCLUDES=-I/opt/X11/include
LDFLAGS=-L/opt/X11/lib -lX11 -lm -ffunction-sections -fopenmp 
FILE=soduku.c

soduku: soduku.c
	gcc -o soduku soduku.c $(CFLAGS) $(LDFLAGS)
clean:
	rm -f ./soduku *.o