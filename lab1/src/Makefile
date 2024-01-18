CC=gcc

MAIN = main.c

CFLAGS = -Wall -Wextra -Werror -std=c11
LDFLAGS = $(shell pkg-config --cflags)
OMPFLAGS = -fopenmp

all: linear parallel clean

linear: $(MAIN)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o linear.out
	./linear.out

parallel: $(MAIN)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(OMPFLAGS) -o parallel.out

clean:
	@rm -rf *.o *.out

codestyle:
	clang-format -i --style=google $(SRC_MAIN_DIR)/*.c
	clang-format -i --style=google $(SRC_MAIN_DIR)/*.h