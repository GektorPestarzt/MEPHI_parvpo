CC=gcc

MAIN = main.c

CFLAGS = -Wall -Wextra -Werror -std=c11
LDFLAGS = $(shell pkg-config --cflags)
OMPFLAGS = -fopenmp

LOG_DIR = ../log_files

LOG_LINEAR = $(LOG_DIR)/log_linear.txt
LOG_PARALLEL = $(LOG_DIR)/log_parallel.txt

all: linear parallel clean

linear: $(MAIN)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) -o linear.out
	@mkdir -p $(LOG_DIR)
	./linear.out $(LOG_LINEAR) 

parallel: $(MAIN)
	$(CC) $(CFLAGS) $^ $(LDFLAGS) $(OMPFLAGS) -o parallel.out
	@mkdir -p $(LOG_DIR)
	./parallel.out $(LOG_PARALLEL)

clean:
	@rm -rf *.o *.out

codestyle:
	clang-format -i --style=google $(SRC_MAIN_DIR)/*.c
	clang-format -i --style=google $(SRC_MAIN_DIR)/*.h