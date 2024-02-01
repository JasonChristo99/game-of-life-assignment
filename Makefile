# Compiler
CC = mpicc

# Compiler flags
CFLAGS = -Wall -Wextra -std=c11 -O3 -fopenmp -g

# Executable name
EXECUTABLE = game_of_life

# Source files
SRC = game_of_life.c

# Object files
OBJ = $(SRC:.c=.o)

# Specify existing header files
HEADERS = beehive.h glider.h grower.h

all: $(EXECUTABLE) $(HEADERS)

$(EXECUTABLE): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) -o $(EXECUTABLE)

# Rule to compile source files
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJ)


