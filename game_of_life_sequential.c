#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include "grower.h"
#include "glider.h"
#include "beehive.h"
#include <time.h>

#define BOARD_SIZE 3000
#define ITERATIONS 1

// Define coordinates for the region to print
#define PRINT_START_X 1495
#define PRINT_START_Y 1580
#define PRINT_END_X 1525
#define PRINT_END_Y 1620

uint8_t board[BOARD_SIZE][BOARD_SIZE];
uint8_t newBoard[BOARD_SIZE][BOARD_SIZE];

void initializeBoard() {
    // Clear the board
    for (int i = 0; i < BOARD_SIZE; i++) {
        for (int j = 0; j < BOARD_SIZE; j++) {
            board[i][j] = 0;
        }
    }

    // Define the 'grower' pattern
    int startRow = 1500;
    int startCol = 1500;

    for (int i = 0; i < GROWER_HEIGHT; i++) {
        for (int j = 0; j < GROWER_WIDTH; j++) {
            board[startRow + i][startCol + j] = grower[i][j];
        }
    }
}

void printBoard(int startX, int startY, int endX, int endY) {
    // Function to print a specific region of the board
    for (int i = startX; i < endX; i++) {
        for (int j = startY; j < endY; j++) {
            printf("%c", (board[i][j] == 1) ? '*' : ' ');
        }
        printf("\n");
    }
    printf("\n");
}

void updateCell(int i, int j) {
    // Function to update the state of a cell based on the specified rules

    int liveNeighbors = 0;

    // Count the number of live neighbors
    for (int x = i - 1; x <= i + 1; x++) {
        for (int y = j - 1; y <= j + 1; y++) {
            if (x >= 0 && x < BOARD_SIZE && y >= 0 && y < BOARD_SIZE && !(x == i && y == j)) {
                liveNeighbors += board[x][y];
            }
        }
    }

    // Apply the specified rules
    if (board[i][j] == 1) {
        if (liveNeighbors < 2 || liveNeighbors > 3) {
            // Rule: Any live cell with fewer than two live neighbors or more than three live neighbors dies
            newBoard[i][j] = 0;
        } else {
            // Rule: Any live cell with two or three live neighbors lives, unchanged, to the next generation
            newBoard[i][j] = 1;
        }
    } else {
        if (liveNeighbors == 3) {
            // Rule: Any dead cell with exactly three live neighbors will come to life
            newBoard[i][j] = 1;
        } else {
            newBoard[i][j] = 0;
        }
    }
}


void updateBoard() {
    // Function to update the entire board based on Game of Life rules
    for (int i = 0; i < BOARD_SIZE; i++) {
        for (int j = 0; j < BOARD_SIZE; j++) {
            updateCell(i, j);
        }
    }

    // Swap the boards after each iteration
    for (int i = 0; i < BOARD_SIZE; i++) {
        for (int j = 0; j < BOARD_SIZE; j++) {
            board[i][j] = newBoard[i][j];
        }
    }
}

int countPopulation() {
    // Function to count the population of the board
    int population = 0;
    for (int i = 0; i < BOARD_SIZE; i++) {
        for (int j = 0; j < BOARD_SIZE; j++) {
            population += board[i][j];
        }
    }
    // printf("Population: %d\n", population);
    return population;
}


void verifyBeehivePopulation(int generation) {
    // Verify that the 'beehive' pattern has a population of exactly 6 in all generations
    if (countPopulation() != 6) {
        printf("Verification failed for Beehive pattern at generation %d\n", generation);
        exit(1);
    }
}

void verifyGliderPopulation(int generation) {
    // Verify that the 'glider' pattern has a population of exactly 5 in all generations until the pattern leaves the board
    if (generation <= ITERATIONS) {
        if (countPopulation() != 5) {
            printf("Verification failed for Glider pattern at generation %d\n", generation);
            exit(1);
        }
    }
}

void verifyGrowerPopulation(int generation) {
    // Verify that the 'grower' pattern has a population of 49 at generation 10, and a population of 138 at generation 100
    int expectedPopulation = (generation == 10) ? 49 : (generation == 100) ? 138 : -1;

    if (expectedPopulation != -1 && countPopulation() != expectedPopulation) {
        printf("Verification failed for Grower pattern at generation %d\n", generation);
        exit(1);
    }
}

int main() {
    initializeBoard();
    clock_t start_time = clock();

    // Main loop for iterations
    for (int iter = 0; iter <= ITERATIONS; iter++) {
        

        if (iter % 1 == 0) {
            // Print the board every 100 iterations
            printf("Iteration %d\n", iter);
            printBoard(PRINT_START_X, PRINT_START_Y, PRINT_END_X, PRINT_END_Y);
        }
        
        printf("Population: %d\n", countPopulation());

        // Verify patterns after each iteration
        // verifyBeehivePopulation(iter);
        // verifyGliderPopulation(iter);
        verifyGrowerPopulation(iter);
        
        updateBoard();
    }

    clock_t end_time = clock();
    double execution_time = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;
    printf("Execution Time: %f seconds\n", execution_time);

    return 0;
}
