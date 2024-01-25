#include <stdio.h>
#include <stdlib.h>
#include "stdint.h"
#include "grower.h"
#include "glider.h"
#include "beehive.h"
#include <omp.h>

#define BOARD_SIZE 3000
#define ITERATIONS 5000

uint8_t board[BOARD_SIZE][BOARD_SIZE];
uint8_t newBoard[BOARD_SIZE][BOARD_SIZE];

void initializeBoard() {
    // Clear the board
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            board[i][j] = 0;
        }
    }

    // Define the 'grower' pattern
    int startRow = 1500;
    int startCol = 1500;

    for (int i = 0; i < GROWER_HEIGHT; ++i) {
        for (int j = 0; j < GROWER_WIDTH; ++j) {
            board[startRow + i][startCol + j] = grower[i][j];
        }
    }
}

void printBoard() {
    // Function to print the current state of the board
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            printf("%c", (board[i][j] == 1) ? '*' : ' ');
        }
        printf("\n");
    }
    printf("\n");
}

void updateCell(int i, int j) {
    // Function to update the state of a cell based on Game of Life rules
    int neighbors = 0;

    // Count the number of live neighbors
    for (int x = i - 1; x <= i + 1; ++x) {
        for (int y = j - 1; y <= j + 1; ++y) {
            if (x >= 0 && x < BOARD_SIZE && y >= 0 && y < BOARD_SIZE && !(x == i && y == j)) {
                neighbors += board[x][y];
            }
        }
    }

    // Apply Game of Life rules
    if (board[i][j] == 1) {
        newBoard[i][j] = (neighbors == 2 || neighbors == 3) ? 1 : 0;
    } else {
        newBoard[i][j] = (neighbors == 3) ? 1 : 0;
    }
}

void updateBoard() {
    // Function to update the entire board based on Game of Life rules in parallel
    #pragma omp parallel for
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            updateCell(i, j);
        }
    }
}

int countPopulation() {
    // Function to count the population of the board
    int population = 0;
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            population += board[i][j];
        }
    }
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
    double start_time = omp_get_wtime();

    // Main loop for iterations
    for (int iter = 0; iter < ITERATIONS; ++iter) {
        updateBoard();

        // Swap the boards after each iteration
        #pragma omp parallel for
        for (int i = 0; i < BOARD_SIZE; ++i) {
            for (int j = 0; j < BOARD_SIZE; ++j) {
                board[i][j] = newBoard[i][j];
            }
        }

        if (iter % 100 == 0) {
            // Print the board every 100 iterations
            printf("Iteration %d\n", iter);
            printBoard();
        }

        // Verify patterns after each iteration
        verifyBeehivePopulation(iter);
        verifyGliderPopulation(iter);
        verifyGrowerPopulation(iter);
    }

    double end_time = omp_get_wtime();
    printf("Execution Time: %f seconds\n", end_time - start_time);

    return 0;
}
