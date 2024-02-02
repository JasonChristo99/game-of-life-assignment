#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>

// Constants defining the size of the board and maximum iterations
#define BOARD_SIZE 3000
#define MAX_ITERATIONS 5000

// Constants defining the starting positions of various patterns on the board
#define GROWER_START_ROW 1500
#define GROWER_START_COL 1500
#define GLIDER_START_ROW 1500
#define GLIDER_START_COL 1500
#define BEEHIVE_START_ROW 1500
#define BEEHIVE_START_COL 1500

// Population expectations for verification at specific iterations
#define GROWER_POPULATION_GEN10 49
#define GROWER_POPULATION_GEN100 138
#define GLIDER_POPULATION 5
#define BEEHIVE_POPULATION 6

// Debug printing controls
#define DEBUG_PRINT_RANK_INFO 0
#define DEBUG_PRINT_REGION 0
#define DEBUG_RUN_VERIFICATION_CHECKS 0
#define DEBUG_PRINT_TOTAL_POPULATION 1
#define DEBUG_PRINT_EXECUTION_TIME 1

// Include patterns
#include "grower.h"
#include "glider.h"
#include "beehive.h"

// Macro to convert 2D indices to a 1D index for a flattened array
#define INDEX(row, col) ((row) * BOARD_SIZE + (col))

// Function to convert local coordinates to global coordinates
void local_coords_to_global(int *global_row, int *global_col, int row, int col, int startRow, int startCol) {
    *global_row = row + startRow;
    *global_col = col + startCol;
}

// Function to print a region of the board
void print_region(uint8_t *board, int startRow, int startCol, int height, int width) {
    for (int i = startRow; i < startRow + height; ++i) {
        for (int j = startCol; j < startCol + width; ++j) {
            printf("%c", board[INDEX(i, j)] ? '#' : '^');
        }
        printf("\n");
    }
}

// Function to initialize the local board with a specific pattern
void initialize_local_board(uint8_t *board, int block_start, int block_size) {
    // Fill the local board with the grower pattern
    // Only the part of the pattern that is within the local board boundaries will be copied
    for (int i = 0; i < block_size + 2; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            if (i == 0 || i == block_size + 1) {
                // Set boundary rows to 0
                board[INDEX(i, j)] = 0;
                continue;
            }
            int row_global, col_global;
            local_coords_to_global(&row_global, &col_global, i, j, block_start, 0);
            if (row_global >= GROWER_START_ROW && row_global < GROWER_START_ROW + GROWER_HEIGHT &&
                col_global >= GROWER_START_COL && col_global < GROWER_START_COL + GROWER_WIDTH) {
                // Copy the grower pattern to the local board
                board[INDEX(i, j)] = grower[row_global - GROWER_START_ROW][col_global - GROWER_START_COL];
            } else
                // Set cells outside the grower pattern to 0
                board[INDEX(i, j)] = 0;
        }
    }
}

// Function to calculate the population of the local board
int local_board_population(uint8_t *board, int block_size) {
    int population = 0;
    for (int i = 1; i <= block_size; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            population += board[INDEX(i, j)];
        }
    }
    return population;
}

// Function to apply the Game of Life rules to a cell
uint8_t apply_rules(int current_state, int neighbors) {
    if (current_state == 1) {
        // Cell is alive
        if (neighbors < 2 || neighbors > 3) {
            // Rule 1 and 3: Cell dies
            return 0;
        } else {
            // Rule 2: Cell survives
            return 1;
        }
    } else {
        // Cell is dead
        if (neighbors == 3) {
            // Rule 4: Cell becomes alive
            return 1;
        } else {
            // Cell remains dead
            return 0;
        }
    }
}

// Function to run one iteration of the Game of Life for the local region of the board
void run_game_of_life(uint8_t *current, uint8_t *next, int block_size) {
#pragma omp parallel for collapse(2)
    for (int i = 1; i <= block_size; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            // Count neighbors for each cell
            int neighbors = 0;
            for (int ni = -1; ni <= 1; ++ni) {
                for (int nj = -1; nj <= 1; ++nj) {
                    if (ni == 0 && nj == 0) continue;

                    int ni_ = i + ni;

                    // Check if the neighbor is within the local board boundaries
                    int nj_ = j + nj;

                    // Check if the neighbor is within the board boundaries
                    if (nj_ >= 0 && nj_ < BOARD_SIZE) {
                        neighbors += current[INDEX(ni_, nj_)];
                    }
                }
            }

            // Apply Game of Life rules
            next[INDEX(i, j)] = apply_rules(current[INDEX(i, j)], neighbors);
        }
    }
}

// Function to calculate the population of the entire board across all processes
int total_board_population(uint8_t *local_board, int block_size) {
    int local_population = local_board_population(local_board, block_size);

    // Use MPI_Allreduce to perform a reduction operation across all processes
    int total_population;
    if (MPI_Allreduce(&local_population, &total_population, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error in MPI_Allreduce operation.\n");
        exit(EXIT_FAILURE);
    }

    return total_population;
}

// Function to perform verification checks at specific iterations
void run_verification(uint8_t *board, int iter, int block_size) {
    // Check total population for specific iterations
    if (iter == 10) {
        int total_population = total_board_population(local_board, block_size);
        printf("Generation 10: Grower pattern population = %d\n", total_population);
        if (total_population != GROWER_POPULATION_GEN10) {
            printf("ERROR: Incorrect population for generation 10\n");
            exit(1);
        }
    }

    if (iter == 100) {
        int total_population = total_board_population(local_board, block_size);
        printf("Generation 100: Grower pattern population = %d\n", total_population);
        if (total_population != GROWER_POPULATION_GEN100) {
            printf("ERROR: Incorrect population for generation 100\n");
            exit(1);
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, size;

    // Start up MPI
    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        fprintf(stderr, "Error initializing MPI.\n");
        exit(EXIT_FAILURE);
    }
    if (MPI_Comm_rank(MPI_COMM_WORLD, &rank) != MPI_SUCCESS) {
        fprintf(stderr, "Error in MPI_Comm_rank.\n");
        exit(EXIT_FAILURE);
    }
    if (MPI_Comm_size(MPI_COMM_WORLD, &size) != MPI_SUCCESS) {
        fprintf(stderr, "Error in MPI_Comm_size.\n");
        exit(EXIT_FAILURE);
    }

    const bool am_master = 0 == rank;

    // Divide the board into blocks for parallel processing (horizontal segments)
    int block_size = BOARD_SIZE / size;
    int block_start = rank * block_size;
    int block_end = block_start + block_size;
    if (DEBUG_PRINT_RANK_INFO) {
        // Print the block start and end indices for each process
        printf("Rank %d: Block start = %d, Block end = %d\n", rank, block_start, block_end);
    }

    // Allocate memory for the local and next generation boards
    uint8_t *local_board = (uint8_t *) malloc(BOARD_SIZE * (block_size + 2) * sizeof(uint8_t));
    uint8_t *next_gen_board = (uint8_t *) malloc(BOARD_SIZE * (block_size + 2) * sizeof(uint8_t));

    if (local_board == NULL || next_gen_board == NULL) {
        fprintf(stderr, "Error allocating memory.\n");
        exit(EXIT_FAILURE);
    }

    // Initialize the local board with the grower pattern
    initialize_local_board(local_board, block_start, block_size);

    // Record the start time for measuring execution time
    double start_time = MPI_Wtime();

    for (int iter = 0; iter <= MAX_ITERATIONS; ++iter) {
        if (DEBUG_PRINT_RANK_INFO) {
            // Print the population of the local board for each iteration
            printf("Rank %d: Iteration %d, Population = %d\n", rank, iter,
                   local_board_population(local_board, block_size));
        }

        if (DEBUG_PRINT_REGION) {
            print_region(local_board, GROWER_START_ROW, GROWER_START_COL, GROWER_HEIGHT, GROWER_WIDTH);
        }

        if (DEBUG_RUN_VERIFICATION_CHECKS) {
            // Verification checks for specific iterations
            run_verification(local_board, iter, block_size);
        }

        if (DEBUG_PRINT_TOTAL_POPULATION && am_master) {
            // Calculate and print the total population of the entire board
            int total_population = total_board_population(local_board, block_size);
            printf("Total population after iteration %d: %d\n", iter, total_population);
        }

        // Communication between neighboring processes using MPI
        // Identify left and right neighbors, set MPI_PROC_NULL if at the boundary
        int left_neighbour = (rank == 0) ? MPI_PROC_NULL : (rank - 1);
        int right_neighbour = (rank == size - 1) ? MPI_PROC_NULL : (rank + 1);

        // Define MPI request objects for non-blocking communication
        MPI_Request left_send_request, right_send_request;
        MPI_Request left_recv_request, right_recv_request;

        // Initiate non-blocking send operations to left and right neighbors
        if (block_start > 0) {
            // Communication with the left neighbor
            if (MPI_Isend(&local_board[INDEX(1, 0)], BOARD_SIZE, MPI_UINT8_T, left_neighbour, 0, MPI_COMM_WORLD,
                          &left_send_request) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Isend to left neighbor.\n");
                exit(EXIT_FAILURE);
            }
            if (MPI_Irecv(&local_board[INDEX(0, 0)], BOARD_SIZE, MPI_UINT8_T, left_neighbour, 0, MPI_COMM_WORLD,
                          &left_recv_request) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Irecv from left neighbor.\n");
                exit(EXIT_FAILURE);
            }
        }

        if (block_end < BOARD_SIZE) {
            // Communication with the right neighbor
            if (MPI_Isend(&local_board[INDEX(block_size, 0)], BOARD_SIZE, MPI_UINT8_T, right_neighbour, 0,
                          MPI_COMM_WORLD,
                          &right_send_request) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Isend to right neighbor.\n");
                exit(EXIT_FAILURE);
            }
            if (MPI_Irecv(&local_board[INDEX(block_size + 1, 0)], BOARD_SIZE, MPI_UINT8_T, right_neighbour, 0,
                          MPI_COMM_WORLD, &right_recv_request) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Irecv from right neighbor.\n");
                exit(EXIT_FAILURE);
            }
        }

        // Wait for completion of all non-blocking communication operations
        if (block_start > 0) {
            if (MPI_Wait(&left_send_request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Wait for left_send_request.\n");
                exit(EXIT_FAILURE);
            }
            if (MPI_Wait(&left_recv_request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Wait for left_recv_request.\n");
                exit(EXIT_FAILURE);
            }
        }
        if (block_end < BOARD_SIZE) {
            if (MPI_Wait(&right_send_request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Wait for right_send_request.\n");
                exit(EXIT_FAILURE);
            }
            if (MPI_Wait(&right_recv_request, MPI_STATUS_IGNORE) != MPI_SUCCESS) {
                fprintf(stderr, "Error in MPI_Wait for right_recv_request.\n");
                exit(EXIT_FAILURE);
            }
        }

        // Run one iteration of the Game of Life for the local region
        run_game_of_life(local_board, next_gen_board, block_size);

        // Swap pointers to update the local board for the next iteration
        uint8_t *temp = local_board;
        local_board = next_gen_board;
        next_gen_board = temp;
    }

    // Record the end time for measuring execution time
    double end_time = MPI_Wtime();

    // Print total execution time
    if (DEBUG_PRINT_EXECUTION_TIME && am_master) {
        printf("Total execution time: %f seconds\n", end_time - start_time);
    }

    // Print the final time : the time taken by the slowest process to complete the execution
    double max_time;
    if (MPI_Reduce(&end_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD) != MPI_SUCCESS) {
        fprintf(stderr, "Error in MPI_Reduce operation.\n");
        exit(EXIT_FAILURE);
    }
    printf("Max execution time: %f seconds\n", max_time - start_time);


    // Free allocated memory
    free(local_board);
    free(next_gen_board);

    // Finalize MPI
    if (MPI_Finalize() != MPI_SUCCESS) {
        fprintf(stderr, "Error finalizing MPI.\n");
        exit(EXIT_FAILURE);
    }

    return 0;
}
