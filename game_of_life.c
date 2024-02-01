#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>
#include <stdbool.h>

// Constants defining the size of the board and maximum iterations
#define BOARD_SIZE 3000
#define MAX_ITERATIONS 110

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

// Number of threads for OpenMP
#define NUM_THREADS 16

// Include patterns
#include "grower.h"
#include "glider.h"
#include "beehive.h"

// Macro to convert 2D indices to a 1D index for a flattened array
#define INDEX(row, col) ((row) * BOARD_SIZE + (col))

int local_board_index_to_global(int row, int col, int startRow, int startCol) {
    return INDEX(row + startRow, col + startCol);
}

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

// Function to print the entire board
void print_board(uint8_t *board) {
    print_region(board, 0, 0, BOARD_SIZE, BOARD_SIZE);
}

// Function to write the board to a file
void write_board(uint8_t *board, char *filename, int iteration) {
    FILE *fp;
    if (iteration == 0)
        fp = fopen(filename, "w");
    else
        fp = fopen(filename, "a");

    fprintf(fp, "--------Iteration %d--------\n", iteration);
    for (int i = 0; i < BOARD_SIZE; ++i) {
        for (int j = 0; j < BOARD_SIZE - 1; ++j) {
            fprintf(fp, "%d,", board[INDEX(i, j)]);
        }
        fprintf(fp, "%d\n", board[INDEX(i, BOARD_SIZE - 1)]);
    }
    fclose(fp);
}

// Function to initialize the board with a specific pattern
//void initialize_board(uint8_t *board) {
//    for (int i = 0; i < GROWER_HEIGHT; ++i) {
//        for (int j = 0; j < GROWER_WIDTH; ++j) {
//            board[INDEX(GROWER_START_ROW + i, GROWER_START_COL + j)] = grower[i][j];
//        }
//    }
//}

// Function to initialize the local board with a specific pattern
void initialize_local_board(uint8_t *board, int block_start, int block_size) {
    // Fill the local board with the grower pattern
    // The grower pattern is larger than the local board, so only a part of it will be copied
    // Use the GROWER_START_ROW and GROWER_START_COL constants to determine the starting indices
    // Only copy the part of the pattern that is within the local board boundaries (applicable)
    for (int i = 0; i < block_size + 2; ++i) {
        for (int j = 0; j < BOARD_SIZE; ++j) {
            if (i == 0 || i == block_size + 1) {
                board[INDEX(i, j)] = 0;
                continue;
            }
            int row_global, col_global;
            local_coords_to_global(&row_global, &col_global, i, j, block_start, 0);
            if (row_global >= GROWER_START_ROW && row_global < GROWER_START_ROW + GROWER_HEIGHT &&
                col_global >= GROWER_START_COL && col_global < GROWER_START_COL + GROWER_WIDTH) {
//                printf("Initializing board[%d][%d] = grower[%d][%d]\n", i, j, row_global - GROWER_START_ROW, col_global - GROWER_START_COL);
                board[INDEX(i, j)] = grower[row_global - GROWER_START_ROW][col_global - GROWER_START_COL];
            } else
                board[INDEX(i, j)] = 0;
        }
    }


}

// Function to calculate the population of the entire board
int board_population(uint8_t *board, int block_size) {
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
//#pragma omp parallel for collapse(2) num_threads(NUM_THREADS)
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

// Function to perform verification checks at specific iterations
//void run_verification(uint8_t *board, int iter) {
//    if (iter == 10) {
//        int population = board_population(board);
//        printf("Generation 10: Grower pattern population = %d\n", population);
//        if (population != GROWER_POPULATION_GEN10) {
//            printf("ERROR: Incorrect population for generation 10\n");
//            exit(1);
//        }
//    }
//
//    if (iter == 100) {
//        int population = board_population(board);
//        printf("Generation 100: Grower pattern population = %d\n", population);
//        if (population != GROWER_POPULATION_GEN100) {
//            printf("ERROR: Incorrect population for generation 100\n");
//            exit(1);
//        }
//    }
//}

int main(int argc, char *argv[]) {
    int rank, size;

    // Start up MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const bool am_master = 0 == rank;

    // Divide the board into blocks for parallel processing (horizontal segments)
    int block_size = BOARD_SIZE / size;
    int block_start = rank * block_size;
    int block_end = block_start + block_size;
    // Print the block start and end indices for each process
    printf("Rank %d: Block start = %d, Block end = %d\n", rank, block_start, block_end);

    // Allocate memory for the local and next generation boards
    uint8_t *local_board = (uint8_t *) malloc(BOARD_SIZE * (block_size + 2) * sizeof(uint8_t));
    uint8_t *next_gen_board = (uint8_t *) malloc(BOARD_SIZE * (block_size + 2) * sizeof(uint8_t));

    // Initialize the local board with a specific pattern
    initialize_local_board(local_board, block_start, block_size);

    // Record the start time for measuring execution time
    double start_time = MPI_Wtime();

    for (int iter = 0; iter < MAX_ITERATIONS; ++iter) {
        printf("Rank %d: Iteration %d, Population = %d\n", rank, iter, board_population(local_board, block_size));

//        print_region(local_board, GROWER_START_ROW, GROWER_START_COL, GROWER_HEIGHT, GROWER_WIDTH);
        // Verification checks for specific iterations
//        run_verification(local_board, iter);

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
            MPI_Isend(&local_board[INDEX(1, 0)], BOARD_SIZE, MPI_UINT8_T, left_neighbour, 0, MPI_COMM_WORLD,
                      &left_send_request);
            MPI_Irecv(&local_board[INDEX(0, 0)], BOARD_SIZE, MPI_UINT8_T, left_neighbour, 0, MPI_COMM_WORLD,
                      &left_recv_request);
        }

        if (block_end < BOARD_SIZE) {
            // Communication with the right neighbor
            MPI_Isend(&local_board[INDEX(block_size, 0)], BOARD_SIZE, MPI_UINT8_T, right_neighbour, 0, MPI_COMM_WORLD,
                      &right_send_request);
            MPI_Irecv(&local_board[INDEX(block_size + 1, 0)], BOARD_SIZE, MPI_UINT8_T, right_neighbour, 0,
                      MPI_COMM_WORLD, &right_recv_request);
        }

        // Wait for completion of all non-blocking communication operations
        if (block_start > 0) {
            MPI_Wait(&left_send_request, MPI_STATUS_IGNORE);
            MPI_Wait(&left_recv_request, MPI_STATUS_IGNORE);
        }
        if (block_end < BOARD_SIZE) {
            MPI_Wait(&right_send_request, MPI_STATUS_IGNORE);
            MPI_Wait(&right_recv_request, MPI_STATUS_IGNORE);
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
//    printf("Total execution time: %f seconds\n", end_time - start_time);

    // Print the final time : the time taken by the slowest process to complete the execution
    double max_time;
    MPI_Reduce(&end_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (am_master) {
        printf("Max execution time: %f seconds\n", max_time - start_time);
    }


    // Free allocated memory
    free(local_board);
    free(next_gen_board);

    // Finalize MPI
    MPI_Finalize();

    return 0;
}
