/*******************************************************************************
Executing game of life for a parallel system.  This is a multi-processor project (i.e. distributed memory).  

During the k th generation, all cells determine their respective states based on the following simple rules: 
i) A living cell that has less than 2 living neighbors dies (as if caused by underpopulation, or loneliness, or simply boredom!); 
ii) A living cell with more than 3 living neighbors also dies (as if caused by overcrowding); 
iii) A living cell with either 2 or 3 living neighbors lives on to the next generation; 
iv) A dead cell with 3 living neighbors will be “born” (or “reborn” if it was ever alive before) and become alive.
********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <string.h>

#define ALIVE 1
#define DEAD 0
#define ROWS 8
#define COLUMNS 8
#define NUM_OF_TURNS 10
#define SAMPLESIZE 1

int m = ROWS;
int n = COLUMNS;

// *******************************************************************************
// Calculate number of rows in rankBoard needed
// NOTE!!!  adds 2 rows for ghost row on top and bottom!!!
// *******************************************************************************
int calcRowRankSize(int rank, int p)
{
    int rowRankSize = ((int)ROWS / p) + 2;
    if (rank == p - 1)
        rowRankSize += ROWS % p;
    return rowRankSize;
}

// *******************************************************************************
// Randomly populate the GoL board.  Obviously, this does not follow the the GOL 
// rules stated above
// *******************************************************************************
void populateGameBoard(char board[][COLUMNS], int p)
{
    int i, j;
    char c;
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            if (rand() % 2 == 0)
                c = ' ';
            else
                c = 'O';
            board[i][j] = c;
        }
    }
}

// *******************************************************************************
// Copy over relevants segments of GoL board into distributed memory
// *******************************************************************************
void populateRankBoard(char board[][COLUMNS], char rankBoard[][COLUMNS], int rank, int p)
{

    int rowRankSize = calcRowRankSize(rank, p);

    int boardRow = rank * ((int)ROWS / p);

    // copying over from board to rankBoard
    int rankRow, j;
    for (rankRow = 1; rankRow < rowRankSize; rankRow++, boardRow++)
    {
        for (j = 0; j < n; j++)
        {
            rankBoard[rankRow][j] = board[boardRow][j];
        }
    }
}

void copyBoard(char rankBoard[][COLUMNS], char tempBoard[][COLUMNS], int rank, int p)
{
    int rowRankSize = calcRowRankSize(rank, p);
    int i, j;
    for (i = 0; i < rowRankSize; i++)
    {
        for (j = 0; j < COLUMNS; j++)
        {
            rankBoard[i][j] = tempBoard[i][j];
        }
    }
}

// *******************************************************************************
// Helper function :: determines if cells lives or dies based on current state and neighbors
// *******************************************************************************
char executeGOLHelper(int liveCells, char center)
{
    char c = ' ';
    if (liveCells == 3 && center == ' ')
        c = 'O';
    else if ((liveCells == 2 || liveCells == 3) && center == 'O')
        c = 'O';

    return c;
}

// *******************************************************************************
// evaluate rankBoard but have the answers stored in tempBoard
// *******************************************************************************
void executeGOL(char rankBoard[][COLUMNS], char tempBoard[][COLUMNS], int rank, int p)
{
    int rowRankSize = calcRowRankSize(rank, p);
    int i, j;
    for (i = 1; i < rowRankSize - 1; i++)
    {
        for (j = 0; j < COLUMNS; j++)
        {
            int liveCells = 0;

            if (rankBoard[i - 1][(j - 1 + COLUMNS) % COLUMNS] == 'O') // upper left
                liveCells++;
            if (rankBoard[i - 1][j] == 'O') // upper center
                liveCells++;
            if (rankBoard[i - 1][(j + 1) % COLUMNS] == 'O') // upper right
                liveCells++;
            if (rankBoard[i][(j - 1 + COLUMNS) % COLUMNS] == 'O') // middle left
                liveCells++;
            if (rankBoard[i][(j + 1) % COLUMNS] == 'O') // middle right
                liveCells++;
            if (rankBoard[i + 1][(j - 1 + COLUMNS) % COLUMNS] == 'O') // lower left
                liveCells++;
            if (rankBoard[i + 1][j] == 'O') // lower center
                liveCells++;
            if (rankBoard[i + 1][(j + 1) % COLUMNS] == 'O') // lower right
                liveCells++;
            tempBoard[i][j] = executeGOLHelper(liveCells, rankBoard[i][j]);
        }
    }
    copyBoard(rankBoard, tempBoard, rank, p);
}

// *******************************************************************************
// print GOL board
// *******************************************************************************
void printBoard(char rankBoard[][COLUMNS], int rank, int p)
{
    int rowRankSize = calcRowRankSize(rank, p);
    int baseRowSize = calcRowRankSize(0, p);
    int i, j;
    for (i = 1; i < rowRankSize - 1; i++)
    {
        printf("row %d \t|", rank * (baseRowSize - 2) + i - 1);
        for (j = 0; j < COLUMNS; j++)
            printf("%c", rankBoard[i][j]);
        printf("|\n");
    }
}

// *******************************************************************************
// Sanity check to make sure initial board seems reasonable
// *******************************************************************************
void printInitBoard(char rankBoard[][COLUMNS], int rowSize)
{
    printf("============================================================\n");
    printf("Init population\n");
    int i, j;

    for (i = 0; i < rowSize; i++)
    {
        printf("row %d \t|", i);
        for (j = 0; j < COLUMNS; j++)
            printf("%c", rankBoard[i][j]);
        printf("|\n");
    }
}

// *******************************************************************************
// Initialize Rank board
// *******************************************************************************
void clearRankBoard(char rankBoard[][COLUMNS], int rank, int p)
{
    int rowRankSize = calcRowRankSize(rank, p);
    int i, j;

    for (i = 0; i < rowRankSize; i++)
    {
        for (j = 0; j < COLUMNS; j++)
            rankBoard[i][j] = ' ';
    }
}

double calculateAverageTime(double timeElapsed[])
{
    int i;
    double sum = 0;
    for (i = 0; i < SAMPLESIZE; i++)
    {
        sum += timeElapsed[i];
    }
    return sum / SAMPLESIZE;
}

int main(int argc, char *argv[])
{
    char c;
    int rank, p;
    int i, j, k;

    srand((unsigned)time(NULL));
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    int rowRankSize = calcRowRankSize(rank, p);

    char board[m][n];
    char rankBoard[rowRankSize][n]; // +2 needed for "ghost" row above and below
    char tempBoard[rowRankSize][n];
    double timeElapsed[SAMPLESIZE];
    double comTimeElapsed[SAMPLESIZE];

    for (k = 0; k < SAMPLESIZE; k++)
    {
        clearRankBoard(rankBoard, rank, p);
        clearRankBoard(tempBoard, rank, p);

        populateGameBoard(board, p);
        populateRankBoard(board, rankBoard, rank, p);

        int right = (rank + 1) % p;
        int left = (rank - 1 + p) % p;

        double totalPrintingTime = 0;
        double totalTime = omp_get_wtime();

        for (i = 0; i < NUM_OF_TURNS; i++)
        {
            // this should synchronize
            MPI_Barrier(MPI_COMM_WORLD);
            
            // sending top row active of rank N to bottom ghost row of rank N-1
            char messageSend[COLUMNS];
            char messageRecv[COLUMNS];
            memcpy(messageSend, rankBoard[1], COLUMNS);
            double printingTime = omp_get_wtime();
            MPI_Sendrecv(messageSend, COLUMNS, MPI_CHAR, left, 123, messageRecv, COLUMNS, MPI_CHAR, right, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printingTime = omp_get_wtime() - printingTime;
            totalPrintingTime += printingTime;
            memcpy(rankBoard[rowRankSize - 1], messageRecv, COLUMNS);

            // sending bottom row active of rank N to top ghost row of rank N+1
            memcpy(messageSend, rankBoard[rowRankSize - 2], COLUMNS);
            printingTime = omp_get_wtime();
            MPI_Sendrecv(messageSend, COLUMNS, MPI_CHAR, right, 123, messageRecv, COLUMNS, MPI_CHAR, left, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printingTime = omp_get_wtime() - printingTime;
            totalPrintingTime += printingTime;
            memcpy(rankBoard[0], messageRecv, COLUMNS);

            // evaluate rankBoard but have the answers stored in tempBoard.  Set tempBoard as new rankBoard
            executeGOL(rankBoard, tempBoard, rank, p);

            // print board and send message to next rank to allow them to print board
            if (rank == 0)
            {
                printBoard(rankBoard, rank, p);
                int dummy = 1;
                MPI_Send(&dummy, 1, MPI_INT, right, 123, MPI_COMM_WORLD);
            }
            else if (rank == p - 1)
            {
                int dummy = 1;
                MPI_Recv(&dummy, 1, MPI_INT, left, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printBoard(rankBoard, rank, p);
            }
            else
            {
                int dummy = 1;
                MPI_Recv(&dummy, 1, MPI_INT, left, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                printBoard(rankBoard, rank, p);
                MPI_Send(&dummy, 1, MPI_INT, right, 123, MPI_COMM_WORLD);
            }
        }

        totalTime = omp_get_wtime() - totalTime - totalPrintingTime;
        if (rank == 0)
        {
            timeElapsed[k] = totalTime;
            comTimeElapsed[k] = totalPrintingTime;
        }
    }
    if (rank == 0)
    {
        double averageTime = calculateAverageTime(timeElapsed);
        double averageComTime = calculateAverageTime(comTimeElapsed);
        printf("Total samples = %d, Procs = %d, size = %d, Average Computational Time = %lf, Average Communation Time = %lf\n",
               SAMPLESIZE, p, COLUMNS, averageTime, averageComTime);
    }
    MPI_Finalize();
}
