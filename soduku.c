#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <omp.h>

// Here the size of the bitmask can be changed
typedef int bitmask;
void lone_rangers(int** board, const int M, const int block_M, bitmask** possible_vals, int n_threads, int solved);
void elimination(int** board, const int M, const int block_M, bitmask** possible_vals, int n_threads, int solved);
void print_board(int** board, const int M, const int block_M);

int solved = 0;
double start_time = 0;
double end_time = 0;
int brute_force_called = 0;
void read_sudoku(char* filename, int** board, const int M){
    FILE *file = fopen(filename, "r");
    if (!file) {
        printf("Unable to open file\n");
        return;
    }
    printf("File opened successfully.\n");
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            int r = fscanf(file, "%d,", &board[i][j]);
        }
    }
    fclose(file);
}

// Function to count the number of values in a bitmask
int num_possible_vals(bitmask vals){
    int count = 0;
    // While there still are values in the bitmask
    while (vals) {
        // Remove the current value and increase count
        vals &= vals - 1;
        count++;
    }
    return count;
}

// Function to check if all cells are filled in
int is_solved(int** board, const int M){
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (board[i][j] == 0)
            {
                return 0;
            }
        }
    }
    return 1;
}

// Function to check if a number is in a soduku minigrid
int in_block(int** board, int i_start, int j_start, const int block_M, int num){
    for (int i = 0; i < block_M; i++)
    {
        for (int j = 0; j < block_M; j++)
        {
           if (board[i_start+i][j_start+j] == num)
           {
            return 1;
           }
        }
    }
    return 0;
    
}

// Function to check if a number is in a row
int in_row(int** board, int i, int num, const int M){
    for (int j = 0; j < M; j++)
    {
        if (board[i][j]==num)
        {
            return 1;
        }
    }
    return 0;
}

// Function to check if a number is in a column
int in_col(int** board, int j, int num, const int M){
    for (int i = 0; i < M; i++)
    {
        if (board[i][j]==num)
        {
            return 1;
        }
    }
    return 0;
}

// Function to compute possible values for a cell and return a bitmask containing 
// which numbers are possible
inline bitmask possible_val(int** board, int i, int j, const int M, const int block_M, int n_threads){
    // Initialize the bitmask to all ones
    bitmask possible_vals = (1ULL << M) - 1;

    // Iterate through each possible number
    for (int k = 1; k <= M; k++)
    {
        // Check if that number is the cells soduku minigrid, row or column
        if (in_block(board, i-i%block_M, j-j%block_M, block_M, k) || in_row(board, i, k, M) || in_col(board, j, k, M))
        {
            // Set the index corresponding to that number to zero
            // This indicates that this number is not possible
            possible_vals &= ~(1ULL << (k - 1));
        }
    }
    return possible_vals;
}

// Simple function that returns if a number is in a bitmask
bitmask is_value_possible(bitmask possible_vals, int num) {
    return (possible_vals & (1ULL << (num - 1)));
}

// Function that updates all bitmask of all cells of the possible values
void update_possible(int** board, const int M, const int block_M, bitmask** possible_vals, int n_threads){
    // Iterate through the board
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            // Skip empty cells
            if (board[i][j]==0)
            {
            // Check which numbers are possible and assign that bitmask to the cell
            possible_vals[i][j] = possible_val(board, i, j, M, block_M, n_threads);
            }
            else{
            possible_vals[i][j] = 0;
            }
        }
    }
}

// Function to remove a number from each row, minigrid and column corresponding
// to a specific cell
void remove_possible(int** board, int i, int j, const int M, const int block_M, bitmask** possible_vals, int n_threads, int num){
    // remove the inserted value from row and column
    possible_vals[i][j]=0;
    for (int index = 0; index < M; index++)
    {
        if (board[i][index] == 0)
        {
            possible_vals[i][index] &= ~(1ULL << (num-1));
        }
        if (board[index][j] == 0)
        {
            possible_vals[index][j] &= ~(1ULL << (num-1));
        }
    
    }
    // remove the inserted value from block
    int istart = i-i%block_M; int jstart = j-j%block_M;
    for (int row = 0; row < block_M; row++)
    {
        for (int col = 0; col < block_M; col++)
        {
            if (board[istart + row][jstart + col] == 0)
            {
                possible_vals[istart + row][jstart + col] &= ~(1ULL << (num-1));
            }  
        }
    }
}

// function to check if a board is solved and if the solution is valid
int check_valid_board(int** board, const int M, const int block_M){
    if (!is_solved(board, M))
    {
        return 0;
    }
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            int num = board[i][j];
            // check the row
            for (int col = 0; col < M; col++){
                if (col != j && board[i][col] == num){
                    return 0;
                }
            }
            // check the column
            for (int row = 0; row < M; row++){
                if (row != i && board[row][j] == num){
                    return 0;
                }
            }
            // check the block
            for (int blocki = 0; blocki < block_M; blocki++){
                for (int blockj = 0; blockj < block_M; blockj++){
                    int indexi = blocki-blocki%block_M+i;
                    int indexj = blockj-blockj%block_M+j;
                    if ((indexi != i && indexj!= j) && board[indexi][indexj] == num){
                        return 0;
                    }
                }
            }

        }
    }
    return 1;
}

// Function to detect the cell with lowest amount of possible values
void lowest_pos_vals(int** board, bitmask** possible_vals, int* row, int* col, const int M){
    // Initialize the count to the max value
    int count = M+1;
    // Iterate through each cell
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            // Skip empty cells
            if (board[i][j]==0)
            {
            // If the count is larger than the current cells number
            // of possible values update the indexes and count
            if (count > num_possible_vals(possible_vals[i][j]))
            {
                count = num_possible_vals(possible_vals[i][j]);
                *row = i;
                *col = j;
            }
            }
        }
    }
}

// Function to count the number of empty cells
int count_empty(int** board, const int M){
    int count = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (board[i][j]==0)
            {
                count++;
            }
        }
    }
    return count;
}

void brute_force(int*** board, const int M, const int block_M, bitmask*** possible_vals, int n_threads, int empty_cells){
    
    if (is_solved(*board, M))
    {
        solved = 1;
        if(check_valid_board(*board, M, block_M)){
        printf("\n\n");
        print_board(*board, M, block_M);
        printf("The solution is valid!\n");
        return;
        }
        else{
        printf("\n\n");
        print_board(*board, M, block_M);
        printf("The solution is not valid!\n");
        return;
        }
    }
    if(solved){
        return;
    }
    brute_force_called = 1;
    // Intialize the indexes of the cell with smallets amount possible values
    int mini = -1; int minj = -1;
    // Assign the cell index to the integers
    lowest_pos_vals(*board, *possible_vals, &mini, &minj, M);
    // Define a taskgroup
    #pragma omp taskgroup
    {
    
    // Iterate through the numbers
    for (int k = 1; k <= M; k++){
        // only proceed if the value is possible to assign
        if(is_value_possible((*possible_vals)[mini][minj],k)){
            // if there is only one option it can be assigned directly
            if (num_possible_vals((*possible_vals)[mini][minj])==1)
            {
                (*board)[mini][minj] = k;
                empty_cells--;
                remove_possible(*board, mini, minj, M, block_M, *possible_vals, n_threads, k);
                brute_force(board, M, block_M, possible_vals, n_threads, empty_cells);
            }
            else{
            // check if the number of empty cells is a high enough number
            // to create new tasks
            if (empty_cells > M/block_M)
            {
            // assign each possible number to explore to a thread
            #pragma omp task
            {
            // ensure no other thread has solved the soduku
            #pragma omp flush(solved)
            // only proceed if no sultion has been found
            if (!solved)
            {
            #pragma omp cancellation point taskgroup

            // create copies of the board and possible values local to each thread
            int** board_copy = (int**)malloc(M * sizeof(int*));
            bitmask** possible_vals_copy = (bitmask**)malloc(M * sizeof(bitmask*));
            for (int i = 0; i < M; i++) {
                board_copy[i] = (int*)malloc(M * sizeof(int));
                memcpy(board_copy[i], (*board)[i], M*sizeof(int));
                possible_vals_copy[i] = (bitmask*)malloc(M * sizeof(bitmask));
                memcpy(possible_vals_copy[i], (*possible_vals)[i], M*sizeof(bitmask));
            }
            // assign the number and update necessary cells
            board_copy[mini][minj] = k;
            empty_cells--;
            remove_possible(board_copy, mini, minj, M, block_M, possible_vals_copy, n_threads, k);
            // recursive call to brute force to explore this branch fully
            brute_force(&board_copy, M, block_M, &possible_vals_copy, n_threads, empty_cells);
            // if solution found cancel the tasks
            if (solved)
            {
                #pragma omp cancel taskgroup
            }
            
            // free the local copies
            for (int i = 0; i < M; i++) {
                    free(board_copy[i]);
                    free(possible_vals_copy[i]);
                }
            free(board_copy);
            free(possible_vals_copy);
            #pragma omp cancellation point taskgroup
            }
            }
            }
            // same brute force function but with no task creation since the tasks
            // are small
            else{
            int** board_copy = (int**)malloc(M * sizeof(int*));
            bitmask** possible_vals_copy = (bitmask**)malloc(M * sizeof(bitmask*));
            for (int i = 0; i < M; i++) {
                board_copy[i] = (int*)malloc(M * sizeof(int));
                memcpy(board_copy[i], (*board)[i], M*sizeof(int));
                possible_vals_copy[i] = (bitmask*)malloc(M * sizeof(bitmask));
                memcpy(possible_vals_copy[i], (*possible_vals)[i], M*sizeof(bitmask));
            }
            
            board_copy[mini][minj] = k;
            empty_cells--;
            remove_possible(board_copy, mini, minj, M, block_M, possible_vals_copy, n_threads, k);
            brute_force(&board_copy, M, block_M, &possible_vals_copy, n_threads, empty_cells);
            for (int i = 0; i < M; i++) {
                    free(board_copy[i]);
                    free(possible_vals_copy[i]);
                }
            free(board_copy);
            free(possible_vals_copy);
            }
            }
        }
        
    }
    }
return;
}

// functions to locate twins in the soduku board and update the board if they are found
void twins(int** board, const int M, const int block_M, bitmask** possible_vals, int n_threads, int solved){
    int change = 0;
    // first check blocks
    // iterate through each block
    for (int blocki = 0; blocki < M; blocki+=block_M){
        for (int blockj = 0; blockj < M; blockj+=block_M){
            for (int i = blocki; i < blocki+block_M; i++){
                for (int j = blockj; j < blockj+block_M; j++){
                    // if there exists exactly two possible values
                    if (num_possible_vals(possible_vals[i][j]) == 2){
                        // store the possible values for that index
                        bitmask vals = possible_vals[i][j];
                        // iterate through the block again
                        for (int k = blocki; k < blocki + block_M; k++){
                            for (int l = blockj; l < blockj + block_M; l++){
                                // exclude the current index
                               if (i != k && j != l){
                               // if the possible values match with another index
                               // they are twins
                               if (vals == possible_vals[k][l]){
                                // iterate through the block again
                                // to find matches of either of the twin values
                                for (int x = blocki; x < blocki + block_M; x++){
                                    for (int y = blockj; y < blockj + block_M; y++)
                                    {
                                        // skip the twins
                                        if ((i == x && j == y) || (k == x && l == y)) continue;

                                        // remove the values from possible values
                                        if (possible_vals[x][y] & vals)
                                        {
                                            possible_vals[x][y] &= ~vals;
                                            change = 1;
                                        }

                                    }
                                }
                               }
                            }
                            if (change){break;}
                            }
                        if (change){break;}
                        }
                    }
                if (change){break;}  
                }
            if (change){break;}   
            }
        if (change){break;}    
        }
    if (change){break;}
    }
    if (check_valid_board(board, M, block_M)){
        return;
    }
    if (change == 1){
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    // check rows
    // iterate through the rows
    for (int i = 0; i < M; i++) {
        // iterate through the cells in that row
        for (int j = 0; j < M; j++) {
           if (num_possible_vals(possible_vals[i][j]) == 2){
            bitmask vals = possible_vals[i][j];
            // search for a possible twin by iterating through the cells in the row
            for (int k = 0; k < M; k++){
                if (k != j) // exclude the stored cell
                {
                
                if (vals == possible_vals[i][k]){ // twin found
                    // iterate throguh the cells in the row again
                    for (int l = 0; l < M; l++){
                        if (l != j && l != k) { // Exclude the twin cells
                        if (possible_vals[i][l] & vals){
                        possible_vals[i][l] &= ~vals;
           
                        change = 1;
                        }
                        }
                    }
                }
                }
            if (change){break;}
            }
            }
        if (change){break;}
        }
    if (change){break;}
    }
    if (check_valid_board(board, M, block_M)){
        return;
    }
    if (change == 1){
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    // check columns
    // iterate through the columns
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < M; i++) {
           if (num_possible_vals(possible_vals[i][j]) == 2){
            bitmask vals = possible_vals[i][j];
            for (int k = 0; k < M; k++){
                if (k != i){ // exclude the stored cell
                if (vals == possible_vals[k][j]){
                    for (int l = 0; l < M; l++){
                        if (l != i && l != k) { // Exclude the twin cells
                        if (possible_vals[l][j] & vals){
                        possible_vals[l][j] &= ~vals;
                        change = 1;
                        }
                        }
                    }
                }
                }
            if (change){break;}
            }
            }
        if (change){break;}
        }
    if (change){break;}
    }
    if (check_valid_board(board, M, block_M)){
        return;
    }
    if (change == 1){
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    update_possible(board, M, block_M, possible_vals, n_threads);
    int count = count_empty(board, M);
    #pragma omp parallel num_threads(n_threads)
    {
    #pragma omp single
    {
    brute_force(&board, M, block_M, &possible_vals, n_threads, count);
    }
    }
}

// functions to identify lone rangers and assign values
void lone_rangers(int** board, const int M, const int block_M, bitmask** possible_vals, int n_threads, int solved){
    int change = 0;

    // temporary arrays to store counts of each number and position if lone ranger found
    int* count = (int*)calloc(M, sizeof(int));
    int* pos = (int*)calloc(M, sizeof(int));
    int* pos2 = (int*)calloc(M, sizeof(int));
    // check blocks
    for (int blocki = 0; blocki < M; blocki+=block_M){
        for (int blockj = 0; blockj < M; blockj+=block_M){
            // initalize arrays
            memset(count, 0, M * sizeof(int));
            memset(pos, -1, M * sizeof(int));
            memset(pos2, -1, M * sizeof(int));
            // iterate through the block
            for (int i = blocki; i < blocki+block_M; i++){
                for (int j = blockj; j < blockj+block_M; j++){
                    // skip empty cells
                    if (board[i][j]==0){
                        // loop to count the occurrence of each number in the block
                        // and store the index for when we only have one count
                        for (int k = 1; k <= M; k++){
                            if (is_value_possible(possible_vals[i][j], k)){
                                count[k-1]++;
                                pos[k-1] = i;
                                pos2[k-1] = j;
                            }  
                        }
                    }  
                }   
            }    
            // if there is only one occurrence of each number in the block
            // that number is lone ranger
            for (int k = 1; k <= M; k++){
                if (count[k-1]==1){
                    board[pos[k-1]][pos2[k-1]] = k;
                    remove_possible(board, pos[k-1], pos2[k-1], M, block_M, possible_vals, n_threads, k);
                    change = 1;
                    break;
                    }
                }
        if (change){break;}
        }
    if (change){break;}
    }
    // if soduku solved we abort the algorithm
    if (check_valid_board(board, M, block_M)){
        free(count);
        free(pos);
        free(pos2);
        return;
    }
    // if a change has been made we restart the algorithm
    if (change == 1){
        free(count);
        free(pos);
        free(pos2);
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    
    // same logic as before but for rows
    for (int i = 0; i < M; i++) {
        memset(count, 0, M * sizeof(int));
        memset(pos, -1, M * sizeof(int));
        for (int j = 0; j < M; j++) {
            if (board[i][j] == 0) {
                for (int k = 1; k <= M; k++) {
                    if (is_value_possible(possible_vals[i][j], k)) {
                        count[k - 1]++;
                        pos[k - 1] = j;
                    }
                }
            }
        }
        for (int k = 1; k <= M; k++) {
            if (count[k - 1] == 1) {
                board[i][pos[k-1]] = k;
                remove_possible(board, i, pos[k-1], M, block_M, possible_vals, n_threads, k);
                change = 1;
                break;
            }
        }
    if (change){break;}
    }
    if (check_valid_board(board, M, block_M)){
        free(count);
        free(pos);
        free(pos2);
        return;
    } 
    if (change == 1){
        free(count);
        free(pos);
        free(pos2);
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    
    // same logic as before but with columns
    for (int j = 0; j < M; j++) {
        memset(count, 0, M * sizeof(int));
        memset(pos, -1, M * sizeof(int));
        for (int i = 0; i < M; i++) {
            if (board[i][j] == 0) {
                for (int k = 1; k <= M; k++) {
                    if (is_value_possible(possible_vals[i][j], k)) {
                        count[k - 1]++;
                        pos[k - 1] = i;
                    }
                }
            }
        }
        for (int k = 1; k <= M; k++) {
            if (count[k - 1] == 1) {
                board[pos[k-1]][j] = k;
                //update_possible(board, M, block_M, possible_vals, n_threads);
                remove_possible(board, pos[k-1], j, M, block_M, possible_vals, n_threads, k);
                change = 1;
                break;
                
            }
        if (change){break;}
        }
    if (change){break;}
    }
    
    if (check_valid_board(board, M, block_M)){
        free(count);
        free(pos);
        free(pos2);
        return;
    }
    if (change == 1){
        free(count);
        free(pos);
        free(pos2);
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    free(count);
    free(pos);
    free(pos2);
    twins(board, M, block_M, possible_vals, n_threads, solved);
}


// function to find cells with only one possible value and assign this value
void elimination(int** board, const int M, const int block_M, bitmask** possible_vals, int n_threads, int solved){
    int change = 0;
    for (int i = 0; i < M; i++)
    {
        for (int j = 0; j < M; j++)
        {
            if (num_possible_vals(possible_vals[i][j]) == 1){
            for (int k = 1; k <= M; k++)
            {
                if (is_value_possible(possible_vals[i][j], k))
                {
                    board[i][j] = k;
                    remove_possible(board, i, j, M, block_M, possible_vals, n_threads, k);
                    change = 1;
                    break;
                }
            }
            }
        if (change){break;}  
        }
   if (change){break;}
    }
    if (check_valid_board(board, M, block_M)){
        return;
    }
    if (change == 1){
        elimination(board, M, block_M, possible_vals, n_threads, solved);
        return;
    }
    lone_rangers(board, M, block_M, possible_vals, n_threads, solved);
    
}
void print_board(int** board, const int M, const int block_M){
    for (int i = 0; i < M; i++)
    {
        if (i%block_M == 0 && i != 0)
        {
            for (int j = 0; j < M*3+(block_M-1)*2-1; j++)
            {
                printf("-");
            }
            printf("\n");
        }
        
        for (int j = 0; j < M; j++)
        {
            if (j%block_M == 0 && j != 0)
            {
                printf("| ");
            }
            printf("%2d ", board[i][j]);
        }
        printf("\n");
    }
}
int main(int argc, char* argv[]){
    if (argc != 4) {
        fprintf(stderr, "Usage: %s M (dimension of MxM soduku board) file (.csv file of soduku board) n_threads \n", argv[0]);
        return 1;
    }
    // read size of soduku board
    const int M = atoi(argv[1]);
    const int block_M = sqrt(M);
    // read name of file containing soduku
    char* filename = argv[2];
    int n_threads = atoi(argv[3]);
    int solved = 0;

    // allocate memory for soduku board
    int** board = (int**)malloc(M * sizeof(int*));
    bitmask** possible_vals = (bitmask**)malloc(M * sizeof(bitmask*));
    for (int i = 0; i < M; i++) {
        board[i] = (int*)malloc(M * sizeof(int));
        possible_vals[i] = (bitmask*)malloc(M * sizeof(bitmask));
    }
    // read the Sudoku board
    read_sudoku(filename, board, M);
    // print the board
    print_board(board, M, block_M);
    // updates possible values
    update_possible(board, M, block_M, possible_vals, n_threads);
    // solve the board
    elimination(board, M, block_M, possible_vals, n_threads, solved);
    if(!brute_force_called){
        if (check_valid_board(board, M, block_M))
        {
            printf("\n\n");
            print_board(board, M, block_M);
            printf("The solution is valid!\n");
        }
        else{
            printf("\n\n");
            print_board(board, M, block_M);
            printf("The solution is not valid!\n");
        }
    }
    
    
    // Free allocated memory
    for (int i = 0; i < M; i++) {
        free(board[i]);
        free(possible_vals[i]);
    }
    free(board);
    free(possible_vals);
    return 0;
}