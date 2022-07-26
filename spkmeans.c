#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

/* returns the distance between two vectors */
static double distance(double vector1[], double vector2[], int length){
    double sum=0.0;
    int i;
    for(i=0; i < length; i++){
        sum += (vector1[i] - vector2[i])*(vector1[i] - vector2[i]);
    }
    return sum;
}

int get_vector_size(FILE * file){
    char line[1000];
    int counter = 1;
    int i;
    fseek(file, 0, SEEK_SET); 
    /* get to the start of the file */
    if(fgets(line, 1000, file)){
        for (i = 0 ; i < 1000 ; i++){
            if (line[i] == ','){
                counter++;
            }
        }
    }
    return counter;
}

int get_num_of_vectors(FILE * file){
    char line[1000];
    int counter = 0;
    fseek(file, 0, SEEK_SET);
    /* get to the start of the file */
    while (fgets(line, 1000, file)){
       counter++;
    }
    return counter;
}

double ** make_vectors_array(FILE * file, int vector_size, int num_of_vectors){
    int i=0,j=0;
    double **vectors_array;
    char trash;
    fseek(file, 0, SEEK_SET);  /* get to the start of the file */
    vectors_array = calloc(num_of_vectors, sizeof(double*));
    for(i=0 ; i < num_of_vectors ; i++){
        vectors_array[i] = calloc(vector_size, sizeof(double));
        for(j=0 ; j<vector_size ; j++){
            fscanf(file, "%lf", &vectors_array[i][j]);
            fscanf(file, "%c", &trash);
        }
    }
    return vectors_array;
}

static double ** make_wam_mat(int n){
    double ** wam_mat;
    int i,j;
    wam_mat = calloc(n+1, sizeof(double*));
    for (i = 0; i < n; i++){
        wam_mat[i] = calloc(n+1, sizeof(double));
        for (j=0 ; j < n ; j++){
            wam_mat[i][j] = 0.0;
        }
    }
    return wam_mat;
}

static void print_mat(double ** mat, int i, int j){
    int x,y;
    for (x = 0 ; x < i ; x++){
        for (y = 0 ; y < j ; y++){
            if (y == j-1){
                printf("%.4f", mat[i][j]);
            }
            else{
                printf("%.4f,", mat[i][j]);
            }
        }
        printf("\n");
    }
}

static double ** wam(double ** vectors, int d, int n, char * goal){
    double ** wam_mat = make_wam_mat(n);
    double weight;
    int i,j;
    for (i = 0; i < n; i++){
        for (j=i+1 ; j < n ; j++){
            weight = exp(-sqrt(distance(vectors[i], vectors[j], d)));
            wam_mat[i][j] = weight;
            wam_mat[j][i] = weight;
        }
    }

    if (goal == "wam"){
        print_mat(wam_mat, n ,n);
    }
    else{
        
    }
}

static double ** ddg(){}

static double ** lnorm(){}

static double ** jacobi(){}


int main(int argc, char *argv[]){
    FILE *input_file;
    int num_of_vectors, vector_size, *cluster_count;
    int i, K, max_iter, iteration, closest_cluster;
    char *input_name, *goal;
    double **vectors_array, **centroids_array, **cluster_sum;

    /* checking correctness of the input */
    if(argc != 3 ){        
        printf("Invalid Input!");
        return 1;
    }

    input_name = argv[1];
    goal = argv[2];
 
    if((input_file = fopen(input_name, "r")) == NULL){
        printf("Invalid Input!");
        return 1;
    }

    vector_size = get_vector_size(input_file);
    num_of_vectors = get_num_of_vectors(input_file);
    vectors_array = make_vectors_array(input_file, vector_size, num_of_vectors);
}