typedef struct Graph{
    double **vertices;
    double **wam_mat;
    double **lnorm_mat;
    double **ddg_mat;
    int n, dim;
} Graph;

typedef enum Goal
{
    e_wam = (int)'w',
    e_ddg = (int)'d',
    e_lnorm = (int)'l',
    e_jacobi = (int)'j',
    e_spk = (int)'s'
} Goal;

double distance(double vector1[], double vector2[], int length);
void assertion_check(int b);
void read_data(Graph * graph, char * file_path);
double ** make_mat(int m, int n);
double ** make_mat_identity(int m, int n);
double * make_vector(int n);
void print_mat(double ** mat, int i, int j);
void print_array(double * array, int n);
void free_mat(double ** mat, int x);
void wam(Graph * graph);
void copy_array(double * old, double * copy, int culs);
void ddg(Graph * graph);
double ** sqrt_diagonal_matrix(double ** diagonal_matrix, int n);
void lnorm(Graph * graph);
void find_largest_element_pivot(double ** matrix, int n, int * piv);
int sign(double num);
double off(double ** matrix, int n);
void iter_jacobi(double ** matrix, double ** matrix_tag, double **eigenvectors,
                        double **eigenvectors_tag, int n, int i, int j);
void update_matrix(double ** matrix, double ** matrix_tag, int n,
                                                        int i, int j);
int check_convergence(double ** matrix, double ** matrix_tag, int n);
void jacobi(double ** matrix, int n, double ** eigenvectors,
                                                double * eigenvalues);
void sort_descending(double * array, int n);
int largest_k_eigenvectors(double * eigenvalues, int n);
void make_U(double ** eigenvectors, double * eigenvalues,
                    double * eigenvalues_sorted, double ** U, int N, int K);
void make_T(double ** U, double ** T, int N, int K);
int * make_cluster_count(int K);
int get_closest_cluster(double *vector, double **centroids_array,
                                            int vector_size, int K);
void add_vector_to_cluster(double *vector, double *cluster_sum,
                                                    int vector_size);
int update_medians(double **cluster_sum, double **centroids_array,
                        int *cluster_count, int vector_size, int K);
int k_means_c(double **vectors_array, double **centroids_array, int dimension,
                                    int num_of_vectors, int k, int max_iter);