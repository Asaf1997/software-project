typedef struct Graph
{
    double **vertices;
    double **wam_mat;
    double **lnorm_mat;
    double *ddg_mat;
    int n, dim;
} Graph;

typedef enum Goal
{
    e_wam = (int)'w',
    e_ddg = (int)'d',
    e_lnorm = (int)'l',
    e_jacobi = (int)'j'
} Goal;

static double distance(double vector1[], double vector2[], int length);