#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include "spkmeans.h"


static void convert_list_to_c(PyObject* pyList, double ** array, int rows, int cull){
    int index_i, index_j;
    PyObject * py_row;
    PyObject * py_num;

    for (index_i = 0 ; index_i < rows ; index_i++){
        py_row = PyList_GetItem(pyList, index_i);
        assertion_check(!PyList_Check(py_row));

        for(index_j=0 ; index_j < cull ; index_j++){
            py_num = PyList_GetItem(py_row, index_j);
            assertion_check(!PyFloat_Check(py_num) && !PyLong_Check(py_num));
            array[index_i][index_j] = PyFloat_AsDouble(py_num);
        }
    }
}


static PyObject* convert_array_to_py(double** array, int rows, int cull){
    int index_i, index_j;

    PyObject * py_table = PyList_New(rows);
    PyObject * py_row;
    PyObject * py_num;

    for (index_i=0 ; index_i<rows ; index_i++){
        py_row = PyList_New(cull);
        for(index_j=0 ; index_j < cull ; index_j++){
            py_num = Py_BuildValue("d", array[index_i][index_j]);
            PyList_SetItem(py_row, index_j, py_num);
        }
        PyList_SetItem(py_table, index_i, py_row);
    }

    return py_table;
}


/* makes a matrix to print according to the jacobi program format as requested in the assignment */
static void make_jacobi_final_array(double ** eigenvectors, double * eigenvalues, double ** final_array, int n){
    int i,j;

    for (i = 0 ; i < n ; i++){
        final_array[0][i] = eigenvalues[i];
    }
    for(i=0; i < n ; i++){
       for (j=0 ; j < n ; j++){
            final_array[i+1][j] = eigenvectors[i][j];
       }
    }
}


/* returns a matrix (wam, ddg or lnorm) according to the goal requested */
static double ** compute_by_goal(Graph * graph, Goal goal){
    wam(graph);
    if (goal == e_wam){ return graph->wam_mat; }

    ddg(graph);
    if (goal == e_ddg){ return graph->ddg_mat; }

    lnorm(graph);
    return graph->lnorm_mat;
}


/* frees graph according to the goal requested */
static void free_graph(Graph * graph, Goal goal, int n){
    free_mat(graph->vertices, n);
    if (goal == e_jacobi){ return; }

    free_mat(graph->wam_mat, n);
    if (goal == e_wam){ return; }

    free_mat(graph->ddg_mat, n);
    if (goal == e_ddg){ return; }

    free_mat(graph->lnorm_mat, n);
}


/* the main function that is used in the python program */
static PyObject* goal_fit(PyObject* self, PyObject* args){
    int dimension, n, goal_num, K, rows, culs;
    double **vectors_array;
    double **final_array;
    double *spk_eigenvalues;
    double *spk_eigenvalues_sorted;
    double **spk_eigenvectors;
    double *jacobi_eigenvalues;
    double **jacobi_eigenvectors;
    double **U;
    double **T;
    Graph graph = {0};
    Goal goal;

    PyObject * py_vectors_array;
    PyObject * py_final_array = NULL;

    /* get arguments to variables */
    assertion_check(!PyArg_ParseTuple(args, "iiiiO", &dimension, &n, &goal_num, &K, &py_vectors_array));

    /* making c arrays with python input */
    assertion_check(!PyList_Check(py_vectors_array));
    vectors_array = make_mat(n, dimension);
    convert_list_to_c(py_vectors_array, vectors_array, n, dimension);

    graph.dim = dimension;
    graph.n = n;
    graph.vertices = vectors_array;

    goal = goal_num;

    rows = n;
    culs = n;

    if (goal != e_jacobi){
        final_array = compute_by_goal(&graph, goal);

        if (goal == e_spk){
            spk_eigenvectors = make_mat_identity(n ,n);
            spk_eigenvalues = make_vector(n);

            /* using jacobi algorithem on lnorm mat */
            jacobi(graph.lnorm_mat, n, spk_eigenvectors, spk_eigenvalues);

            /* making a sorted eigenvalues array */
            spk_eigenvalues_sorted = make_vector(n);
            copy_array(spk_eigenvalues, spk_eigenvalues_sorted, n);
            sort_descending(spk_eigenvalues_sorted, n);
                    
            /* if K == 0 then needs to check the ideal number of clusters */
            if (K == 0){
               K = largest_k_eigenvectors(spk_eigenvalues_sorted, n);
            }

            U = make_mat(n, K);
            T = make_mat(n, K);
            make_U(spk_eigenvectors, spk_eigenvalues, spk_eigenvalues_sorted, U, n, K);
            make_T(U, T, n, K);
            
            /* sending T matrix back to python to continue the process of kmeans */
            py_final_array = convert_array_to_py(T, rows, culs);
            culs = K;

            free_mat(U, n);
            free_mat(T, n);
            free_mat(spk_eigenvectors, n);
            free(spk_eigenvalues);
            free(spk_eigenvalues_sorted);
        }
        else{
            /* final array is wam, ddg or lnorm */
            py_final_array = convert_array_to_py(final_array, rows, culs);
        }
    }
    else{ /* goal == jacobi */
        jacobi_eigenvectors = make_mat_identity(n ,n);
        jacobi_eigenvalues = make_vector(n);

        jacobi(graph.vertices, n, jacobi_eigenvectors, jacobi_eigenvalues);
        final_array = make_mat(n+1, n);

        make_jacobi_final_array(jacobi_eigenvectors, jacobi_eigenvalues, final_array, n);
        rows = n+1;

        py_final_array = convert_array_to_py(final_array, rows, culs);

        free_mat(final_array, rows);
        free_mat(jacobi_eigenvectors, n);
        free(jacobi_eigenvalues);
    }

    free_graph(&graph, goal, n);
    return py_final_array;
}


/* kmeans algorithem for python program if goal is 'skp' */
static PyObject* fit_kmeans(PyObject* self, PyObject* args){
    int dimension, num_of_vectors, k, max_iter, i;
    double **vectors_array, **centroids_array;

    PyObject * py_vectors_array;
    PyObject * py_centroids_array;
    PyObject * py_final_centroids;

    /* get arguments to variables */
    if (!PyArg_ParseTuple(args, "iiiiOO", &dimension, &num_of_vectors, &k, &max_iter, &py_vectors_array, &py_centroids_array)){
        return NULL;
    }

    /* making c arrays with python input */
    vectors_array = make_mat(num_of_vectors, dimension);
    centroids_array = make_mat(k, dimension);
    convert_list_to_c(py_vectors_array, vectors_array, num_of_vectors, dimension);
    convert_list_to_c(py_centroids_array, centroids_array, k, dimension);

    /* initialize kmeans algorithem */
    assertion_check(!k_means_c(vectors_array, centroids_array, dimension, num_of_vectors, k, max_iter));

    py_final_centroids = convert_array_to_py(centroids_array, k, dimension);

    for(i=0 ; i < num_of_vectors ; i++){
        free(vectors_array[i]);
    }
    for(i=0 ; i < k ; i++){
        free(centroids_array[i]);
    }
    free(vectors_array);
    free(centroids_array);

    return py_final_centroids;
}


static PyMethodDef mykmeanssp_Methods[] = {
    {"goal_fit",
    (PyCFunction) goal_fit,
    METH_VARARGS,
    PyDoc_STR("if goal != spk")},
    {"fit_kmeans",
    (PyCFunction) fit_kmeans,
    METH_VARARGS,
    PyDoc_STR("if goal == spk")},
    {NULL,NULL,0,NULL}
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    mykmeanssp_Methods
};


PyMODINIT_FUNC
PyInit_mykmeanssp(void){
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m){
        return NULL;
    }
    return m;
}