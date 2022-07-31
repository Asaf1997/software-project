import sys
import numpy as np
import pandas as pd


def error(b):
    sys.exit("Invalid Input!") if b else sys.exit("An Error Has Occurred")


def get_args():
    # checks if the arguments are valid
    argc = len(sys.argv)
    if argc != 4:
        error(True)
        return
    try:
        k = int(sys.argv[1])
        goal_str = str(sys.argv[2])
        input_name = str(sys.argv[3])
        if goal_str != "spk" and goal_str != "wam" and goal_str != "ddg" and goal_str != "lnorm" and goal_str != "jacobi":
            error(True)
            return

        return k, goal_str, input_name

    except:
        error(True)



# write the numpy format
float_formatter = "{:.4f}".format
np.set_printoptions(formatter={'float_kind': float_formatter})

# checking arguments and putting them into variables
k, goal_str, input_name = get_args()
vectors = "" # using C func

# if k == 0 then we need to use C func to bring k from L-Norm
if k == 0:
    pass

if k >= len(vectors):
    error(True)

initial_centroids = k_means_plus_plus(vectors, num_of_clusters)
observ_index = initial_centroids[:, 0:1].astype('int32')
initial_centroids = initial_centroids[:, 1:].tolist()
vectors = vectors[:,1:].tolist()

final_centroids = km.fit(len(vectors[0]), len(vectors), num_of_clusters, max_iter, eps, vectors, initial_centroids)

printSolution(final_centroids, observ_index)