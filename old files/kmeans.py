import sys
import math

def error(mes):
    sys.exit("Invalid Input!") if mes == 1 else sys.exit("An Error Has Occurred")

def distance(vector1, vector2):
    sum = 0.0
    for i in range(vector_size):
        sum += (vector1[i] - vector2[i])*(vector1[i] - vector2[i])
    return sum

def get_closest_cluster(vector, centroids_list):
    closest_cluster = -1
    min_dist = float('inf')
    for j in range(K):
        dist = distance(vector, centroids_list[j])
        if min_dist > dist:
            min_dist = dist
            closest_cluster = j
    return closest_cluster

def add_vector_to_cluster(vector, cluster_sum):
    for j in range(vector_size):
        cluster_sum[j] += vector[j]

def update_medians(cluster_sum, centroids_list, cluster_count):
    epsilon = 0.001
    epsilon_bool = 0
    for i in range(K):
        norma = 0
        for j in range(vector_size):
            if cluster_count[i] == 0:
                error(2)
            median = cluster_sum[i][j]/cluster_count[i]
            norma += (median - centroids_list[i][j]) * (median - centroids_list[i][j])
            centroids_list[i][j] = median
            cluster_sum[i][j] = 0
        cluster_count[i] = 0
        if math.sqrt(norma) >= epsilon:
            epsilon_bool = 1
    return epsilon_bool

def write_on_file(centroids_list, output_name):
    file = ""
    try:
        file = open(output_name, "w")
        for i in range(K):
            for j in range(vector_size):
                if j < vector_size - 1:
                    file.write(str(format(centroids_list[i][j], ".4f")) + ",")
                else:
                    file.write(str(format(centroids_list[i][j], ".4f")) + "\n")
    except:
        error(2)
    file.close()




K = int(sys.argv[1]) if sys.argv[1].isdigit() else error(1)
if len(sys.argv) <= 3 or len(sys.argv) >=6:
    error(1)
elif len(sys.argv) == 5:
    max_iter = int(sys.argv[2]) if sys.argv[2].isdigit() else error(1)
    input_name = sys.argv[3]
    output_name = sys.argv[4]
else:
    max_iter = 200
    input_name = sys.argv[2]
    output_name = sys.argv[3]

input_file, vector_size, num_of_vectors, vectors_list, centroids_list = "", 0, 0, [], []
try:
    input_file = open(input_name, "r")
    for line in input_file:
        if num_of_vectors == 0:
            vector_size = line.count(",") + 1 if line else error(2)
        line = line[:-1].split(sep=",")
        line = [float(item) for item in line]
        num_of_vectors += 1
        if num_of_vectors <= K:
            centroids_list.append(line.copy())
        vectors_list.append(line)

except:
    error(2)

input_file.close()

cluster_sum = [[0 for j in range(vector_size)] for i in range(K)]
cluster_count = [0 for i in range(K)]

iteration = 0
while iteration < max_iter:
    for i in range(num_of_vectors):
        closest_cluster = get_closest_cluster(vectors_list[i], centroids_list)
        add_vector_to_cluster(vectors_list[i], cluster_sum[closest_cluster])
        cluster_count[closest_cluster] += 1
    if not update_medians(cluster_sum, centroids_list, cluster_count):
        break
    iteration += 1

if write_on_file(centroids_list, output_name):
    error(2)
