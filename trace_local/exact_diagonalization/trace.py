from numpy import dot, genfromtxt, copy, array, eye, trace
from numpy.linalg import eig
from math import exp
from scipy.linalg import expm

def mergesort_list(a, b, a_type, b_type):
    i, j = [0, 0]
    la = len(a)
    lb = len(b)
    merged_list = []
    merged_type = []
    if(la == 0):
        merged_list = copy(b)
        merged_type = copy(b_type)
        return merged_list, merged_type
    if(lb == 0):
        merged_list = copy(a)
        merged_type = copy(b_type)
        return merged_list, merged_type

    while (i < la or j < lb):
        if(i == la):
            merged_list.append(b[j])
            merged_type.append(b_type[j])
            j += 1
        elif(j == lb):
            merged_list.append(a[i])
            merged_type.append(a_type[i])
            i += 1
        elif (a[i] > b[j]):
            merged_list.append(b[j])
            merged_type.append(b_type[j])
            j += 1
        elif(a[i] <= b[j]):
            merged_list.append(a[i])
            merged_type.append(a_type[i])
            i += 1
    return merged_list, merged_type


def mergesort_multiple_lists(lists, types):
    merged_list = []
    merged_type = []
    for i in range(len(lists)):
        merged_list, merged_type = mergesort_list(merged_list, lists[i],
                                                  merged_type, types[i])
    return merged_list, merged_type

def expv(v):
    ev=eye(8)
    for i in range(8):
        ev[i][i]=exp(v[i])
    return ev

m = genfromtxt("matrix.txt")
v, u = eig(m)
ut = u.transpose()

beta = 10

cd_u = [1, 5]
c_u = [3, 7]
cd_d = [2]
c_d = [6]

cd_u_type = [0, 0]
c_u_type = [1, 1]
cd_d_type = [2]
c_d_type = [3]

lists = [cd_u, c_u, cd_d, c_d]
types = [cd_u_type, c_u_type, cd_d_type, c_d_type]
merged_list, merged_type = mergesort_multiple_lists(lists, types)

cd_up_matrix = array([[0]*8]*8)
cd_dn_matrix = array([[0]*8]*8)

cd_up_matrix[1][0] = 1
cd_up_matrix[3][2] = 1
cd_up_matrix[5][4] = 1
cd_up_matrix[7][6] = 1

cd_dn_matrix[2][0] = 1
cd_dn_matrix[3][1] = 1
cd_dn_matrix[6][4] = 1
cd_dn_matrix[7][5] = 1

c_up_matrix = cd_up_matrix.transpose()
c_dn_matrix = cd_dn_matrix.transpose()

cd_up_new = dot(ut, dot(cd_up_matrix, u))
cd_dn_new = dot(ut, dot(cd_dn_matrix, u))
c_up_new = dot(ut, dot(c_up_matrix, u))
c_dn_new = dot(ut, dot(c_dn_matrix, u))

operators_eig = [cd_up_new, c_up_new, cd_dn_new, c_dn_new]

local_term = dot(expv(-1.0 * v * merged_list[0]), ut)
local_term = dot(operators_eig[merged_type[0]], local_term)

n = len(merged_list)
for i in range(1, n):
    t = merged_list[i]-merged_list[i-1]
    local_term = dot(expv(-1.0 * v * t), local_term)
    local_term = dot(operators_eig[merged_type[i]], local_term)

local_term = dot(expv(-1.0 * v * (beta - merged_list[n-1])), local_term)
local_term = dot(u, local_term)


operators_occ = [cd_up_matrix, c_up_matrix, cd_dn_matrix, c_dn_matrix]

local_brute = expm(-1.0 * m * merged_list[0])
local_brute = dot(operators_occ[merged_type[0]], local_brute)
for i in range(1, n):
    t = merged_list[i]-merged_list[i-1]
    local_brute = dot(expm(-1.0 * m * t), local_brute)
    local_brute = dot(operators_occ[merged_type[i]], local_brute)

local_brute = dot(expm(-1.0 * m * (beta - merged_list[n-1])), local_brute)

