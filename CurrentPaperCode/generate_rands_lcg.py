##################################
### Helper function for        ###
### generating random sequences###
### using Linear Congruence    ###
### implemented in Python      ###
##################################

# Generate sequence with same parameters
def generate_lcg(seed, n, a, c, m):
    s = seed
    rand_seq = []
    for _ in range(n):
        s = (a * s + c) % m
        #rand_seq.append(s / m)  ## If normalizing to [0,1)
        rand_seq.append(s)
    return rand_seq

# #### Example usage
# # LCG parameters (must match R)
# a = 1664525
# c = 1013904223
# m = 2**32
# seed = 1234
# n = 16383  # same as in R
# 
# sequence_py = lcg(seed, n, a, c, m)
