# This script finds the largest number of processors and nodes
# you can use, based on the number of grid points in the i/j directions
# on your domain. 

# Note: The largest number may not decompose the best way. If you want
# additional values, set some print statements in the code below


# enter the namelist values of e_we and e_sn
e_we = 3551
e_sn = 2441
 
# number of cores you want to use per node (Cheyenne has a max of 36/node)
cores = 128

# The value for 'cores' gets incremented later, so we want a static variable for the original value 
cores_orig = cores

# set upper limit of nodes - the max you want to loop through
node_max = 200 

# This is the least number of grid points allowed for each processor. 
# Dont' change this value.
smallest_size = 10

x = 1
while x <= node_max:
       
# finds the factor pairs for the total number of cores
    def f(cores):
        factors = []
        for i in range(1, int(cores**0.5)+1):
            if cores % i == 0:
                factors.append((i, cores/i ))
        return factors

    factors = f(cores)
    print(factors)

# Of the factor pairs, this finds the closest values (pair) in that array
    closest_factors = factors[-1]

# Of the set of closest values, assign the i and j values
    i_array_value = closest_factors[0]
    j_array_value = closest_factors[-1]

# Calculate how the domain will be decomposed
    e_we_decomp = int(e_we / i_array_value )
    e_sn_decomp = int(e_sn / j_array_value )

# Once the decomposition becomes smaller than the least number of grid points
# allowed for each processor, the loop will quit and display the max 
# number of processors and nodes you can use for your domain.
    if ((e_sn_decomp < smallest_size) or (e_we_decomp < smallest_size)):

# test to see if the max number of processors allowed is within the number for a single node 
        initial_factor_pair = factors[0]
        initial_factor = initial_factor_pair[-1]
        if initial_factor == cores_orig:

# start with value of cores_orig and decrease by 1 for each iteration
# until the value is allowed
           y = cores_orig
           while y >= 1:
                processors = y

# finds the factor pairs for the total number of processors
# still testing processor values for a single node
                def f(processors):
                    factors = []
                    for i in range(1, int(processors**0.5)+1):
                        if processors % i == 0:
                            factors.append((i, processors/i ))
                    return factors

                factors = f(processors)

# Of the factor pairs, this finds the closest values (pair) in that array
# still testing processor values for a single node
                closest_factors = factors[-1]

# Of the set of closest values, assign the i and j values
# still testing processor values for a single node
                i_array_value = closest_factors[0]
                j_array_value = closest_factors[-1]

# Calculate how the domain will be decomposed
# still testing processor values for a single node
                e_we_decomp = int(e_we / i_array_value )
                e_sn_decomp = int(e_sn / j_array_value )

# Once the decomposition becomes larger or equal to the least number of grid points
# allowed for each processor, the loop will quit and display the max 
# number of processors and nodes you can use for your domain.
                if ((e_sn_decomp >= smallest_size) and (e_we_decomp >= smallest_size)): 
                    max_procs = (i_array_value * j_array_value)
                    print( "max # of processors that can be used is: ", max_procs)
                    print( "max # of nodes that can be used is 1 ")
                    break

# if you haven't reached your limit, the loop continues
# still testing processor values for a single node   
                else:
                    y -= 1

# if the size of the domain allows multiple nodes
        else:
            max_procs = (i_array_value * j_array_value) - cores_orig
            max_nodes = (max_procs / cores_orig)
            print( "max # of processors that can be used is: ", max_procs)
            print( "max # of nodes that can be used is: ", max_nodes)
        break

# If you haven't reached your limit, the loop continues    
    x += 1
    cores = (cores+cores_orig)


