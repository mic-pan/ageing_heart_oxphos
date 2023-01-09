#=
Author: Giovanni Guglielmi and Michael Pan
Description: sim2PowerGaussian generates realisations of 
             2^X, where X~MV(m, Cov). 
             Since the algorithm is based on pseudo-random generation of 
             numbers, please set the seed before launching the function for
             the sake of reproducibility. 
=#
module DataDistributions

using DataFrames, CSV, LinearAlgebra, Distributions
export read_data, filter_names, split, sim2Gaussian

"""
Reads in data from a CSV file.
"""
function read_data(filename)
    data = CSV.read(filename,DataFrame)
    rename!(data,:Column1 => :name)
    return data
end

"""
Splits dataset into two groups. The idx argument is the number of individuals in group 1.
"""
function split(data,idx)
    col1 = data[:,[1]]
    group1 = data[:,2:idx+1]
    group2 = data[:,idx+2:end]

    df1 = [col1 group1]
    df2 = [col1 group2]

    return (df1,df2)
end

"""
Filters the data by entries with name in the list given by the argument names.
"""
filter_names(data,names) = subset(data, :name => ByRow(n -> n in names))


"""
'''
Description:
-----------
    Generate nSim (default: 10000) of random variable X, 
    where X ~ Gaussian(mu_vec, sigma_mat). 
    The estimation of mu_vec and sigma_mat is done through MLE procedure,
    using the analytical formula

Input:
-----
    log2Exprs <- N*K real-valued matrix (or DataFrame). 
                 N = observation (samples), K = variables (analytes) 
                 This matrix should be in log2 scale.
    nSim <- real-valued scalar. 
            It gives us the number of patients.
    covBool <- Boolean
               Will remove correlation structure if set to false
    power <- Boolean
             Will output 2^X instead of X if set to true
    zero_mean <- Boolean
                 Will set mean to zero if set to true
               
Output:
------
    <- real-valued matrix nSim*K. 
       nSim simulated samples, K analytes with estimated average, and 
       covariance structure. This is 2^X random variable
"""
function sim2Gaussian(log2Exprs, nSim=10000; covBool=true, power=false, zero_mean=true)
    matrix = Matrix(log2Exprs[:,2:end])
    n_names = size(matrix,1)

    if zero_mean
        μ_mle = zeros(n_names)
    else
        μ_mle = reshape(mean(matrix, dims=2),n_names)
    end

    σ_mle = cov(matrix, dims=2, corrected=false)
    if covBool
        σ_mle = Diagonal(σ_mle)
    end

    dist = MvNormal(μ_mle,σ_mle)
    x = rand(dist,nSim)
    
    if power
        return 2 .^ x
    else
        return x
    end
end

end