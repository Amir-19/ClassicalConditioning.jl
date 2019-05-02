using LinearAlgebra
using Statistics
X_prime = [3,3,3]
X = copy((X_prime - minimum(X_prime))/(maximum(X_prime)-minimum(X_prime)))
println(X)
