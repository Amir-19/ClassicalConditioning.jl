
#=
    Julia code for the Temporal-Difference (TD) model of classical conditioning.
    As specified in:

        Sutton, R.S., Barto, A.G. (1990) "Time-Derivative Models of Pavlovian
        Reinforcement," in Learning and Computational Neuroscience: Foundations
        of Adaptive Networks, M. Gabriel and J. Moore, Eds., pp. 497--537.
        MIT Press. http://incompleteideas.net/papers/sutton-barto-90.pdf

    This code was written by Amir Samani, based on the C++ program written by
    Rich Sutton for the same purposes.  April 4, 2019
=#

using LinearAlgebra

mutable struct Experiment
    num_stimuli::UInt32             # number of stimuli, including US
    α::Float16                      # TD model of Classical Conditioning params
    β::Float16
    γ::Float16
    δ::Float16
    t::UInt32                       # current time step in the experiment
    Vbar_prev_t::Float64
    V::Array
    Z::Array                        # Trace vector
end
function v_bar(V, X)
    value = dot(V', X)
    return value >= 0 ? value : 0
end
function steps(num_steps, X, λ, ep::Experiment)

    Vbar_t = 0
    alpha_beta_error = 0

    for i = 1:num_steps
        Vbar_t = v_bar(ep.V, X)
        alpha_beta_error = ep.α * ep.β * (λ + ep.γ*Vbar_t - ep.Vbar_prev_t)
        ep.t += 1
        ep.V += alpha_beta_error * ep.Z
        ep.Z += ep.δ * (X - ep.Z)
        ep.Vbar_prev_t = v_bar(ep.V, X)
    end
end

function fig_tes()
    background = zeros(2,1)
    background[1] = 1.0
    CS_and_background = ones(2,1)
    ep = Experiment(2,0.1,1.0,0.95,0.2,0,0,zeros(2,1),zeros(2,1))
    steps(100,background,0.0,ep)
    for i = 1:20
        steps(1,background,1.0,ep)
        steps(2-1,background,0.0,ep)
        steps(4,CS_and_background,0.0,ep)
        steps(100,background,0.0,ep)
        print("v_back: ");print(ep.V[1]);print(" v_cs: ");print(ep.V[2]);
        println(" ")
    end
end
function figure_19()
    e = Experiment(2,0.1,1.0,0.95,0.2,0,[Stimulus(1,2), Stimulus(2,3)])
end
fig_tes()
