
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
end
function v_bar(V, X)
    value = dot(V', X)
    return value >= 0 ? value : 0
end
function steps(num_steps, X, λ, exp_detail::Experiment)

end

function fig_tes()
    e = Experiment(2,0.1,1.0,0.95,0.2,0,[Stimulus(1,2), Stimulus(2,3)])
end
function figure_19()
    e = Experiment(2,0.1,1.0,0.95,0.2,0,[Stimulus(1,2), Stimulus(2,3)])
end
