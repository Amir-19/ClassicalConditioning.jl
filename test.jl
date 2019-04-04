
#=  Julia code for the Temporal-Difference (TD) model of classical conditioning.
    As specified in:

        Sutton, R.S., Barto, A.G. (1990) "Time-Derivative Models of Pavlovian
        Reinforcement," in Learning and Computational Neuroscience: Foundations
        of Adaptive Networks, M. Gabriel and J. Moore, Eds., pp. 497--537.
        MIT Press. http://incompleteideas.net/papers/sutton-barto-90.pdf

    This code was written by Amir Samani, inspired by the C++ program written by
    Rich Sutton.  April 4, 2019

=#
struct Stimulus
    onset::UInt32                   # onset time step of the stimulus
    offset::UInt32                  # offset time step of the stimulus
end

mutable struct Experiment
    num_stimuli::UInt32             # number of stimuli, including US
    α::Float16                      # TD model of Classical Conditioning params
    β::Float16
    γ::Float16
    δ::Float16
    t::UInt32                       # current time step in the experiment
    stimuli::Array                  # list of stimuli, with US as last index
end

e = Experiment(2,0.1,1.0,0.95,0.2,0,[Stimulus(1,2), Stimulus(2,3)])
e.num_stimuli = 234
print(e.num_stimuli)
