#=
    Julia code for the Temporal-Difference (TD) model of classical conditioning.
    Including state augmentation to fill the trace interval gap.

    This code was written by Amir Samani, inspired by TD model of classical
    conditioning code written by Rich Sutton.  April 12, 2019
=#
using LinearAlgebra
using Plots

mutable struct ExperimentSettings
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

mutable struct ExperimentData
    CS::Array                       # [time_step+1][episode][id1]
    US::Array                       # [time_step+1][episode][id2]
    feature::Array                  # [time_step+1][episode][id3]
    Z::Array                        # [time_step+1][episode][id3]
    td_error::Array                 # [time_step+1][episode][id3]
    actual_predition::Array         # [time_step+1][episode][id3]
end

function v_bar(V, X)
    value = dot(V', X)
    return value >= 0 ? value : 0
end

function steps(num_steps, X, λ, ep::ExperimentSettings)

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

function first_experiment()
    background = zeros(2,1)
    background[1] = 1.0
    CS_and_background = ones(2,1)

    ep = ExperimentSettings(2,0.1,1.0,0.95,0.2,0,0,zeros(2,1),zeros(2,1))

    plot_x_time_steps = []
    plot_y_V_background = []
    plot_y_V_cs = []

    num_episodes = 80
    num_time_steps_per_episode = 107

    
    for i = 1:num_episodes
        steps(4,CS_and_background,0.0,ep) # present CS with background
        steps(1,background,0.0,ep)        # trace interval
        steps(2,background,1.0,ep)        # US/reward
        steps(100,background,0.0,ep)      # inter-trial-interval

        ep.t = 0

        append!(plot_x_time_steps,i)
        append!(plot_y_V_cs,ep.V[2])
        append!(plot_y_V_background,ep.V[1])
    end
    plot(plot_x_time_steps,[plot_y_V_background,plot_y_V_cs],
        label=["V background" "V cs"])
end

first_experiment()
