#=
    Julia code for the Temporal-Difference (TD) model of classical conditioning.
    Including state augmentation to fill the trace interval gap.

    This code was written by Amir Samani, inspired by TD model of classical
    conditioning code written by Rich Sutton.  April 12, 2019
=#
using LinearAlgebra
using Plots

num_traces = 6
from_trace = 2

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
    capture::UInt32
end
mutable struct GenSet
    num::UInt32
    from::UInt32
end

function v_bar(V, X)
    value = dot(V', X)
    return value >= 0 ? value : 0
end

function steps(num_steps, X, λ, ep::ExperimentSettings, ep_data::ExperimentData, is_onset::Bool)

    Vbar_t = 0
    alpha_beta_error = 0
    gen = true
    decay = 0.6

    for i = 1:num_steps
        if gen == true && num_traces != 0
            if is_onset == true
                X[from_trace + 1] = X[from_trace + 1]*decay + X[from_trace]
                is_onset = false
            elseif is_onset == false
                X[from_trace + 1] = X[from_trace + 1]*decay + 0
            end

            for i = num_traces:-1:2
                X[from_trace + i] = X[from_trace + i]*decay + X[from_trace + i - 1]
            end
        end
        ep_data.CS[ep.t+1] = X[2]
        ep_data.US[ep.t+1] = λ
        ep_data.feature[ep.t+1,:] = copy(X)
        ep_data.Z[ep.t+1] = ep.Z[ep_data.capture]

        Vbar_t = v_bar(ep.V, X)
        alpha_beta_error = ep.α * ep.β * (λ + ep.γ*Vbar_t - ep.Vbar_prev_t)
        ep_data.td_error[ep.t+1] = alpha_beta_error
        ep.V += alpha_beta_error * ep.Z
        ep.Z += ep.δ * (X - ep.Z)
        ep.Vbar_prev_t = v_bar(ep.V, X)
        xcopy = copy(X)
        xcopy[1] = 0
        ep_data.actual_predition[ep.t+1] = v_bar(ep.V,xcopy)
        ep.t += 1
    end
    return X
end

function first_experiment()
    println("------------------------------")
    background = zeros(2,1)
    background[1] = 1.0
    CS_and_background = ones(2,1)

    num_episodes = 80
    num_time_step = 200 # +1
    background = zeros(2+num_traces,1)
    background[1] = 1.0
    CS_and_background = ones(2+num_traces,1)
    representation = zeros(2+num_traces,1)

    ep = ExperimentSettings(2,0.1,1.0,0.95,0.2,0,0,zeros(2+num_traces,1)
                                                  ,zeros(2+num_traces,1))
    ep_data = ExperimentData(zeros(num_time_step,1),zeros(num_time_step,1)
                            ,zeros(num_time_step,num_traces+2),zeros(num_time_step,1)
                            ,zeros(num_time_step,1),zeros(num_time_step,1),2)
    plot_x_time_steps = []
    plot_y_V_background = []
    plot_y_V_cs = []
    cap_time = 0
    for i = 1:8000
        representation = copy(CS_and_background)
        representation = steps(1,representation,0.0,ep,ep_data,true) # present CS with background
        representation = steps(3,representation,0.0,ep,ep_data,false)
        representation[2] = 0
        representation = steps(6,representation,0.0,ep,ep_data,false)    # trace interval
        representation = steps(2,representation,1.0,ep,ep_data,false)        # US/reward
        representation = steps(10,representation,0.0,ep,ep_data,false)      # inter-trial-interval
        cap_time = ep.t
        ep.t = 0

        append!(plot_x_time_steps,i)
        append!(plot_y_V_cs,ep.V[2])
        append!(plot_y_V_background,ep.V[1])
    end
    #plot(plot_x_time_steps,[plot_y_V_background,plot_y_V_cs],label=["V background" "V cs"])
    data = [ep_data.CS[1:cap_time],ep_data.US[1:cap_time],ep_data.actual_predition[1:cap_time]]
    labelo = ["CS" "US" "prediction"]
    for i=from_trace+1:num_traces+from_trace
        data = [data, ep_data.feature[1:cap_time,i]]
        labelo = [labelo (i-from_trace)]
    end
    plot(data, layout = (3+num_traces,1), label = labelo)
end

first_experiment()
