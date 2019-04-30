#=
    Julia code for the Temporal-Difference (TD) model of classical conditioning.
    Including state augmentation to fill the trace interval gap.

    This code was written by Amir Samani, inspired by TD model of classical
    conditioning code written by Rich Sutton.  April 12, 2019
=#
using LinearAlgebra
using Plots


mutable struct ExperimentSettings
    num_stimuli::UInt32          # number of stimuli, including US
    α::Float16                   # TD model of Classical Conditioning params
    β::Float16
    γ::Float16
    δ::Float16
    trace_decay::Float16         # decay param for the trace features
    t::UInt32                    # current time step in the experiment
    Vbar_prev_t::Float64
    V::Array
    Z::Array                     # eligibility trace vector
    ExperimentSettings(;num_stimuli,α,β,γ,δ,trace_decay,t,Vbar_prev_t,V,Z) = new(num_stimuli,α,β,γ,δ,trace_decay,t,Vbar_prev_t,V,Z)
end

mutable struct ExperimentData
    CS::Array                    # [time_step+1][episode][id1]
    US::Array                    # [time_step+1][episode][id2]
    feature::Array               # [time_step+1][episode][id3]
    Z::Array                     # [time_step+1][episode][id3]
    td_error::Array              # [time_step+1][episode][id3]
    actual_predition::Array      # [time_step+1][episode][id3]
    ExperimentData(;CS,US,feature,Z,td_error,actual_predition) = new(CS,US,feature,Z,td_error,actual_predition)
end

function v_bar(V, X)
    value = dot(V', X)
    return value >= 0 ? value : 0
end

function steps(num_steps, X, λ, ep::ExperimentSettings, ep_data::ExperimentData, is_onset::Bool, src_vector::Array)

    Vbar_t = 0
    alpha_beta_error = 0
    X_prime = zeros(size(X))
    for i = 1:num_steps
        for j = 1:size(src_vector)[1]
            if src_vector[j]!=0
                if src_vector[src_vector[j]]!=0
                    X_prime[j] = ep.trace_decay * X[j] + X[src_vector[j]]
                else
                    if is_onset == true
                        X_prime[j] = ep.trace_decay * X[j] + X[src_vector[j]]
                        is_onset = false
                    else
                        X_prime[j] = ep.trace_decay * X[j]
                    end
                end
            else
                X_prime[j] = X[j]
            end
        end
        X = copy(X_prime)
        Vbar_t = v_bar(ep.V, X)
        alpha_beta_error = ep.α * ep.β * (λ + ep.γ*Vbar_t - ep.Vbar_prev_t)
        ep.t += 1
        ep.V += alpha_beta_error * ep.Z
        ep.Z += ep.δ * (X - ep.Z)
        ep.Vbar_prev_t = v_bar(ep.V, X)
    end
    return X_prime
end

function experiment_test_traces()
    println("-------------------------------------------------------------------------")
    m = 5 # size of the feature vector [background,CSs,Traces]

    src_vector = collect(0:m-1)
    src_vector[2] = 0

    ep = ExperimentSettings(num_stimuli=m,α=0.1,β=1.0,γ=0.95,δ=0.2, trace_decay = 0.1,t=0,Vbar_prev_t=0,V=zeros(m,1),Z=zeros(m,1))


    # forming the representation vectors
    background = zeros(m,1)
    background[1] = 1.0
    CS_and_background = zeros(m,1)
    CS_and_background[1] = 1.0
    CS_and_background[2] = 1.0
    feature_vector = zeros(m,1)

    num_episodes = 100
    max_time_step = 200 # for storage
    ep_data = ExperimentData(CS=zeros(num_episodes,max_time_step,1),US=zeros(num_episodes,max_time_step,1),feature=zeros(num_episodes,max_time_step,1),Z=zeros(num_episodes,max_time_step,1),td_error = zeros(num_episodes,max_time_step,1), actual_predition=zeros(num_episodes,max_time_step,1))

    cap_time = 0

    for i = 1:num_episodes
        feature_vector = copy(CS_and_background)
        feature_vector = steps(4,feature_vector,0.0,ep,ep_data,true,src_vector) #CS+BG
        feature_vector[2] = 0
        feature_vector = steps(0,feature_vector,0.0,ep,ep_data,true,src_vector) #Trace Interval
        feature_vector = steps(1,feature_vector,1.0,ep,ep_data,true,src_vector) #US
        feature_vector = steps(100,feature_vector,0.0,ep,ep_data,true,src_vector) #inter-trail
    end
    sum_value = 0
    for i=2:m
        sum_value+=ep.V[i]
    end
    print(sum_value)
end
experiment_test_traces()
