#=
    Julia code for the Temporal-Difference (TD) model of classical conditioning.
    Including state augmentation to fill the trace interval gap.

    This code was written by Amir Samani, inspired by TD model of classical
    conditioning code written by Rich Sutton.  April 12, 2019
=#
using LinearAlgebra
using Statistics
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
    #return value
end

function calculate_traces(X,ep::ExperimentSettings,is_onset::Bool, src_vector::Array)
    X_prime = zeros(size(X))
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
    return X_prime
end

function save_data(X,Z,λ,t,episode,ep_data::ExperimentData)
    for i=1:size(X)[1]
        ep_data.feature[episode,t+1,i]=X[i]
        ep_data.Z[episode,t+1,i]=Z[i]
    end
    ep_data.US[episode,t+1,1]=λ
end

function steps(num_steps, X, λ, ep::ExperimentSettings, ep_data::ExperimentData, is_onset::Bool, src_vector::Array,episode_index)
    normalization = "sumtoone"
    Vbar_t = 0
    alpha_beta_error = 0
    X_prime = zeros(size(X))
    for i = 1:num_steps
        save_data(X,ep.Z,λ,ep.t,episode_index,ep_data)
        X_prime = calculate_traces(X,ep,is_onset,src_vector)
        Vbar_t = v_bar(ep.V, X)
        alpha_beta_error = ep.α * ep.β * (λ + ep.γ*Vbar_t - ep.Vbar_prev_t)
        ep.t += 1
        ep.V += alpha_beta_error * ep.Z
        ep.Z += ep.δ * (X - ep.Z)
        ep.Vbar_prev_t = v_bar(ep.V, X)
        if normalization == "scaling"
            X = copy((X_prime .- minimum(X_prime))/(maximum(X_prime)-minimum(X_prime)))
        elseif normalization == "normalization"
            X = copy((X_prime .- mean(X_prime))/(maximum(X_prime)-minimum(X_prime)))
        elseif normalization == "standardaztion"
            X = copy((X_prime .- minimum(X_prime))/(stdm(X_prime,mean(X_prime))))
        elseif normalization == "sumtoone"
            X = copy(X_prime/sum(X))
        else
            X = copy(X_prime)
        end
    end
    return X_prime
end

function experiment_test_traces()
    println("-------------------------------------------------------------------------")
    m = 15 # size of the feature vector [background,CSs,Traces]

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
    max_time_step = 300 # for storage
    ep_data = ExperimentData(CS=zeros(num_episodes,max_time_step,1),US=zeros(num_episodes,max_time_step,1),feature=zeros(num_episodes,max_time_step,m),Z=zeros(num_episodes,max_time_step,m),td_error = zeros(num_episodes,max_time_step,1), actual_predition=zeros(num_episodes,max_time_step,1))

    cap_time = 0

    for i = 1:num_episodes
        # presenting background
        feature_vector = steps(100,background,0.0,ep,ep_data,false,src_vector,i)
        feature_vector = copy(CS_and_background)
        # presenting CS and background
        feature_vector = steps(4,feature_vector,0.0,ep,ep_data,true,src_vector,i)
        feature_vector[2] = 0
        # trace Interval
        feature_vector = steps(10,feature_vector,0.0,ep,ep_data,false,src_vector,i)
        # presenting the US
        feature_vector = steps(10,feature_vector,1.0,ep,ep_data,false,src_vector,i)
        # inter-trail
        feature_vector = steps(100,feature_vector,0.0,ep,ep_data,false,src_vector,i)

        cap_time = ep.t
        ep.t = 0
    end
    # plotting code
    prediction_data = []
    CS_data = []
    US_data = []
    selected_feature_data = []
    #for i=1:cap_time
    for i=90:140
        prediction_data = [prediction_data;dot(ep.V[2:end]',ep_data.feature[num_episodes,i,2:end])]
        CS_data = [CS_data;ep_data.feature[num_episodes,i,2]]
        US_data = [US_data;ep_data.US[num_episodes,i]]
    end
    data = [CS_data,US_data,prediction_data]
    plot(data,layout = (3,1))
end
experiment_test_traces()
