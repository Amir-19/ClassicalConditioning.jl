
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
using Plots

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

function figure_test()
    background = zeros(2,1)
    background[1] = 1.0
    CS_and_background = ones(2,1)
    ep = Experiment(2,0.1,1.0,0.95,0.2,0,0,zeros(2,1),zeros(2,1))
    for i = 1:20
        steps(4,CS_and_background,0.0,ep)
        steps(10,background,0.0,ep)
        steps(1,background,1.0,ep)
        steps(100,background,0.0,ep)
        print("v_back: ");print(ep.V[1]);print(" v_cs: ");print(ep.V[2]);
        println(" ")
    end
end

function figure_19()
    background = zeros(3,1)
    background[1] = 1.0
    CSp_and_background = ones(3,1)
    CSp_and_background[3] = 0.0
    CSn_and_background = ones(3,1)
    CSn_and_background[2] = 0.0
    CSp_and_CSn_and_background = ones(3,1)

    ep = Experiment(3,0.1,1.0,0.95,0.2,0,0,zeros(3,1),zeros(3,1))
    plot_x_time_steps = []
    plot_y_V_csp = []
    plot_y_V_csn = []

    for i = 1:80
        steps(100,background,0.0,ep)
        steps(4,CSp_and_background,0.0,ep)
        steps(2,background,1.0,ep)
        steps(100,background,0.0,ep)
        steps(4,CSp_and_CSn_and_background,0.0,ep)
        steps(2,background,0.0,ep)

        append!(plot_x_time_steps,i)
        append!(plot_y_V_csp,ep.V[2])
        append!(plot_y_V_csn,ep.V[3])
    end
    for i = 81:130
        steps(100,background,0.0,ep)
        steps(4,CSp_and_background,0.0,ep)
        steps(2,background,0.0,ep)
        steps(100,background,0.0,ep)
        steps(4,CSn_and_background,0.0,ep)
        steps(2,background,0.0,ep)

        append!(plot_x_time_steps,i)
        append!(plot_y_V_csp,ep.V[2])
        append!(plot_y_V_csn,ep.V[3])
    end
    plot(plot_x_time_steps,[plot_y_V_csp,plot_y_V_csn],label=["CS+" "CS-"])
end

function figure_20()
    background = zeros(3,1)
    background[1] = 1.0
    CSA_and_background = ones(3,1)
    CSA_and_background[3] = 0.0
    CSB_and_background = ones(3,1)
    CSB_and_background[2] = 0.0
    CSA_and_CSB_and_background = ones(3,1)

    ep_present = Experiment(3,0.1,1.0,0.95,0.2,0,0,zeros(3,1),zeros(3,1))
    ep_absent = Experiment(3,0.1,1.0,0.95,0.2,0,0,zeros(3,1),zeros(3,1))
    plot_x_time_steps = []
    plot_y_V_csa_bpresent = []
    plot_y_V_csa_babsent = []

    for i = 1:80
        steps(100,background,0.0,ep_absent)
        steps(4,CSA_and_background,0.0,ep_absent)
        steps(4,background,0.0,ep_absent)
        steps(2,background,1.0,ep_absent)

        append!(plot_x_time_steps,i)
        append!(plot_y_V_csa_babsent,ep_absent.V[2])
    end

    for i = 1:80
        steps(100,background,0.0,ep_present)
        steps(4,CSA_and_background,0.0,ep_present)
        steps(4,CSB_and_background,0.0,ep_present)
        steps(2,background,1.0,ep_present)

        append!(plot_y_V_csa_bpresent,ep_present.V[2])
    end
    plot(plot_x_time_steps,[plot_y_V_csa_babsent,plot_y_V_csa_bpresent],
        label=["CSB Absent" "CSB Present"])
end

function figure_21()
    background = zeros(3,1)
    background[1] = 1.0
    CSA_and_background = ones(3,1)
    CSA_and_background[3] = 0.0
    CSB_and_background = ones(3,1)
    CSB_and_background[2] = 0.0
    CSA_and_CSB_and_background = ones(3,1)

    ep_present = Experiment(3,0.1,1.0,0.95,0.2,0,0,zeros(3,1),zeros(3,1))
    ep_absent = Experiment(3,0.1,1.0,0.95,0.2,0,0,zeros(3,1),zeros(3,1))
    plot_x_time_steps = []
    plot_y_V_csb_apresent = []
    plot_y_V_csb_aabsent = []

    for i = 1:80
        steps(100,background,0.0,ep_absent)
        steps(4,background,0.0,ep_absent)
        steps(4,CSB_and_background,0.0,ep_absent)
        steps(2,background,1.0,ep_absent)

        append!(plot_x_time_steps,i)
        append!(plot_y_V_csb_aabsent,ep_absent.V[3])
    end

    for i = 1:80
        steps(100,background,0.0,ep_present)
        steps(4,CSA_and_background,0.0,ep_present)
        steps(4,CSA_and_CSB_and_background,0.0,ep_present)
        steps(2,background,1.0,ep_present)

        append!(plot_y_V_csb_apresent,ep_present.V[3])
    end
    plot(plot_x_time_steps,[plot_y_V_csb_aabsent,plot_y_V_csb_apresent],
        label=["CSA Absent" "CSA Present"])
end

figure_test()
