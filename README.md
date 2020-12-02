ClassicalConditioning.jl
For more information about TD model of classical conditioning read [sutton-barto-90](http://incompleteideas.net/papers/sutton-barto-90.pdf)
The code for TD model of classical coditioning. `TDmodel.jl`
Examples of runs:

CS A presentfrom time step 0 to 2, CS B present from time step 2 to 4 and US occurs at time step 4 with the duration of 1 time step.

![alt text](https://github.com/Amir-19/ClassicalConditioning.jl/blob/master/Fig20.JPG)

CS A presentfrom time step 0 to 4, CS B present from time step 2 to 4 and US occurs at time step 4 with the duration of 1 time step

![alt text](https://github.com/Amir-19/ClassicalConditioning.jl/blob/master/Fig21.JPG)


The code for deep traces to be added to TD model. `state_construction.jl`
Using deep traces we can make features for the trace interval gap. Starting from presense representation and makeing traces of available features (including traces) to fill the gap of trace interval. This is an example for trace inteval of 10 time steps:

![alt text](https://github.com/Amir-19/ClassicalConditioning.jl/blob/master/deep_traces.JPG)
