/* C++ code for the Temporal-Difference (TD) model of classical conditioning.
   As specified in:

        Sutton, R.S., Barto, A.G. (1990) "Time-Derivative Models of Pavlovian 
        Reinforcement," in Learning and Computational Neuroscience: Foundations 
        of Adaptive Networks, M. Gabriel and J. Moore, Eds., pp. 497--537. 
        MIT Press.  ftp://ftp.cs.umass.edu/pub/anw/pub/sutton/sutton-barto-90.ps

   This code was written by Rich Sutton, a poor C programmer, based on a
   Lisp program for the same purposes.  June 3, 1996

   New experiments can easily be created by adding additional stimulus vectors, etc.
*/

#include <string.h>
#include <stdio.h>
#include <stream.h>

#define N 2                               // number of stimuli, length of arrays
#define alpha 0.1                         // Basic parameters of model
#define beta 1.0                          // See Sutton & Barto, 1990
#define delta 0.2
#define gamma 0.95

#define ISI 2

double V[N];                              // array of associative strengths
double trace[N];                          // array of stimulus traces
double old_Vbar;                          // current V dotted with last X
int time;                                 // steps since beginning of experiment

double background[N] =          {1.0, 0.0};    // Stimulus vectors
double CS_and_background[N] =   {1.0, 1.0};    // Make more for other experiments


// setup for an experiment

void setup()
{
  time = 0;
  old_Vbar = 0;
  for (int i=0; i<N; i++)
    {
      V[i] = 0.0;
      trace[i]=0;
    }
}


// Computes prediction Vbar for inputs X given associative strengths V

double Vbar (double* V, double* X)
{
  double value = 0;

  for (int i=0; i<N; i++)
    value += V[i]*X[i];

  if (value < 0.0)
    value = 0.0;

  return(value);
}


// Runs TD model for num_steps time steps with CSs=X (a list) and US=lambda

void steps (int num_steps, double* X, double lambda)
{
  double new_Vbar, alpha_beta_error;

  for (int k=0; k<num_steps; k++)
    {
      new_Vbar = Vbar(V,X);
      alpha_beta_error = alpha * beta * (lambda + gamma*new_Vbar - old_Vbar);
      time++;
      
      for (int i=0; i<N; i++)
        {
          V[i] += alpha_beta_error * trace[i];
          trace[i] += delta * (X[i]-trace[i]);
//        cout << time << ": i:" << i << " trace:" << trace[i] << " V:" << V[i] 
//             << " a_b_error:" << alpha_beta_error << endl;
        }
      old_Vbar = Vbar(V,X);
    }
}

// Does trace conditioning experiment for ISIs > CS duration (in number of steps).
// CS duration is 4, US duration is 1

void trace_conditioning()
{
  if (ISI < 4)
    {
      cout << "This program won't work for ISIs < CS duration" << endl;
    }

  setup();
  steps(100,background,0.0);              // inter-trial-interval

  for (int trial=0; trial<20; trial++)
    {
      steps(4,CS_and_background,0.0);        // present CS with background
      steps(ISI-4,background,0.0);           // trace interval
      steps(1,background,1.0);               // US/reward
      steps(100,background,0.0);             // inter-trial-interval
      
      cout << "V_background: " << V[0] << "   V_CS: " << V[1] << endl;

    }
}



// Does backwards conditioning experiment for ISIs >= US duration (in number of steps).
// CS duration is 4. US duration is 1"

void backward_conditioning()
{
  if (ISI < 1)
    {
      cout << "This program won't work for ISIs < CS duration" << endl;
    }

  setup();
  steps(100,background,0.0);              // inter-trial-interval

  for (int trial=0; trial<20; trial++)
    {
      steps(1,background,1.0);               // US/reward
      steps(ISI-1,background,0.0);           // trace interval
      steps(4,CS_and_background,0.0);        // present CS with background
      steps(100,background,0.0);             // inter-trial-interval
      
      cout << "V_background: " << V[0] << "   V_CS: " << V[1] << endl;

    }
}

void main (void)
{
  backward_conditioning();
}

