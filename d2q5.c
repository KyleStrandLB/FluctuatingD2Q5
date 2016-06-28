/****************************************************************/
/*        D2Q5 Fluctuating Lattice Boltzmann Simulation         */
/*                      Diffusive System                        */
/*            Kyle Strand:  kyle.t.strand@ndsu.edu              */
/*               North Dakota State University                  */
/*                     24 June 2016                          */
/****************************************************************/

/* New lattice Boltzmann algorithm for fluctuating diffusion.

   Requires Alexander Wagner's Graph Library.
	https://www.ndsu.edu/pubweb/~carswagn/GUI/index.html

   The algorithm will follow the following methodology:
	1) Forward matrix transformation
	2) Collision step - Noise added here
	3) Backward matrix transformation
	4) Streaming step 
*/ 

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mygraph.h>
#include <unistd.h>
#include <time.h>
#include <complex.h>

#define xdim 100                                   // Number of x lattice points
#define ydim 100                                   // Number of y lattice points

double f[5][xdim][ydim], n[xdim][ydim]; 
double noise[5], tau[5]={1,0.6,0.6,0.6,0.6};
double n0[2]={30,30}, T0 = 1./3.;
int next = 0, Pause = 1, done = 0, repeat = 1, iterations, writedata = 0, interval=1, intervalcnt=0;
char filename[500]="data.dat";

void init() {                                       // Initializing Eq. Dists

  iterations = 0;
  for (int i = 1; i < xdim-1; i++) {
    for (int j = 1; j < ydim-1; j++) {                  
      n[i][j]=n0[0];
      f[0][i][j] = n[i][j] * (1 - 2*T0);
      for (int a=1; a<5; a++) f[a][i][j]=n[i][j]/2. * T0;
    }
  }
}

void iteration() {                                  // Iteration step
  
  double M[5];                              

  double m[5][5] = {
                     {1.,1.,1.,1.,1.},
                     {0.,sqrt(1./T0),-sqrt(1./T0),0.,0.},
                     {0,0,0,sqrt(1./T0),-sqrt(1./T0)},
                     {0,sqrt(1./(2.*T0)),sqrt(1./(2.*T0)),-sqrt(1./(2.*T0)),-sqrt(1./(2.*T0))},
                     {-sqrt(2.*T0/(1.-2.*T0)),sqrt((1.-2.*T0)/(2.*T0)),sqrt((1.-2.*T0)/(2.*T0)),sqrt((1.-2.*T0)/(2.*T0)),sqrt((1.-2.*T0)/(2.*T0))}
                   };                                   // Transformation matrix

  for (int i = 1; i < xdim-1; i++) {                                     
      for (int j = 1; j < ydim-1; j++) {
        n[i][j] = f[0][i][j] + f[1][i][j] + f[2][i][j] + f[3][i][j] + f[4][i][j];                   // Sum of Eq dists = density

        for (int a=0; a<5; a++) M[a]=0;

        for (int a=0; a<5; a++) {
          for (int b=0; b<5; b++) {         //Forward transformation
            M[a] += m[a][b] * f[b][i][j];   // the matrix multiplication was backwards before i think
          }
        }

      noise[0]=0;
      for (int a=1; a<5; a++) noise[a]=sqrt(n[i][j]*12)*((double)rand()/RAND_MAX - 0.5);     //local density noise

      M[1] = ((1.-1./tau[1]) * M[1] + 1./tau[1]*(sqrt(2 * tau[1] - 1.) * noise[1]));  
      M[2] = ((1.-1./tau[2]) * M[2] + 1./tau[2]*(sqrt(2 * tau[2] - 1.) * noise[2]));
      M[3] = ((1.-1./tau[3]) * M[3] + 1./tau[3]*(sqrt(2 * tau[3] - 1.) * noise[3]));
      M[4] = ((1.-1./tau[4]) * M[4] + 1./tau[4]*(sqrt(2 * tau[4] - 1.) * noise[4]));

      f[0][i][j] = f[1][i][j] = f[2][i][j] = f[3][i][j] = f[4][i][j] = 0;

      for (int a=0; a<5; a++) {          // Back transform
          for (int u=0; u<5; u++) {
            if (a == 0) f[a][i][j] += m[u][a] *(1-2*T0)*M[u];                      //m[b][a] * M[b]
            else f[a][i][j] += m[u][a] * (T0/2) *M[u];
          }
      }
    }
  } 

    for (int y=0;y<ydim;y++) {      // periodic boundary conditions
      f[1][0][y]=f[1][xdim-2][y];
      f[2][xdim-1][y]=f[2][1][y];
    }
    for (int x=0;x<xdim;x++) {
      f[3][x][0]=f[3][x][ydim-2];
      f[4][x][ydim-1]=f[4][x][1];
    }

      memmove(&f[1][1][0],&f[1][0][0],(xdim-1)*ydim*sizeof(double));       //Check this
      memmove(&f[2][0][0],&f[2][1][0],(xdim-1)*ydim*sizeof(double));
      memmove(&f[3][0][1],&f[3][0][0],(xdim*ydim-1)*sizeof(double));
      memmove(&f[4][0][0],&f[4][0][1],(xdim*ydim-1)*sizeof(double));        // i had paraenths on xdim*(ydim-1)

  iterations++;
}

void GUI() {
  static int Xdim = xdim;
  static int Ydim = ydim;                                  //

  DefineGraphNxN_R("f0",&f[0][0][0], &Xdim, &Ydim, NULL);
  DefineGraphNxN_R("f1",&f[1][0][0], &Xdim, &Ydim, NULL);
  DefineGraphNxN_R("f2",&f[2][0][0], &Xdim, &Ydim, NULL);
  DefineGraphNxN_R("f3",&f[3][0][0], &Xdim, &Ydim, NULL);
  DefineGraphNxN_R("f4",&f[4][0][0], &Xdim, &Ydim, NULL);
  DefineGraphNxN_R("n",&n[0][0], &Xdim, &Ydim, NULL);
  NewGraph();
  StartMenu("Fluctuating D2Q5",1);
  DefineInt("Iterations",&iterations);
  StartMenu("Parameters",0);
    DefineDouble("n0_x",&n0[0]);
    DefineDouble("n0_y",&n0[1]);
    DefineDouble("T0",&T0);
    StartMenu("Relaxation Times",0);
      DefineDouble("Tau_0",&tau[0]);
      DefineDouble("Tau_1",&tau[1]);
      DefineDouble("Tau_2",&tau[2]);
      DefineDouble("Tau_3",&tau[3]);
      DefineDouble("Tau_4",&tau[4]);
    EndMenu();    
  EndMenu();
  DefineFunction("init",&init);
  SetActiveGraph(0);
  DefineGraph(contour2d_,"Graphs");
  DefineInt("Repeat",&repeat);
  DefineBool("Next",&next);
  DefineBool("Pause",&Pause);
  DefineBool("Close",&done);
  EndMenu();
}

int main(int argc, char *argv[]) {
  int newdata = 1;
  int i;

  init();
  GUI();

  while (done == 0) {
    Events(newdata);
    //GetData();
    DrawGraphs();
    if (next || !Pause) {
      newdata = 1;
      next = 0;
      for (i = 0; i < repeat; i++) {
        iteration();
      }
    }
    else sleep(1);
  }

  return 0;
}



