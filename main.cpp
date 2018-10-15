//Propag-Bin.c
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>



double fon_d2g(double t, double fc)
//*** Seismic source function***
//Source: Ricker(Gaussian second derivative )
{
double amp;
double pi=3.141592653589793238460;
amp = - pi * (pi * fc * t) * (pi * fc * t);
amp = exp(amp);
amp *= 1.0 - 2.0 * pi * (pi * fc * t) * (pi * fc * t);
return (amp);
}

int main()
{
//Variables declaration
int Nx, Nz, Nt, Ns, ixf, izf, iprof, k_snap, k_sismo, dt_sismo, dt_snap, cerjan_n, i, j, k;
float v;
double Dt, h, fcorte, pi, nf, tf, fc, cput, cerjan_fat, *g, **vel, **sismo, **p1, **p2, **p3;
char *modelo, *sismo_arq;
FILE *parametros; //data input file name
//FILE *locais; //file that indicates where to save the files created
FILE *vel_arq; //velocity files name
FILE *varsnap; //snapshots file name
FILE *seismogram; //seismogram file name
clock_t start, end;
//Data input
parametros = fopen ("Parametros.txt", "r");
fscanf (parametros, "%d %d %lf %d %d %d %lf %d", &Nx, &Nz, &h, &ixf, &izf, &iprof, &Dt, &Nt);
//Nx x Nz = grid size
//h = receivers spacing
//ixf x izf = seismic source position
//iprof = seismogram reading depth
//Dt = delta_t
//Nt = time increments
//modelo = velocity file name
//sismo_arq = seismogram file name
pi=3.141592653589793238460; //Pi
fcorte = 40.0; //cutting frequency
cerjan_n = 100;
cerjan_fat = 0.00105; //AB Cerjan
//Number of time steps for**********************:
dt_snap = 400; //saves snapshots every 400 iterations
dt_sismo = 10; //Writes one seismogram sample every 10 iterations
Ns = Nt / dt_sismo; //Number of Seismogram samples
g = (double *) malloc (cerjan_n * sizeof(double));
vel = (double **) malloc (Nz * sizeof(double *));
for (i = 0; i < Nz; i++)
vel[i] = (double *) malloc (Nx * sizeof(double));
sismo = (double **) malloc (Ns * sizeof(double *));
for (i = 0; i < Ns; i++)
sismo[i] = (double *) malloc (Nx * sizeof(double));
p1 = (double **) malloc (Nz * sizeof(double *));
for (i = 0; i < Nz; i++)
p1[i] = (double *) malloc (Nx * sizeof(double));
p2 = (double **) malloc (Nz * sizeof(double *));
for (i = 0; i < Nz; i++)
    p2[i] = (double *) malloc (Nx * sizeof(double));
p3 = (double **) malloc (Nz * sizeof(double *));
for (i = 0; i < Nz; i++)
p3[i] = (double *) malloc (Nx * sizeof(double));
modelo = (char *) malloc (256 * sizeof(char));
sismo_arq = (char *) malloc (256 * sizeof(char));
modelo = "vp_model_400x300.bin";
sismo_arq = "Seismogram.bin";
fclose(parametros);
//p1(Nx,Nz) corresponds to time k-1
//p2(Nx,Nz) corresponds to time k
//p3(Nx,Nz) corresponds to time k+1

//Data Processing

//Processing time
start = clock();
//Source term calculation
nf = 4 * sqrt(pi) / (fcorte * Dt);
tf = 2 * sqrt(pi) / fcorte;
fc = fcorte / (3.0 * sqrt(pi));
//Dumping factors calculation
for (i = 0; i < cerjan_n; i++)
g[i] = exp( - pow(cerjan_fat * (cerjan_n - i), 2) );
printf("Source application length: Nf = %.4f\n", nf);
printf("Seismogram temporal samples number: Ns = %d\n", Ns);
//reading velocity grid
printf("\nSaving velocity model...\n");
vel_arq = fopen(modelo, "rb");
for (j = 0; j < Nx; j++)
for (i = 0; i < Nz; i++)
{
fread (&v, sizeof(float), 1, vel_arq);
vel[i][j] = (double) v;
}
fclose(vel_arq);
//Wavefield initialization (initial condition)
for (i = 0; i < Nz; i++)
for (j = 0; j < Nx; j++)
{
p1[i][j] = 0;
p2[i][j] = 0;
}
k_snap = 0; //Snapshots counter
k_sismo = 0; //Seismogram sample counter
printf("\nStart\n");
for (k = 0; k < Nt; k++)//Start time loop
{
//Printing the time loop on the screen
if (k % 50 == 0)
printf("n = %d\n", k);
//Applying the source
if ((double) k <= nf)
p1[izf][ixf] -= fon_d2g((float) (k - 1) * Dt - tf, fc);
//FDM operator (involving the whole grid, except boundaries)
//Fourth order spatial operator
for (i = 2; i < Nz - 2; i++)
    for (j = 2; j < Nx - 2; j++)
p3[i][j] = pow(vel[i][j] * (Dt/12*h), 2)  * (-(p2[i-2][j] + p2[i][j-2] + p2[i+2][j] + p2[i][j+2]) + 16 * (p2[i-1][j] + p2[i][j-1] + p2[i+1][j] + p2[i][j+1]) - 60 * p2[i][j]) + 2 * p2[i][j] - p1[i][j];
//Second order spatial operator
//p3[i][j] = pow(vel[i][j] * (Dt/h), 2)  * ( p2[i+1][j] - 2*p2[i][j] + p2[i-1][j] + p2[i][j+1] - 2*p2[i][j] + p2[i][j-1]) + 2*p2[i][j] - p1[i][j];

//Superior and inferior
for (j = 2; j < Nx - 2; j++)
{
p3[1][j] = (pow(vel[1][j] * Dt / h, 2)) * (p2[2][j] + p2[1][j+1] + p2[0][j] + p2[1][j-1] - 4 * p2[1][j]) + 2 * p2[1][j] - p1[1][j];
p3[Nz-2][j] = (pow(vel[Nz-2][j] * Dt / h, 2)) * (p2[Nz-1][j] + p2[Nz-2][j+1] + p2[Nz-3][j] + p2[Nz-2][j-1] - 4 * p2[Nz-2][j]) + 2 * p2[Nz-2][j] - p1[Nz-2][j];
}
//Left and Right
for (i = 1; i < Nz - 1; i++)
{
p3[i][1] = (pow(vel[i][1] * Dt / h, 2)) * (p2[i+1][1] + p2[i][2] + p2[i-1][1] + p2[i][0] - 4 * p2[i][1]) + 2 * p2[i][1] - p1[i][1];
p3[i][Nx-2] = (pow(vel[i][Nx-2] * Dt / h, 2)) * (p2[i+1][Nx-2] + p2[i][Nx-1] + p2[i-1][Nx-2] + p2[i][Nx-3] - 4 * p2[i][Nx-2]) + 2 * p2[i][Nx-2] - p1[i][Nx-2];
}
//Applying non-reflexive conditions on the boundaries
//Superior
for (j = 1; j < Nx - 1; j++)
p3[0][j] = 0; //Free surface condition
//Left
for (i = 0; i < Nz - 1; i++)
p3[i][0] = p2[i][0] + (Dt * vel[i][0] / h) * (p2[i][1] - p2[i][0]);
//Right
for (i = 0; i < Nz - 1; i++)
p3[i][Nx-1] = p2[i][Nx-1] + (Dt * vel[i][Nx-1] / h) * (p2[i][Nx-2] - p2[i][Nx-1]);
//Inferior
for (j = 0; j < Nx; j++)
p3[Nz-1][j] = p2[Nz-1][j] + (Dt * vel[Nz-1][j] / h) * (p2[Nz-2][j] - p2[Nz-1][j]);
//Applying dumping zone close to the boundaries
//Left
for (i = 0; i < Nz; i++)
for (j = 0; j < cerjan_n; j++)
{
p3[i][j] *= g[j];
p2[i][j] *= g[j];
}
//Right
for (i = 0; i < Nz; i++)
    for (j = 0; j < cerjan_n; j++)
{
p3[i][Nx - j - 1] *= g[j];
p2[i][Nx - j - 1] *= g[j];
}
//Inferior
for (i = 0; i < cerjan_n; i++)
for (j = 0; j < Nx; j++)
{
p3[Nz - i - 1][j] *= g[i];
p2[Nz - i - 1][j] *= g[i];
}
//Saves a seismogram sample every "dt_sismo" time steps

if (k % dt_sismo == 0)
{
for (j = 0; j < Nx; j++)
sismo[k_sismo][j] = p3[iprof][j];
k_sismo++;
}

//Snapshots print every "dt_snap" time intervals
if (k % dt_snap == 0 && k > 0)
{
char *nome;
nome = (char *) malloc (strlen(sismo_arq) * sizeof(char));
for (i = 0; i < strlen(sismo_arq)-14; i++)
nome[i] = sismo_arq[i];
nome[i] = 'S';
nome[i+1] = 'n';
nome[i+2] = 'a';
nome[i+3] = 'p';
nome[i+4] = '_';
nome[i+5] = '\0';
char snap[2];
itoa(k / dt_snap, snap, 10);
strcat(nome, snap);
strcat(nome, ".bin");
varsnap = fopen(nome, "wb");
printf("Snapshot = %s\n", nome);
for (j = 0; j < Nx; j++)
for (i = 0; i < Nz; i++)
{
v = (float) p3[i][j];
fwrite(&v, sizeof(float), 1, varsnap);
}
fclose(varsnap);
free (nome);
k_snap++;
}

//Wavefield update for the next time loop
for (i = 0; i < Nz; i++)
for (j = 0; j < Nx; j++)
p1[i][j] = p2[i][j];
for (i = 0; i < Nz; i++)
for (j = 0; j < Nx; j++)
p2[i][j] = p3[i][j];
}
//Data output
//Calculate processing time and prints it
end = clock();
cput = ((double) (end - start)) / CLOCKS_PER_SEC;
printf("\nExecution time: %.4f s.", cput);
seismogram = fopen(sismo_arq, "wb");
for(j = 0; j < Nx; j++)
for(i = 0; i < Ns; i++)
{
v = (float) sismo[i][j];
fwrite (&v, sizeof(float), 1, seismogram);
}
fclose(seismogram);
free(g);
free(vel);
free(sismo);
free(p1);
free(p2);
free(p3);
free(modelo);
free(sismo_arq);
return 0;
}
