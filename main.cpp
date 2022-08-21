<<<<<<< HEAD:src/main.cpp

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

//using namespace std;

//Domain size
const int NX = 40;
const int NY = 40;

//Initialize weights
const double w0 = 4.0/9.0;  // zero weight
const double ws = 1.0/9.0;  // adjacent weight
const double wd = 1.0/36.0; // diagonal weight

//initialize direction and weight vectors
const int ndir = 9;
const double wi[] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};
const int dirx[] = {0,1,0,-1, 0,1,-1,-1, 1};
const int diry[] = {0,0,1, 0,-1,1, 1,-1,-1};

/*

## Direction Scheme ##

    6 2 5
    3 0 1
    7 4 8

*/

const int opposite[] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; //For storing opposite directions used in bounce back

//const double nu = 1.0/6.0;
const double tau = 0.6;
const double nu = (tau - 0.5)/3;
double tauinv = 1/tau;
const double u_max = 0.02;


const double rho0 = 1.0;    //Initial density
const double mu = rho0*nu;

//Number of timesteps
const unsigned int NSTEPS = 20000;
const unsigned int NSAVE  =  50;


//For Poseuille's flow
const double dpdx = 1.0e-3;

int field_index(int x, int y, int d){
    
    return (NX*(NY*d+y)+x);

}

int scalar_index(int x, int y){
    
    return (NX*y+x);

}

void taylor_green(int t, double *rho, double *ux, double *uy){

    double kx = (2*3.14159)/NX;
    double ky = (2*3.14159)/NY;
    double td = 1/(nu*(kx*kx + ky*ky));

    for(int x = 0; x < NX; x++){

        for(int y = 0; y < NY; y++){
            
            double X = x+0.5;
            double Y = y+0.5;

            ux[scalar_index(x, y)] = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1*t/td);
            uy[scalar_index(x, y)] = -u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1*t/td);

            double P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2*kx*X) + (kx/ky)*cos(2.0*ky*Y))*exp(-2*t/td);

            rho[scalar_index(x, y)] = rho0 + 3*P;

        }
    }
}

void poiseuille_flowAnalytical(double *ux){

    double velConst = -(dpdx)/(2*mu);

    //printf("%lf", velConst);

    for(int x = 0; x < NX; x++){

        for(int y = 0; y < NY; y++){

            int r = y - NY + 1;
            ux[scalar_index(x, y)] = velConst*(y)*r;

        }

    }
}


//Initialization function used for Poseuille's flow
void init_distribution(double *f, double *r, double *u, double *v){
    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){

            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            double rho = r[scalar_index(x, y)];


            for(int i = 0; i<ndir; i++){

                double cidotu = dirx[i]*ux + diry[i]*uy;
                f[field_index(x,y,i)] = wi[i]*rho0;

            }
        }

    }
}


void init_wall_distribution(double *f_wall){

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<2; y++){
            for(int i = 0; i<ndir; i++){
                f_wall[field_index(x, y, i)] = wi[i]*rho0;   
            }
        }
    }
}


void init_Geo(double *geo){
    for(int y = 0; y<NY; y++){
        for(int x = 0; x<NX; x++){
            if(y == 0 || y == NY-1){
                //Boundary nodes are labelled as 1
                geo[scalar_index(x, y)] = 1;
            }

            else{
                //Fluid nodes are labelled as 0
                geo[scalar_index(x, y)] = 0;
            }

            //std::cout << geo[scalar_index(x,y)];
            
        }

        //std::cout << std::endl;

        
    }
}


void compute_rho_u(double *f, double *r, double *u, double *v)
{
    for(int y = 0; y < NY; ++y)
    {
        for(int x = 0; x < NX; ++x)
        {
            double rho = 0.0;
            double ux  = 0.0;
            double uy  = 0.0;
            
            for(unsigned int i = 0; i < ndir; ++i)
            {
                rho += f[field_index(x,y,i)];
                //std::cout << rho << std::endl;
                ux  += dirx[i]*f[field_index(x,y,i)];
                uy  += diry[i]*f[field_index(x,y,i)];

            }

            ux += dpdx/2;   //Add this term for poseuille flow (pressure difference)
            
            r[scalar_index(x,y)] = rho;
            u[scalar_index(x,y)] = ux/rho;
            v[scalar_index(x,y)] = uy/rho;

            //std::cout << ux << std :: endl;

        }
    }
}



void collison(double *f, double *r, double *u, double *v, double *source, double *feq, double *ft){

    double omtau = (1 - (1/tau)); 

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            double rho = r[scalar_index(x, y)];

            double u2 = ux*ux + uy*uy;

            for(int i = 0; i<ndir; i++){

                double t1 = ux*dirx[i] + uy*diry[i];
                double t2 = t1*t1;

                source[i] = (1 - 0.5/tau)*wi[i]*(3*(dirx[i] - ux) + 9*(t1)*dirx[i])*dpdx;   //Added for Poseuille flow (External force)

                feq[i] = wi[i]*rho*(1 + 3*t1 + 4.5*t2 - 1.5*u2);
                ft[field_index(x, y, i)] = f[field_index(x, y, i)]*omtau  + feq[i]*tauinv + source[i];

                
            }

            
        }
            
    }

}


void stream(double *f, double *ft){

    //Function streams periodically

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){
            for(int i = 0; i < ndir; i++){

                unsigned int xmd = (NX+x+dirx[i])%NX;
                unsigned int ymd = (NY+y+diry[i])%NY;
                
                f[field_index(xmd,ymd,i)] = ft[field_index(x,y,i)];
                
            }
        }
    }
}




void bounce_back(double *geo, double *f, double *ft, double *f_wall){

    for(int y = 0; y<NY; y++){
        for(int x = 0; x<NX; x++){

            if(geo[scalar_index(x, y)] == 0){

                for(int i = 0; i<ndir; i++){

                    int xmd = (NX + x + dirx[i])%NX;
                    int ymd =  (NY + y + diry[i])%NY;


                    if(geo[scalar_index(xmd, ymd)] == 1){

                        f[field_index(x, y, opposite[i])] = ft[field_index(x, y, i)];

                    }
                }
            }
        }
    }
}



void save_scalar(const char *name, double *scalar, unsigned int n)
{
    const size_t mem_size_scalar = sizeof(double) * NX * NY;

    // assume reasonably-sized file names
    char filename[128];
    char format[16];

    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS) + 1.0);

    // file name format is name0000nnn.bin
    sprintf(format, "%%s%%0%dd.txt", ndigits);
    sprintf(filename, format, name, n); 

    FILE *output = fopen(filename, "w");
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            fprintf(output, "%.10f\t", scalar[scalar_index(x, y)]);
        }

        fprintf(output, "\n");
    }
    fclose(output);
}





int main(){

    //Allocate memory

    double *f  = (double*) malloc(sizeof(double)*NX*NY*ndir);
    double *f_wall  = (double*) malloc(sizeof(double)*NX*2*ndir);
    double *geo = (double*) malloc(sizeof(double)*NX*NY);
    double *ux  = (double*) malloc(sizeof(double)*NX*NY);
    double *uy  = (double*) malloc(sizeof(double)*NX*NY);
    double *rho  = (double*) malloc(sizeof(double)*NX*NY);
    double *source  = (double*) malloc(sizeof(double)*ndir);
    double *feq  = (double*) malloc(sizeof(double)*ndir);
    double *ft  = (double*) malloc(sizeof(double)*NX*NY*ndir);
    double *ux_star  = (double*) malloc(sizeof(double)*NX*NY);
    //double *uy_star  = (double*) malloc(sizeof(double)*NX*NY);
    //double *rho_star  = (double*) malloc(sizeof(double)*NX*NY);

    std::cout << "-----SIMULATING POSEUILLE FLOW-----" << std::endl;


    poiseuille_flowAnalytical(ux_star);
    save_scalar("UxTrial", ux_star, 0);
    //taylor_green(0, rho, ux, uy);
    //Initialization functions
    init_distribution(f, rho, ux, uy);
    init_wall_distribution(f_wall);
    init_Geo(geo);
    //save_scalar("geo", geo, 0);

    //Time loop
    for(int i = 0; i<NSTEPS; i++){

        compute_rho_u(f, rho, ux, uy);
        
        collison(f, rho, ux, uy, source, feq, ft);

        stream(f, ft);

        bounce_back(geo, f, ft, f_wall);

    }

    
    //Save results
    save_scalar("Ux", ux, NSTEPS);
    save_scalar("Uy", uy, NSTEPS);
    //save_scalar("Ux", ux, NSTEPS);


    std::cout << "DOMAIN SIZE: " << NX << " X " << NY << std::endl;
    std::cout << "MAX VELOCITY IN X DIRECTION: " << ux[scalar_index(0, 20)] << " m/s" << std::endl;
    std::cout << "MAX VELOCITY IN X DIRECTION (Analytical): " << ux_star[scalar_index(0, 20)] << " m/s" << std::endl;


    //std :: cout << f[field_index(5, 6, 5)] << std:: endl;


}






=======

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

//using namespace std;

//Domain size
const int NX = 40;
const int NY = 40;

//Initialize weights
const double w0 = 4.0/9.0;  // zero weight
const double ws = 1.0/9.0;  // adjacent weight
const double wd = 1.0/36.0; // diagonal weight

//initialize direction and weight vectors
const int ndir = 9;
const double wi[] = {w0,ws,ws,ws,ws,wd,wd,wd,wd};
const int dirx[] = {0,1,0,-1, 0,1,-1,-1, 1};
const int diry[] = {0,0,1, 0,-1,1, 1,-1,-1};

/*

## Direction Scheme ##

    6 2 5
    3 0 1
    7 4 8

*/

const int opposite[] = {0, 3, 4, 1, 2, 7, 8, 5, 6}; //For storing opposite directions used in bounce back

//const double nu = 1.0/6.0;
const double tau = 0.6;
const double nu = (tau - 0.5)/3;
double tauinv = 1/tau;
const double u_max = 0.02;


const double rho0 = 1.0;    //Initial density
const double mu = rho0*nu;

//Number of timesteps
const unsigned int NSTEPS = 20000;
const unsigned int NSAVE  =  50;


//For Poseuille's flow
const double dpdx = 1.0e-3;

int field_index(int x, int y, int d){
    
    return (NX*(NY*d+y)+x);

}

int scalar_index(int x, int y){
    
    return (NX*y+x);

}

void taylor_green(int t, double *rho, double *ux, double *uy){

    double kx = (2*3.14159)/NX;
    double ky = (2*3.14159)/NY;
    double td = 1/(nu*(kx*kx + ky*ky));

    for(int x = 0; x < NX; x++){

        for(int y = 0; y < NY; y++){
            
            double X = x+0.5;
            double Y = y+0.5;

            ux[scalar_index(x, y)] = -u_max*sqrt(ky/kx)*cos(kx*X)*sin(ky*Y)*exp(-1*t/td);
            uy[scalar_index(x, y)] = -u_max*sqrt(kx/ky)*sin(kx*X)*cos(ky*Y)*exp(-1*t/td);

            double P = -0.25*rho0*u_max*u_max*((ky/kx)*cos(2*kx*X) + (kx/ky)*cos(2.0*ky*Y))*exp(-2*t/td);

            rho[scalar_index(x, y)] = rho0 + 3*P;

        }
    }
}

void poiseuille_flowAnalytical(double *ux){

    double velConst = -(dpdx)/(2*mu);

    //printf("%lf", velConst);

    for(int x = 0; x < NX; x++){

        for(int y = 0; y < NY; y++){

            int r = y - NY + 1;
            ux[scalar_index(x, y)] = velConst*(y)*r;

        }

    }
}


//Initialization function used for Poseuille's flow
void init_distribution(double *f, double *r, double *u, double *v){
    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){

            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            double rho = r[scalar_index(x, y)];


            for(int i = 0; i<ndir; i++){

                double cidotu = dirx[i]*ux + diry[i]*uy;
                f[field_index(x,y,i)] = wi[i]*rho0;

            }
        }

    }
}


void init_wall_distribution(double *f_wall){

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<2; y++){
            for(int i = 0; i<ndir; i++){
                f_wall[field_index(x, y, i)] = wi[i]*rho0;   
            }
        }
    }
}


void init_Geo(double *geo){
    for(int y = 0; y<NY; y++){
        for(int x = 0; x<NX; x++){
            if(y == 0 || y == NY-1){
                //Boundary nodes are labelled as 1
                geo[scalar_index(x, y)] = 1;
            }

            else{
                //Fluid nodes are labelled as 0
                geo[scalar_index(x, y)] = 0;
            }

            //std::cout << geo[scalar_index(x,y)];
            
        }

        //std::cout << std::endl;

        
    }
}


void compute_rho_u(double *f, double *r, double *u, double *v)
{
    for(int y = 0; y < NY; ++y)
    {
        for(int x = 0; x < NX; ++x)
        {
            double rho = 0.0;
            double ux  = 0.0;
            double uy  = 0.0;
            
            for(unsigned int i = 0; i < ndir; ++i)
            {
                rho += f[field_index(x,y,i)];
                //std::cout << rho << std::endl;
                ux  += dirx[i]*f[field_index(x,y,i)];
                uy  += diry[i]*f[field_index(x,y,i)];

            }

            ux += dpdx/2;   //Add this term for poseuille flow (pressure difference)
            
            r[scalar_index(x,y)] = rho;
            u[scalar_index(x,y)] = ux/rho;
            v[scalar_index(x,y)] = uy/rho;

            //std::cout << ux << std :: endl;

        }
    }
}



void collison(double *f, double *r, double *u, double *v, double *source, double *feq, double *ft){

    double omtau = (1 - (1/tau)); 

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){
            double ux = u[scalar_index(x, y)];
            double uy = v[scalar_index(x, y)];
            double rho = r[scalar_index(x, y)];

            double u2 = ux*ux + uy*uy;

            for(int i = 0; i<ndir; i++){

                double t1 = ux*dirx[i] + uy*diry[i];
                double t2 = t1*t1;

                source[i] = (1 - 0.5/tau)*wi[i]*(3*(dirx[i] - ux) + 9*(t1)*dirx[i])*dpdx;   //Added for Poseuille flow (External force)

                feq[i] = wi[i]*rho*(1 + 3*t1 + 4.5*t2 - 1.5*u2);
                ft[field_index(x, y, i)] = f[field_index(x, y, i)]*omtau  + feq[i]*tauinv + source[i];

                
            }

            
        }
            
    }

}


void stream(double *f, double *ft){

    //Function streams periodically

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){
            for(int i = 0; i < ndir; i++){

                unsigned int xmd = (NX+x+dirx[i])%NX;
                unsigned int ymd = (NY+y+diry[i])%NY;
                
                f[field_index(xmd,ymd,i)] = ft[field_index(x,y,i)];
                
            }
        }
    }
}




void bounce_back(double *geo, double *f, double *ft, double *f_wall){

    for(int y = 0; y<NY; y++){
        for(int x = 0; x<NX; x++){

            if(geo[scalar_index(x, y)] == 0){

                for(int i = 0; i<ndir; i++){

                    int xmd = (NX + x + dirx[i])%NX;
                    int ymd =  (NY + y + diry[i])%NY;


                    if(geo[scalar_index(xmd, ymd)] == 1){

                        f[field_index(x, y, opposite[i])] = ft[field_index(x, y, i)];

                    }
                }
            }
        }
    }
}



void save_scalar(const char *name, double *scalar, unsigned int n)
{
    const size_t mem_size_scalar = sizeof(double) * NX * NY;

    // assume reasonably-sized file names
    char filename[128];
    char format[16];

    // compute maximum number of digits
    int ndigits = floor(log10((double)NSTEPS) + 1.0);

    // file name format is name0000nnn.bin
    sprintf(format, "%%s%%0%dd.txt", ndigits);
    sprintf(filename, format, name, n); 

    FILE *output = fopen(filename, "w");
    for (int y = 0; y < NY; y++)
    {
        for (int x = 0; x < NX; x++)
        {
            fprintf(output, "%.10f\t", scalar[scalar_index(x, y)]);
        }

        fprintf(output, "\n");
    }
    fclose(output);
}





int main(){

    //Allocate memory

    double *f  = (double*) malloc(sizeof(double)*NX*NY*ndir);
    double *f_wall  = (double*) malloc(sizeof(double)*NX*2*ndir);
    double *geo = (double*) malloc(sizeof(double)*NX*NY);
    double *ux  = (double*) malloc(sizeof(double)*NX*NY);
    double *uy  = (double*) malloc(sizeof(double)*NX*NY);
    double *rho  = (double*) malloc(sizeof(double)*NX*NY);
    double *source  = (double*) malloc(sizeof(double)*ndir);
    double *feq  = (double*) malloc(sizeof(double)*ndir);
    double *ft  = (double*) malloc(sizeof(double)*NX*NY*ndir);
    double *ux_star  = (double*) malloc(sizeof(double)*NX*NY);
    //double *uy_star  = (double*) malloc(sizeof(double)*NX*NY);
    //double *rho_star  = (double*) malloc(sizeof(double)*NX*NY);

    std::cout << "-----SIMULATING POSEUILLE FLOW-----" << std::endl;


    poiseuille_flowAnalytical(ux_star);
    save_scalar("UxTrial", ux_star, 0);
    //taylor_green(0, rho, ux, uy);
    //Initialization functions
    init_distribution(f, rho, ux, uy);
    init_wall_distribution(f_wall);
    init_Geo(geo);
    //save_scalar("geo", geo, 0);

    //Time loop
    for(int i = 0; i<NSTEPS; i++){

        compute_rho_u(f, rho, ux, uy);
        
        collison(f, rho, ux, uy, source, feq, ft);

        stream(f, ft);

        bounce_back(geo, f, ft, f_wall);

    }

    
    //Save results
    save_scalar("Ux", ux, NSTEPS);
    save_scalar("Uy", uy, NSTEPS);
    //save_scalar("Ux", ux, NSTEPS);


    std::cout << "DOMAIN SIZE: " << NX << " X " << NY << std::endl;
    std::cout << "NUMBER OF ITERATIONS: " << NSTEPS << std::endl;
    std::cout << "MAX VELOCITY IN 'X' DIRECTION: " << ux[scalar_index(0, 20)] << " m/s" << std::endl;
    std::cout << "MAX VELOCITY IN 'X' DIRECTION (Analytical): " << ux_star[scalar_index(0, 20)] << " m/s" << std::endl;
    std::cout << "BOUNDARY CONDITIONS USED: Half-way Bounce-back for walls, Periodic for inlet and outlet" << std::endl;


    //std :: cout << f[field_index(5, 6, 5)] << std:: endl;


}






>>>>>>> d155345a6f23836c65eedf1dabc02505f9df7864:main.cpp
