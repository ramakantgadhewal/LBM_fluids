#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>

//using namespace std;


//Domain size
const int NX = 200;
const int NY = 120;

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
const double tau = 0.55;
const double nu = (tau - 0.5)/3;
double tauinv = 1/tau;
const double u_max = 0.02;


const double rho0 = 1.0;    //Initial density
const double mu = rho0*nu;


//Number of timesteps
const unsigned int NSTEPS = 5000;
const unsigned int NSAVE  =  2;


//For Poseuille's flow
const double dpdx = 1.0e-3;


//Define problem parameters for flow past cylinder simulation

double Re = 100;
//Inlet = 15D => NY/30 = cylinder Radius
int cylinderRadius = NY / 30;
double Lc = 2*cylinderRadius;
double u_inlet = Re*nu/Lc;
double rho_out = 1.0;
int Xc = NX / 5;
int Yc = NY / 2;




int field_index(int x, int y, int d){
    
    return (NX*(NY*d+y)+x);

}

int scalar_index(int x, int y){
    
    return (NX*y+x);

}

double distance(int x, int y)
{
    double xdist = pow((x - Xc), 2);
    double ydist = pow((y - Yc), 2);
    double dist = sqrt(xdist + ydist);
    return dist;
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

            else if(x == 0){

                geo[scalar_index(x, y)] = 8;

            }

            else if(x == NX-1){

                geo[scalar_index(x, y)] = 9;

            }

            else if (distance(x, y) <= cylinderRadius) {
                geo[scalar_index(x, y)] = 4;
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


void compute_rho_u(double *f, double *r, double *u, double *v, double *geo, double *v_net){

    for(int y = 0; y < NY; ++y){

        for(int x = 0; x < NX; ++x){

            double rho = 0.0;
            double ux  = 0.0;
            double uy  = 0.0;

            for(unsigned int i = 0; i < ndir; ++i) {
                rho += f[field_index(x,y,i)];
                ux  += dirx[i]*f[field_index(x,y,i)];
                uy  += diry[i]*f[field_index(x,y,i)];

            }
            


            //ux += dpdx/2;   //Add this term for poseuille flow (pressure difference)
            
            r[scalar_index(x,y)] = rho;
            u[scalar_index(x,y)] = ux/rho;
            v[scalar_index(x,y)] = uy/rho;

            double u2 = (u[scalar_index(x, y)])*(u[scalar_index(x, y)]);
            double v2 = (v[scalar_index(x, y)])*(v[scalar_index(x, y)]);

            v_net[scalar_index(x, y)] = sqrt(u2 + v2);

            //std::cout << ux << std :: endl;

        }
    }
}


void collison(double *f, double *r, double *u, double *v, double *feq, double *ft){

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

                //source[i] = (1 - 0.5/tau)*wi[i]*(3*(dirx[i] - ux) + 9*(t1)*dirx[i])*dpdx;   //Added for Poseuille flow (External force)

                feq[i] = wi[i]*rho*(1 + 3*t1 + 4.5*t2 - 1.5*u2);
                ft[field_index(x, y, i)] = f[field_index(x, y, i)]*omtau  + feq[i]*tauinv ; //+ source[i];

                
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


void mid_bounce_back(double *geo, double *f, double *ft){

    for(int y = 0; y<NY; y++){
        for(int x = 0; x<NX; x++){

            if(geo[scalar_index(x, y)] == 0){
                //Check for fluid nodes

                for(int i = 0; i<ndir; i++){

                    int xmd = (NX + x + dirx[i])%NX;
                    int ymd =  (NY + y + diry[i])%NY;


                    if(geo[scalar_index(xmd, ymd)] == 1 || geo[scalar_index(xmd, ymd)] == 4){
                        //If adjacent node is on wall

                        f[field_index(x, y, opposite[i])] = ft[field_index(x, y, i)];

                    }

                }
            }
        }
    }
}


void compute_forces(double *geo, double *f, double *ft, double *r, double *v_net, double *Cd, int i){

    double px, py, cd, cl;

    for(int y = 0; y<NY; y++){
        for(int x = 0; x<NX; x++){

            if(geo[scalar_index(x, y)] == 0){
                //Check for fluid nodes

                for(int i = 0; i<ndir; i++){

                    int xmd = x + dirx[i];
                    int ymd = y + diry[i];

                    if(geo[scalar_index(xmd, ymd)] == 4){
                        //If adjacent node is on wall

                        //double px, py;

                        px += f[field_index(x, y, i)]*dirx[i] - f[field_index(x, y, i)]*dirx[opposite[i]];
                        py += f[field_index(x, y, i)]*diry[i] - f[field_index(x, y, i)]*diry[opposite[i]];

                        //double v_net2 = v_net[scalar_index(x, y)]*v_net[scalar_index(x, y)];

                    
                    }

                }
            }
        }
    }


    cd = (px)/(0.5*rho0*(2*cylinderRadius)*u_inlet*u_inlet);
    cl = (py)/(0.5*rho0*(2*cylinderRadius)*u_inlet*u_inlet);
    //std::cout << "DRAG COEFFICIENT :  "<< cd << std::endl;
    //std::cout << "LIFT COEFFICIENT :  "<< cl << std::endl;

    //Cd[i] = cd;
    Cd[i] = cd;

}


void zou_inlet(double *geo, double *f, double *ft, double *rho, double u_inlet){

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){
            
            if(geo[scalar_index(x, y)] == 8){

                double rhod = (1/(1-u_inlet))*(f[field_index(x, y, 0)] + f[field_index(x, y, 2)] + f[field_index(x, y, 4)] + 2*(f[field_index(x, y, 3)] + f[field_index(x, y, 6)] + f[field_index(x, y, 7)]));

                f[field_index(x, y, 1)] = f[field_index(x, y, 3)] + (2.0/3.0)*rhod*u_inlet;
                f[field_index(x, y, 5)] = f[field_index(x, y, 7)] - (1.0/2.0)*(f[field_index(x, y, 2)] - f[field_index(x, y, 4)]) + (1.0/6.0)*rhod*u_inlet;
                f[field_index(x, y, 8)] = f[field_index(x, y, 6)] + (1.0/2.0)*(f[field_index(x, y, 2)] - f[field_index(x, y, 4)]) + (1.0/6.0)*rhod*u_inlet;

            }
        }
    }
}


void zou_outlet(double *geo, double *f, double *ft, double *rho, double u_inlet){

    for(int x = 0; x<NX; x++){
        for(int y = 0; y<NY; y++){
            if(geo[scalar_index(x, y)] == 9){

                double ux = 1 - ((2*(f[field_index(x, y, 1)] + f[field_index(x, y, 1)] + f[field_index(x, y, 1)]) + f[field_index(x, y, 0)] + f[field_index(x, y, 2)] + f[field_index(x, y, 4)])/rho_out);

                f[field_index(x, y, 3)] = f[field_index(x, y, 1)] - (2.0/3.0)*rho_out*ux;
                f[field_index(x, y, 7)] = f[field_index(x, y, 5)] + (1.0/2.0)*(f[field_index(x, y, 2)] - f[field_index(x, y, 4)]) - (1.0/6.0)*rho_out*ux;
                f[field_index(x, y, 6)] = f[field_index(x, y, 8)] - (1.0/2.0)*(f[field_index(x, y, 2)] - f[field_index(x, y, 4)]) - (1.0/6.0)*rho_out*ux;

            }
        }
    }
}


void outflow(double *geo, double *f, double *ft){

    for(int x = 0; x<NX; x++){
        for(int y = 0; y < NY; y++){
            if(geo[scalar_index(x, y)] == 9){

                f[field_index(x, y, 3)] = f[field_index(x-1, y, 3)];
                f[field_index(x, y, 6)] = f[field_index(x-1, y, 6)];
                f[field_index(x, y, 7)] = f[field_index(x-1, y, 7)];


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


void save_scalar_linear(const char *name, double *scalar, unsigned int n)
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
    for (int y = 0; y < NSTEPS; y++)
    {
        fprintf(output, "%.10f\t", scalar[y]);
         
    }
    fclose(output);
}




int main(){

    //Allocate memory

    double *f  = (double*) malloc(sizeof(double)*NX*NY*ndir);
    //double *f_wall  = (double*) malloc(sizeof(double)*NX*2*ndir);
    double *geo = (double*) malloc(sizeof(double)*NX*NY);
    double *ux  = (double*) malloc(sizeof(double)*NX*NY);
    double *uy  = (double*) malloc(sizeof(double)*NX*NY);
    double *rho  = (double*) malloc(sizeof(double)*NX*NY);
    //double *source  = (double*) malloc(sizeof(double)*ndir);
    double *feq  = (double*) malloc(sizeof(double)*ndir);
    double *ft  = (double*) malloc(sizeof(double)*NX*NY*ndir);
    //double *ux_star  = (double*) malloc(sizeof(double)*NX*NY);
    double *v_net  = (double*) malloc(sizeof(double)*NX*NY);
    double *Cd  = (double*) malloc(sizeof(double)*NSTEPS);
    //double *Cl  = (double*) malloc(sizeof(double)*NSTEPS);
    //double *uy_star  = (double*) malloc(sizeof(double)*NX*NY);
    //double *rho_star  = (double*) malloc(sizeof(double)*NX*NY);

    std::cout << "-----SIMULATING FLOW PAST CYLINDER-----" << std::endl;


    //taylor_green(0, rho, ux, uy);
    //Initialization functions
    init_distribution(f, rho, ux, uy);
    //init_wall_distribution(f_wall);
    init_Geo(geo);
    //save_scalar("geo", geo, 0);


    //Time loop
    for(int i = 0; i<NSTEPS; i++){

        compute_rho_u(f, rho, ux, uy, geo, v_net);

        /*if(((i+1)%NSAVE == 0) && (i >= 2000)) {

            save_scalar("VNet", v_net, i+1);
        
        }*/
        
        collison(f, rho, ux, uy,feq, ft);

        stream(f, ft);

        zou_inlet(geo, f, ft, rho, u_inlet);
        //zou_outlet(geo, f, ft, rho, u_inlet);
        
        outflow(geo, f, ft);

        mid_bounce_back(geo, f, ft);

        compute_forces(geo, f, ft, rho, v_net, Cd, i);

        if(i >= 1800) {

            save_scalar_linear("Cd", Cd, 0);
        
        }


    }


    //compute_forces(geo, f, ft, rho, v_net);

    
    //Save results
    //save_scalar("Ux", ux, NSTEPS);
    //save_scalar("Uy", uy, NSTEPS);
    //std:: cout << rho[scalar_index(39, 8)] << std::endl;
    save_scalar("VNet", v_net, NSTEPS);


    std::cout << "DOMAIN SIZE: " << NX << " X " << NY << std::endl;
    std::cout << "NUMBER OF ITERATIONS: " << NSTEPS << std::endl;
    std::cout << "REYNOLDS NUMBER: " << Re << std::endl;
    std::cout << "INLET VELOCITY: " << u_inlet << " m/s" << std::endl;
    //std::cout << "# # # # # # # # # # # # # # # # # # # #" << std::endl;
    //std::cout << "BOUNDARY CONDITIONS USED: Half-way Bounce-back for walls, Periodic for inlet and outlet" << std::endl;



    //std :: cout << f[field_index(5, 6, 5)] << std:: endl;


}


/*

   % momentum exchange
    net_del_Px = 2*sum(repmat(cx,length(fluid_boundaryNodes),1).*((f(:,fluid_boundaryNodes)').*lat_sol_nodes_fluidNodeNe_logical),2);
    net_del_Py = 2*sum(repmat(cy,length(fluid_boundaryNodes),1).*((f(:,fluid_boundaryNodes)').*lat_sol_nodes_fluidNodeNe_logical),2);

%     neigh_x = repmat(cx,length(boundarySolidNodes),1);
%     neighxindicator_x = neigh_x.*nodeNe_boundarySolidNodes_logical;
%     P_x = 2*neighxindicator_x.*f(:,boundarySolidNodes)';
%     net_del_Px = sum(P_x,2);
%     
%     neigh_y = repmat(cy,length(boundarySolidNodes),1);
%     neighxindicator_y = neigh_y.*nodeNe_boundarySolidNodes_logical;
%     P_y = 2*neighxindicator_y.*f(:,boundarySolidNodes)';
%     net_del_Py = sum(P_y,2);
    
    

    % forces in x,y
    Fx = [Fx,sum(net_del_Px)];
    Fy = [Fy,sum(net_del_Py)];
    
    % conversion of u_star to u and v_star to v
    u = u_star*delx/delt;
    v = v_star*delx/delt;
    
    % pressure calculation
    p_star = cs_star^2*rho_star;
    
    % pressure conversion
    p = p_star*rho_phy*delx^2/delt^2;
    toc
    
    % force conversion
    Fx_phy = Fx*rho_phy*delx^4/delt^2;
    Fy_phy = Fy*rho_phy*delx^4/delt^2;
    
    % drag and lift coef
%     cd = 2*Fx_phy/((u0^2)*rho_phy*d);
%     cl = 2*Fy_phy/((u0^2)*rho_phy*d);

    cd = 2*Fx_phy/((((Uave_star*delx/delt)^2)*rho_phy*d_star*delx));
    cl = 2*Fy_phy/((((Uave_star*delx/delt)^2)*rho_phy*d_star*delx));





*/