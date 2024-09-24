/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
#include<iostream>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<sstream>
#include<string>
#include<chrono>
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
using namespace std;
int i, j;
const int     nx = 100;
const int     ny = 100;
double  u[nx][ny+1], uc[nx][ny], ustar[nx][ny+1],
        v[nx+1][ny], vc[nx][ny], vstar[nx+1][ny],
        p[nx+1][ny+1], pc[nx][ny],
        m[nx+1][ny+1];
double unew, vnew, pnew;
double ue,uw,un,us,ve,vw,vn,vs, ute, utw, vtn, vts;
int        step    = 1;
double        X       = 1.0;
double        Y       = 1.0;
double        dx      = X/(nx-1);
double        dy      = Y/(ny-1);
double        dt      = 0.001;
double        delta   = 4.5;
double        error   = 1.0;
double        Re      = 100.0;
double        uResidual = 1e-12;
double        vResidual = 1e-12;
double        uChangeMax = 1.0;
double        vChangeMax = 1.0;
double        pChangeMax = 1.0;
double        continuity = 1e-08;
double        running_time = 0.0;
double        simulation_time = 0.0;
string fileName = "Result_for_100_";
string printStep;
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void initial_conditions ()
{
    //initial condition for u
    for (i=0; i<=(nx-1); i++)             //i = west to east
    {
        for (j=0; j<=(ny); j++)           //j = south to north
        {
            u[i][j]         = 0.0;
            u[i][ny]        = 1.0;
            u[i][ny-1]      = 1.0;
        }
    }
    //initial condition for v
    for (i=0; i<=(nx); i++)
    {
        for (j=0; j<=(ny-1); j++)
        {
            v[i][j]         = 0.0;
        }
    }
    //initial condition for p
    for (i=0; i<=(nx); i++)
    {
        for (j=0; j<=(ny); j++)
        {
            p[i][j]         = 1.0;
        }
    }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void velocity_BCs ()
{
     //boundary condition for u velocity
    for (j=1; j<=(ny-1); j++)          //node inside domain
    {
        u[0][j]        = 0.0;         //u at west (dirichlet)
        u[nx-1][j]     = 0.0;        //u at east (dirichlet)
    }
    for (i=0; i<=(nx-1); i++)
    {
        u[i][0]        = -u[i][1] + 2.0*0.0; //u at south wall (dirichlet)
        u[i][ny]       = -u[i][ny-1] + 2.0*1.0; //u at north wall = 1 ((x+y)/2 = 1 >> x = -y + 2*1)
    }
    //boundary condition for v velocity
    for (j=1; j<=(ny-2); j++) //only inside node (exclude bottom and top)
    {
        v[0][j]        = -v[1][j] + 2.0*0.0;  //west
        v[nx][j]       = -v[nx-1][j] + 2.0*0.0;   //east
    }
    for (i=0; i<=(nx); i++)      //all node south to north
    {
        v[i][0]        = 0.0;                //south
        v[i][ny-1]     = 0.0;               //north
    }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void pressure_BCs ()
{
     //boundary condition for pressure
    for (i=1; i<=(nx-1); i++)
    {
        p[i][0]        = p[i][1];             //west
        p[i][ny]       = p[i][ny-1];          //east
    }
    for (j=0; j<=(ny); j++)
    {
        p[0][j]        = p[1][j];             //south
        p[nx][j]       = p[nx-1][j];          //north
    }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void velocity_schemer_u (double& ue, double& uw, double& un, double& us, double& ute, double& utw, double& vtn, double& vts)
{
    ute = (0.5*(u[i][j] + u[i+1][j]));
    utw = (0.5*(u[i-1][j] + u[i][j]));
    vtn = (0.5*(v[i][j] + v[i+1][j]));
    vts = (0.5*(v[i][j-1] + v[i+1][j-1]));
//quick scheme
    if (ute >= 0){ ue = 0.75*u[i][j]   + 0.375*u[i+1][j] - 0.125*u[i-1][j]; }
    else         { ue = 0.75*u[i+1][j] + 0.375*u[i][j]   - 0.125*u[i+2][j]; }
    if (utw >= 0){ uw = 0.75*u[i-1][j] + 0.375*u[i][j]   - 0.125*u[i-2][j]; }
    else         { uw = 0.75*u[i][j]   + 0.375*u[i-1][j] - 0.125*u[i+1][j]; }
    if (vtn >= 0){ un = 0.75*u[i][j]   + 0.375*u[i][j+1] - 0.125*u[i][j-1]; }
    else         { un = 0.75*u[i][j+1] + 0.375*u[i][j]   - 0.125*u[i][j+2]; }
    if (vts >= 0){ us = 0.75*u[i][j-1] + 0.375*u[i][j]   - 0.125*u[i][j-2]; }
    else         { us = 0.75*u[i][j]   + 0.375*u[i][j-1] - 0.125*u[i][j+1]; }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void velocity_schemer_v (double& ve, double& vw, double& vn, double& vs, double& ute, double& utw, double& vtn, double& vts)
{
    ute = (0.5*( u[i][j] + u[i][j+1]));
    utw = (0.5*( u[i-1][j] + u[i-1][j+1]));
    vtn = (0.5*( v[i][j] + v[i][j+1]));
    vts = (0.5*( v[i][j-1] + v[i][j]));
//qiuck scheme
    if (ute >= 0){ ve = 0.75*v[i][j]   + 0.375*v[i+1][j] - 0.125*v[i-1][j]; }
    else         { ve = 0.75*v[i+1][j] + 0.375*v[i][j]   - 0.125*v[i+2][j]; }
    if (utw >= 0){ vw = 0.75*v[i-1][j] + 0.375*v[i][j]   - 0.125*v[i-2][j]; }
    else         { vw = 0.75*v[i][j]   + 0.375*v[i-1][j] - 0.125*v[i+1][j]; }
    if (vtn >= 0){ vn = 0.75*v[i][j]   + 0.375*v[i][j+1] - 0.125*v[i][j-1]; }
    else         { vn = 0.75*v[i][j+1] + 0.375*v[i][j]   - 0.125*v[i][j+2]; }
    if (vts >= 0){ vs = 0.75*v[i][j-1] + 0.375*v[i][j]   - 0.125*v[i][j-2]; }
    else         { vs = 0.75*v[i][j]   + 0.375*v[i][j-1] - 0.125*v[i][j+1]; }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void uv_calculation ()
{
        //calculate ustar (u star/intermediate u velocity)
        for (i=1; i<=(nx-2); i++)         //calculation only inside domain, exclude boundary node
        {
            for (j=1; j<=(ny-1); j++)
            {
                ue=uw=un=us=ute=utw=vtn=vts=0.0;
                velocity_schemer_u(ue, uw, un, us, ute, utw, vtn, vts);
                unew = u[i][j] - (dt/dx)*(ue*ute-uw*utw) - (dt/dy)*(un*vtn-us*vts) - dt/dx*(p[i+1][j]-p[i][j]) + dt/Re*((u[i+1][j]-2*u[i][j]+u[i-1][j])/dx/dx) + dt/Re*((u[i][j+1]-2*u[i][j]+u[i][j-1])/dy/dy);
                if (abs(unew - u[i][j]) > uChangeMax)
                {uChangeMax = abs(unew - u[i][j]);}
                u[i][j] = unew;
            }
        }
        //calculate vstar (v star/intermediate v velocity)
        for (i=1; i<=(nx-1); i++)
        {
            for (j=1; j<=(ny-2); j++)
            {
                ve=vw=vn=vs=ute=utw=vtn=vts=0.0;
                velocity_schemer_v(ve, vw, vn, vs, ute, utw, vtn, vts);
                vnew = v[i][j] - (dt/dx)*(ute*ve-utw*vw) - (dt/dy)*(vn*vtn-vs*vts) - dt/dy*(p[i][j+1]-p[i][j]) + dt/Re*((v[i+1][j]-2*v[i][j]+v[i-1][j])/dx/dx) + dt/Re*((v[i][j+1]-2*v[i][j]+v[i][j-1])/dy/dy);
                if (abs(vnew - v[i][j]) > vChangeMax)
                {vChangeMax = abs(vnew - v[i][j]);}
                v[i][j] = vnew;
            }
        }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void pnew_calculation ()
{
        //calculate new pressure pnew
        for (i=1; i<=(nx-1); i++)
        {
            for (j=1; j<=(ny-1); j++)
            {
               pnew = p[i][j] - dt*delta*( (u[i][j]-u[i-1][j])/dx + (v[i][j]-v[i][j-1])/dy );
               if (abs(pnew - p[i][j]) > pChangeMax)
               {pChangeMax = abs(pnew - p[i][j]);}
               p[i][j] = pnew;
            }
        }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void continuity_results ()
{
        //displaying error (div vel. vector = 0)
        error = 0.0;
        for (i=1; i<=(nx-1); i++)
        {
            for (j=1; j<=(ny-1); j++)
            {
                m[i][j] = (u[i][j]-u[i-1][j])/dx + (v[i][j]-v[i][j-1])/dy  ;
                if (abs(m[i][j]) > error)
                {error = abs(m[i][j]);}
            }
        }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void monitor ()
{if (step == 1 || step % 1000 == 0) {cout << "\n" << setw(12) << "Step" << setw(17) << "Continuity" << setw(17) << "Max. uChange"
                << setw(17) << "Max. vChange"
                << setw(17) << "Max. pChange" << "\n" << "\n";}
        if (step % 100 == 0) { cout << setw(12) << step << setw(17) << error << setw(17) << uChangeMax << setw(17) << vChangeMax
            << setw(17) << pChangeMax << "\n";}
        step = step + 1;
        printStep = to_string(step);
        if ( step % 10000 == 0)
        {for (i=0; i<=(nx-1); i++)
            {for (j=0; j<=(ny-1); j++)
                    {uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
                     vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
                     pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);}}
            ofstream output_inter;
            output_inter.open( fileName + printStep + ".plt");
            if(!output_inter.fail())
            {output_inter << "VARIABLES = X, Y, U, V, P\n"
               "ZONE F=POINT\n"
               "I = "<<nx << " " << "J = " << ny << "\n";
                for (i=0; i<=(nx-1); i++)
                {for (j=0; j<=(nx-1); j++)
                    {double xpos,ypos;
                        xpos = i*dx;
                        ypos = j*dy;
                        output_inter << xpos << " " << ypos << " " << uc[i][j] << " " << vc[i][j] << " " << pc[i][j] << "\n";}}
            output_inter.close();}}}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void final_results ()
{for (i=0; i<=(nx-1); i++)
    { for (j=0; j<=(ny-1); j++)
            {uc[i][j] = 0.5*(u[i][j]+u[i][j+1]);
             vc[i][j] = 0.5*(v[i][j]+v[i+1][j]);
             pc[i][j] = 0.25*(p[i][j]+p[i+1][j]+p[i][j+1]+p[i+1][j+1]);}}
//output data
ofstream output;
output.precision(8);
output.setf(ios::fixed);
output.open("Results_for_final_step.plt");
if(!output.fail())
{output << "VARIABLES = X, Y, U, V, P\n" "ZONE F=POINT\n" "I = "<<nx << " " << "J = " << ny << "\n";
    for (i=0; i<=(nx-1); i++)
    {for (j=0; j<=(ny-1); j++)
        { double xpos,ypos;
          pos = i*dx;
          ypos = j*dy;
          output << xpos << " " << ypos << " " << uc[i][j] << " " << vc[i][j] << " " << pc[i][j] << "\n"; }}
    output.close();
    cout << "Writing file is succesfully" << "\n";}
else
{ cout << "File cannot writed" << "\n";}}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
void case_info ()
{
    ofstream info_end;
    info_end.open("case_info.txt");
    if(!info_end.fail())
    {
        info_end << "Case info" << "\n" << "\n"
                 << setw(20) << "simulation time = " << simulation_time << ' '<< 's' << "\n"
                 << setw(20) << "Running time = " << running_time << ' ' << 's'<< "\n";
        info_end.close();
    }
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
int main ()
{
    chrono::steady_clock::time_point tStart = chrono::steady_clock::now();
    initial_conditions();
    while ( (uChangeMax > uResidual) || (vChangeMax > vResidual) )
    {
        uChangeMax = 0.0;
        vChangeMax = 0.0;
        pChangeMax = 0.0;
        uv_calculation();
        velocity_BCs();
        pnew_calculation();
        pressure_BCs();
        continuity_results();
        monitor();
    }
    chrono::steady_clock::time_point tEnd = chrono::steady_clock::now();
    running_time = chrono::duration_cast<chrono::seconds>(tEnd - tStart).count();
    simulation_time = step*dt;
    cout << "Simulation time = " << simulation_time << ' ' << "s" << "\n";
    cout << "Running time = " << running_time << ' ' << 's' << "\n";
    final_results();
    case_info();
    return 0;
}
/*--------------------------------------Program CFD Berbasis Artificial Compressibility Method-----------------------------------*/
