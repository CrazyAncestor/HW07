
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <mpi.h>
#include<chrono>
#include "vector"

using namespace std;

#define background_pot 10.
#define PI 3.14159265359


double f( const double x, const double y )
{   
    return sin(x) * sin(y) + background_pot;
} // FUNCTION : f

double absolute(double x){
    if(x>=0)return x;
    return -x;
}
void solved( double ** a,int dim,double h )
{

    const int dim_2 = dim*dim;
    
    for ( int idx = 0; idx < dim_2; idx++ )
    {
        const int i = idx / dim;
        const int j = idx % dim;
        
        const double x = h*i;
        const double y = h*j;
      
        a[i][j] =  f(x, y) ;
    } // for ( int idx = 0; idx < dim_2; idx++ )

} // FUNCTION : solved

//Initialize Density Array
void init_density(double **dens,double h,int dim)
{
    for( int i = 0; i < dim; i++ )
    {
        for( int j = 0; j < dim; j++ )
        {
            double x = h*i;
            double y = h*j;

            dens[i][j] = -2.*sin(x)*sin(y);
        } // for( int j = 0; j < dim; j++ )
    } // for( int i = 0; i < dim; i++ )

} //FUNCTION : init_density

//Initialize Potential Array
void init_potential(double **pot,double h,int dim)
{
    for( int i = 0; i < dim; i++ )
    {
        for( int j = 0; j < dim; j++ )
        {
            pot[i][j] = background_pot;
        } // for( int j = 0; j < dim; j++ )
    } // for( int i = 0; i < dim; i++ )

} //FUNCTION : init_density

void SOR(double **pot_com,double **pot_ref, double **rho,int dim,double h, double omega, bool odd )
{

    if(odd){
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if ((i + j) % 2 == 1) {
                    if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1) {
                        pot_com[i][j] +=  omega * 0.25 * (pot_ref[i + 1][j] + pot_ref[i - 1][j] + pot_ref[i][j + 1] + pot_ref[i][j - 1] - pot_com[i][j] * 4 - h * h * rho[i][j]);
                    } //if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1)
                } //if ((i + j) % 2 == 0)
            } //for (int j = 0; j < dim; j++)
        } //for (int i = 0; i < dim; i++) 
    }              

    else{
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if ((i + j) % 2 == 0) {
                    if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1) {
                        pot_com[i][j] += omega * 0.25 * (pot_ref[i + 1][j] + pot_ref[i - 1][j] + pot_ref[i][j + 1] + pot_ref[i][j - 1] - pot_com[i][j] * 4 - h * h * rho[i][j]);
                    } //if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1)
                } //if ((i + j) % 2 == 0)
            } //for (int j = 0; j < dim; j++)
        } //for (int i = 0; i < dim; i++) 
    }    
}
void SOR(double **pot_com,double **pot_ref, double **rho,int dim,double h, double omega, bool odd,bool upper )
{
    int half = int(dim/2);
    int i_low,i_high;
    if(upper){
        i_low = 0;
        i_high = half;
    }
    else{
        i_low = half;
        i_high = dim;
    }
    if(odd){
        for (int i = i_low; i < i_high; i++) {
            for (int j = 0; j < dim; j++) {
                if ((i + j) % 2 == 1) {
                    if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1) {
                        pot_com[i][j] +=  omega * 0.25 * (pot_ref[i + 1][j] + pot_ref[i - 1][j] + pot_ref[i][j + 1] + pot_ref[i][j - 1] - pot_com[i][j] * 4 - h * h * rho[i][j]);
                    } //if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1)
                } //if ((i + j) % 2 == 0)
            } //for (int j = 0; j < dim; j++)
        } //for (int i = 0; i < dim; i++) 
    }              

    else{
        for (int i = i_low; i < i_high; i++) {
            for (int j = 0; j < dim; j++) {
                if ((i + j) % 2 == 0) {
                    if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1) {
                        pot_com[i][j] += omega * 0.25 * (pot_ref[i + 1][j] + pot_ref[i - 1][j] + pot_ref[i][j + 1] + pot_ref[i][j - 1] - pot_com[i][j] * 4 - h * h * rho[i][j]);
                    } //if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1)
                } //if ((i + j) % 2 == 0)
            } //for (int j = 0; j < dim; j++)
        } //for (int i = 0; i < dim; i++) 
    }    
}

void filling(double **pot, double **pot_rec,int dim, int row )
{
for(int j=0;j<dim;j++){
    pot_rec[0][j] = pot[row][j];
    //cout<<pot_rec[0][j]<<" ";
        
    }
    //cout<<endl;
}
void refilling(double **pot, double **pot_rec,int dim, int row )
{
for(int j=0;j<dim;j++) pot[row][j] = pot_rec[0][j];
}

void Error( double **a,double **b,int dim )
{
    double sum = 0;
    double ave = 0;
    
    for( int i = 0; i < dim; i++ )
    {
        for( int j = 0; j < dim; j++ )
        {
            sum += absolute(a[i][j]-b[i][j]);
            ave += absolute(b[i][j]);
        }
        
    } // for( int i = 0; i < dim; i++ )

    cout << sum/ave << endl;

} //FUNCTION : Error
double **alloc_2d_init(int rows, int cols) {
    double *data = (double *)malloc(rows*cols*sizeof(double));
    double **array= (double **)malloc(rows*sizeof(double*));
    for (int i=0; i<rows; i++){
        array[i] = &(data[cols*i]);}

    return array;
}
void init_rec(double** pot_rec,int dim){
    //pot_rec = new double[dim];
    for( int i = 0; i < dim; i++ )
    {
        pot_rec[0][i]=0.;
        //cout<<pot_rec[0][i]<<" ";
    }
    //cout<<endl;
}
int main( int argc, char *argv[] )
{
    int NRank, MyRank;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );
    MPI_Comm_size( MPI_COMM_WORLD, &NRank );
    printf( "Hello World on rank %d/%d\n", MyRank, NRank );
    

    //Set Physical Parameters
    int N_steps = 100;
    int iter = 1;
    int N = 2000;
    int half = int(N/2);
    double h = PI/(N-1);
    double SOR_omega = 1.9;
    
    bool option_parallel;
    if(NRank==1)option_parallel = 0;
    else option_parallel = 1;
    
    //Set CLock
    auto start = chrono::steady_clock::now();

    //Parallelized
    if(option_parallel)
    {
        //Initialize Potential and Density
        double **pot,**potential_receive,**density,**ans;
        pot = alloc_2d_init(N,N);
        potential_receive = alloc_2d_init(1,N);
        density   = alloc_2d_init(N,N);
        ans         = alloc_2d_init(N,N);

        init_potential(pot,h,N);
        init_rec(potential_receive,N);
        init_density(density,h,N);
        solved(ans,N,h);
        
        //Set Parallelization parameter
        bool odd=1,even=0;
        bool upper =1, lower =0;

        const int Count=N, TargetRank=(MyRank+1)%2, Tag=123;
        const int NReq= 2;
        MPI_Request Req[NReq];

        
        //Begin SOR iteration
        
        for(int steps = 0;steps<N_steps;steps++){
            //First compute odd cases
            if(MyRank==0) SOR(pot,pot,density,N,h,SOR_omega,odd,upper);
            else if(MyRank==1) SOR(pot,pot,density,N,h,SOR_omega,odd,lower);

            //send & receive boundary

            filling(pot,potential_receive,N,half-1);            
            if(MyRank==1) MPI_Recv( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
            else if(MyRank==0) MPI_Ssend( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );  
            refilling(pot,potential_receive,N,half-1);

            filling(pot,potential_receive,N,half);
            if(MyRank==0) MPI_Recv( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
            if(MyRank==1) MPI_Ssend( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );
            refilling(pot,potential_receive,N,half);
            

            //Then compute even cases
            if(MyRank==0) SOR(pot,pot,density,N,h,SOR_omega,even,upper);
            else if(MyRank==1) SOR(pot,pot,density,N,h,SOR_omega,even,lower);

            //send & receive boundary

            filling(pot,potential_receive,N,half-1);            
            if(MyRank==1) MPI_Recv( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
            else if(MyRank==0) MPI_Ssend( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );  
            refilling(pot,potential_receive,N,half-1);

            filling(pot,potential_receive,N,half);
            if(MyRank==0) MPI_Recv( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
            if(MyRank==1) MPI_Ssend( &(potential_receive[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );
            refilling(pot,potential_receive,N,half);

        }

        double **pot_lowerfinal = alloc_2d_init(N,N);
        
        if(MyRank==0) MPI_Recv( &(pot_lowerfinal[0][0]), N*N, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE );
        if(MyRank==1) MPI_Ssend( &(pot[0][0]), N*N, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD );
        if(MyRank==0) {
            for(int i=half;i<N;i++){
            for(int j=0;j<N;j++){
                pot[i][j] = pot_lowerfinal[i][j];
        }
        }
        }

        if(MyRank==0)Error(pot,ans,N);
    }
    
    
    
    else
    {//Initialize Potential and Density
        double **pot,**potential_receive,**density,**ans;
        pot = alloc_2d_init(N,N);
        potential_receive = alloc_2d_init(N,N);
        density   = alloc_2d_init(N,N);
        ans         = alloc_2d_init(N,N);

        init_potential(pot,h,N);
        init_potential(potential_receive,h,N);
        init_density(density,h,N);
        solved(ans,N,h);
        
        //Set Parallelization parameter
        bool odd,even;
        odd = 1;
        even = 0;

        //Begin SOR iteration
        for(int steps = 0;steps<N_steps;steps++){
            SOR(pot,pot,density,N,h,SOR_omega,odd);
            SOR(pot,pot,density,N,h,SOR_omega,even);
        }
        
       
        if(MyRank==0)Error(pot,ans,N);
    }


    //Output Time
    auto elapsed = chrono::steady_clock::now() - start;
    auto sec_double = chrono::duration<double>(elapsed);     // double
    cout<<"Total Time:"<<sec_double.count()<<"s"<<endl;
    MPI_Finalize();
} // FUNCTION : main 
