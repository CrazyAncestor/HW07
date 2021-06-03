#include "classes.h"
#define PI 3.14159265359


double f( const double x, const double y )
{   
    return sin(x) * sin(y) + background_pot;
} // FUNCTION : f

extern double absolute(double x);
void solved( matrix &m )
{
    const int dim = m.get_dim();
    const int dim_2 = dim*dim;
    
    for ( int idx = 0; idx < dim_2; idx++ )
    {
        const int i = idx / dim;
        const int j = idx % dim;
        
        const double h = m.get_h();
        const double x = h*i;
        const double y = h*j;
      
        m.input_answer( i, j, f(x, y) );
    } // for ( int idx = 0; idx < dim_2; idx++ )

} // FUNCTION : solved

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

void SOR(double **pot, double **rho,int dim,double h, double omega, bool odd )
{

    if(odd){
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                if ((i + j) % 2 == 1) {
                    if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1) {
                        pot[i][j] += omega * 0.25 * (pot[i + 1][j] + pot[i - 1][j] + pot[i][j + 1] + pot[i][j - 1] - pot[i][j] * 4 - h * h * rho[i][j]);
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
                        pot[i][j] += omega * 0.25 * (pot[i + 1][j] + pot[i - 1][j] + pot[i][j + 1] + pot[i][j - 1] - pot[i][j] * 4 - h * h * rho[i][j]);
                    } //if (i != 0 && i != dim - 1 && j != 0 && j != dim - 1)
                } //if ((i + j) % 2 == 0)
            } //for (int j = 0; j < dim; j++)
        } //for (int i = 0; i < dim; i++) 
    }    
}

void display(double **matrix,int dim)
{
    for( int i=0; i<dim; i++ )
    {
        for( int j=0; j<dim; j++ )
        {
            cout << matrix[i][j] << " ";
        }
        cout << endl;
    }
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
int main( int argc, char *argv[] )
{
    int NRank, MyRank;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &MyRank );
    MPI_Comm_size( MPI_COMM_WORLD, &NRank );
    printf( "Hello World on rank %d/%d\n", MyRank, NRank );
    

    //Set Physical Parameters
    int N_steps = 1000;
    int N = 10;
    double h = PI/(N-1);
    double SOR_omega = 1.9;
    
    bool option_parallel = 0;
    
    //Parallelized
    if(option_parallel)
    {
        //Initialize Potential and Density
        double **potential01,**potential02,**density,**ans;
        potential01 = alloc_2d_init(N,N);
        potential02 = alloc_2d_init(N,N);
        density   = alloc_2d_init(N,N);
        ans         = alloc_2d_init(N,N);

        init_potential(potential01,h,N);
        init_potential(potential02,h,N);
        init_density(density,h,N);
        solved(ans,N,h);
        
        //Set Parallelization parameter
        bool odd;
        if(MyRank==0)odd = 0;
        else         odd = 1;

        const int Count=N*N, TargetRank=(MyRank+1)%2, Tag=123;\
        const int NReq= 2;
        MPI_Request Req[NReq];

        //Begin SOR iteration
        for(int steps = 0;steps<N_steps;steps++){
            SOR(potential01,density,N,h,SOR_omega,odd);
            MPI_Irecv(  &(potential02[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, &Req[0] );
            MPI_Isend( &(potential01[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, &Req[1] );
            
            MPI_Waitall( NReq, Req, MPI_STATUSES_IGNORE );
           
            SOR(potential02,density,N,h,SOR_omega,odd);

            MPI_Irecv( &(potential01[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, &Req[0] );
            MPI_Isend( &(potential02[0][0]), Count, MPI_DOUBLE, TargetRank, Tag, MPI_COMM_WORLD, &Req[1] );
            MPI_Waitall( NReq, Req, MPI_STATUSES_IGNORE );
        }
        
       
        if(MyRank==0)Error(potential02,ans,N);
    }
    
    
    
    else
    {//Initialize Potential and Density
        double **potential01,**potential02,**density,**ans;
        potential01 = alloc_2d_init(N,N);
        potential02 = alloc_2d_init(N,N);
        density   = alloc_2d_init(N,N);
        ans         = alloc_2d_init(N,N);

        init_potential(potential01,h,N);
        init_potential(potential02,h,N);
        init_density(density,h,N);
        solved(ans,N,h);
        
        //Set Parallelization parameter
        bool odd,even;
        odd = 0;
        even = 1;

        //Begin SOR iteration
        for(int steps = 0;steps<N_steps*2;steps++){
            SOR(potential01,density,N,h,SOR_omega,odd);
            SOR(potential01,density,N,h,SOR_omega,even);
        }
        
       
        if(MyRank==0)Error(potential01,ans,N);
    }
    MPI_Finalize();
} // FUNCTION : main 
