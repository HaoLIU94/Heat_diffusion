# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>

/*
[x] = read("/home/hao/ensiieS3/PAP/project/2dimension/square/square_u0001.txt",1000,1000)
set(gcf(),"color_map",[jetcolormap(64);hotcolormap(64)])
zm = min(x); zM = max(x);
colorbar(zm,zM,[1,64])
Sgrayplot(1:1000,1:1000,x,strf="041",rect=[0,0,1000,1000],colminmax=[1,64])
*/
using namespace std;

double *dirichlet_condition ( int node_num, double node_xy[], double time );
double *initial_condition ( int node_num, double node_xy[], double time );
double k_coef ( int node_num, double node_xy[], double time );
double rhs ( int node_num, double node_xy[], double time );

//****************************************************************************80

double *dirichlet_condition ( int node_num, double node_xy[], double time )

//****************************************************************************80
//
//  Purpose:
//
//    DIRICHLET_CONDITION sets the value of a Dirichlet boundary condition.
//
//  Discussion:
//
//    The input points (X,Y) are assumed to lie on the boundary of the
//    region.
//
//    This routine is set for the unit square.
//
//    We assume that the equation to be solved is
//
//      dUdT - Laplacian U + K * U = F
//
//    with K = 0, and F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
//
//    The exact solution is:
//
//      U = sin(pi*x) * sin(pi*y) * exp(-t)
//
//  Modified:
//
//    04 December 2013
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of points.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
//
//    Input, double TIME, the current time.
//
//    Output, double DIRICHLET_CONDITION[NODE_NUM], the value of 
//    the solution at the the point.
//
{
  int node;
  double *u;

  u = new double[node_num];

  for ( node = 0; node < node_num; node++ )
  {
    u[node] = 250.15;
  }

  return u;
}
//****************************************************************************80

double *initial_condition ( int node_num, double node_xy[], double time )

//****************************************************************************80
//
//  Purpose:
//
//    INITIAL_CONDITION sets the initial condition.
//
//  Discussion:
//
//    The input value TIME is assumed to be the initial time.
//
//    We assume that the equation to be solved is
//
//      dUdT - Laplacian U + K * U = F
//
//    with K = 0, and F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
//
//    The exact solution is:
//
//      U = sin(pi*x) * sin(pi*y) * exp(-t)
//
//  Modified:
//
//    08 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of points.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
//
//    Input, double TIME, the current time.
//
//    Output, double INITIAL_CONDITION[NODE_NUM], the value of the 
//    solution at the the initial time.
//
{
  int node;
  double *u;

  u = new double[node_num];

  for ( node = 0; node < node_num; node++ )
  {
    u[node] = 250.15;
  }

  return u;
}
//****************************************************************************80

double k_coef ( int node_num, double node_xy[], double time )

//****************************************************************************80
//
//  Purpose:
//
//    K_COEF evaluates the coefficient K(X,Y,T) function.
//
//  Discussion:
//
//    Right now, we are assuming that NODE_NUM is always 1!
//
//    We assume that the equation to be solved is
//
//      dUdT - Laplacian U + K * U = F
//
//    with K = 0, and F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
//
//    The exact solution is:
//
//      U = sin(pi*x) * sin(pi*y) * exp(-t)
//
//  Modified:
//
//    08 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of points.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
//
//    Input, double TIME, the current time.
//
//    Output, double K_COEF, the value of the coefficient.
//
{
  double k;
  //fer={80.2, 7874.0, 440.0}
  k = 80.2/(7874*440);

  return k;
}
//****************************************************************************80

double rhs ( int node_num, double node_xy[], double time )

//****************************************************************************80
//
//  Purpose:
//
//    RHS gives the right-hand side of the differential equation.
//
//  Discussion:
//
//    Right now, we are assuming that NODE_NUM is always 1!
//
//    We assume that the equation to be solved is
//
//      dUdT - Laplacian U + K * U = F
//
//    with 
//
//      K = 0, 
//
//    and 
//
//      F = (2*pi*pi-1)*sin(pi*x)*sin(pi*y)*exp(-t).
//
//    The exact solution is:
//
//      U = sin(pi*x) * sin(pi*y) * exp(-t)
//
//  Modified:
//
//    08 January 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int NODE_NUM, the number of points.
//
//    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
//
//    Input, double TIME, the current time.
//
//    Output, double RHS, the value of the right hand side.
//
{
  int node;
  double f;

  f = 250.15;
  for ( node = 0; node < node_num; node++ )
  {
    u[node] = sin ( pi * node_xy[0+node*2] ) 
            * sin ( pi * node_xy[1+node*2] ) * exp ( - time );
  }
  return f;
}
