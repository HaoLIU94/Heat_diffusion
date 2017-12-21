# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;

int main ( );
void dtable_data_write ( ofstream &output, int m, int n, double table[] );
void dtable_write ( string output_filename, int m, int n, double table[], 
  bool header );
void f ( double a, double b, double t0, double t, int n, double x[], 
  double value[] );
int decomposition ( int n, double a[] );
double *resolution ( int n, double a_lu[], double b[], int job );
void timestamp ( );
void u0 ( double a, double b, double t0, int n, double x[], double value[] );
double ua ( double a, double b, double t0, double t );
double ub ( double a, double b, double t0, double t );

//****************************************************************************80

int main ( )
{
  double *a;
  double *b;
  double *fvec;
  bool header;
  int i;
  int info;
  int j;
  int job;
  double k;
  double *t;
  double t_delt;
  string t_file;
  double t_max;
  double t_min;
  int t_num;
  double *u;
  string u_file;
  double w;
  double *x;
  double x_delt;
  string x_file;
  double x_max;
  double x_min;
  int x_num;

  timestamp ( );

  struct material
   {
	double lamda;
 	double den;
 	double c;
   }; 
  	
  material cuivre={389.0, 8940.0, 380.0};
  material fer={80.2, 7874.0, 440.0};
  material verre={1.2, 2530.0, 840.0};
  material ploystyrene={0.1, 1040.0, 1200.0};
  
  cout << " Choose between \n cuivre,\n fer,\n verre,\n ploystyrene\n";
  cout << " Please input your choice\n";
  string choice;
  cin >> choice;
  if (choice.compare("cuivre"))
  	k = cuivre.lamda/(cuivre.den*cuivre.c);
  if (choice.compare("fer"))
  	k = fer.lamda/(fer.den*fer.c);
  if (choice.compare("verre"))
  	k = verre.lamda/(verre.den*cuivre.c);
  if (choice.compare("ploystyrene"))
  	k = ploystyrene.lamda/(ploystyrene.den*ploystyrene.c);
//
//  Set X values.
//
  x_min = 0.0;
  x_max = 1.0;
  x_num = 1001;
  x_delt = ( x_max - x_min ) / ( double ) ( x_num - 1 );

  x = new double[x_num];
// 
//  Set T values.
//
  t_min = 0.0;
  t_max = 16000.0;
  t_num = 100;
  t_delt = ( t_max - t_min ) / ( double ) ( t_num - 1 );

  t = new double[t_num];
//
//  Set the initial data, for time T_MIN.
//
  u = new double[x_num*t_num];

  u0 ( x_min, x_max, t_min, x_num, x, u );
//
//  The matrix A does not change with time.  We can set it once,
//  factor it once, and solve repeatedly.
//
  w = k * t_delt / x_delt / x_delt;

  a = new double[3*x_num];

  a[0+0*3] = 0.0;

  a[1+0*3] = 1.0;
  a[0+1*3] = 0.0;

  for ( i = 1; i < x_num - 1; i++ )
  {
    a[2+(i-1)*3] =           - w;
    a[1+ i   *3] = 1.0 + 2.0 * w;
    a[0+(i+1)*3] =           - w;
  }

  a[2+(x_num-2)*3] = 0.0;
  a[1+(x_num-1)*3] = 1.0;

  a[2+(x_num-1)*3] = 0.0;
//
//  Factor the matrix.
//
  info = decomposition ( x_num, a );

  b = new double[x_num];
  fvec = new double[x_num];

  for ( j = 1; j < t_num; j++ )
  {
//
//  Set the right hand side B.
//
    b[0] = ua ( x_min, x_max, t_min, t[j] );

    f ( x_min, x_max, t_min, t[j-1], x_num, x,fvec );

    for ( i = 1; i < x_num - 1; i++ )
    {
      b[i] = u[i+(j-1)*x_num] + t_delt * fvec[i]/(1040.0*1200.0);
    }

    b[x_num-1] = ub ( x_min, x_max, t_min, t[j] );

    delete [] fvec;

    job = 0;
    fvec = resolution( x_num, a, b, job );

    for ( i = 0; i < x_num; i++ )
    {
      u[i+j*x_num] = fvec[i];
    }
  }

  x_file = "x.txt";
  header = false;
  dtable_write ( x_file, 1, x_num, x, header );

  cout << "\n";
  cout << "  X data written to \"" << x_file << "\".\n";

  t_file = "t.txt";
  header = false;
  dtable_write ( t_file, 1, t_num, t, header );

  cout << "  T data written to \"" << t_file << "\".\n";

  u_file = "u.txt";
  header = false;
  dtable_write ( u_file, x_num, t_num, u, header );

  cout << "  U data written to \"" << u_file << "\".\n";

  cout << "\n";
  cout << "FD1D_HEAT_IMPLICIT\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  delete [] a;
  delete [] b;
  delete [] fvec;
  delete [] t;
  delete [] u;
  delete [] x;

  return 0;
}

//****************************************************************************80

void f ( double a, double b, double t0, double t, int n, double x[],
  double value[] )
{
  int i;
  for ( i = n/10; i < n/5; i++ )
  {
    value[i] = 16*(80+237.15)*(80+237.15);
  }
  for ( i = n/2; i < 3*n/5; i++ )
  {
    value[i] = 12*(80+237.15)*(80+237.15);
  }
  return;
}
//****************************************************************************80

int decomposition ( int n, double a[] )
{
  int i;

  for ( i = 1; i <= n-1; i++ )
  {
//
//  Store the multiplier in L.
//
    a[2+(i-1)*3] = a[2+(i-1)*3] / a[1+(i-1)*3];
//
//  Modify the diagonal entry in the next column.
//
    //a[0+i*3] = a[2+(i-1)*3];

    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3];
  }

  return 0;
}
//****************************************************************************80

double *resolution ( int n, double a_lu[], double b[], int job )
{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    for ( i = 1; i < n; i++ )
    {
      x[i] = x[i] - a_lu[2+(i-1)*3] * x[i-1];
    }
//
//  Solve U * X = Y.
//
    for ( i = n; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( 1 < i )
      {
        x[i-2] = x[i-2] - a_lu[0+(i-1)*3] * x[i-1];
      }
    }
  }
  else
  {
//
//  Solve U' * Y = B
//
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( i < n )
      {
        x[i] = x[i] - a_lu[0+i*3] * x[i-1];
      }
    }
//
//  Solve L' * X = Y.
//
    for ( i = n-1; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] - a_lu[2+(i-1)*3] * x[i];
    }
  }

  return x;
}
//****************************************************************************80

void u0 ( double a, double b, double t0, int n, double x[], double value[] )
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    value[i] = 250.15;
  }
  return;
}
//****************************************************************************80

double ua ( double a, double b, double t0, double t )
{
  double value;

  value = 250.15;

  return value;
}
//****************************************************************************80

double ub ( double a, double b, double t0, double t )
{
  double value;

  value = 250.15;

  return value;
}


//****************************************************************************80

void dtable_data_write ( ofstream &output, int m, int n, double table[] )
{
  int i;
  int j;

  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < m; i++ )
    {
      output << table[i+j*m] << "\t";
    }
    output << "\n";
  }

  return;
}
//****************************************************************************80

void dtable_write ( string output_filename, int m, int n, double table[], 
  bool header )
{
  ofstream output;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "DTABLE_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }

  if ( header )
  {
//  dtable_header_write ( output_filename, output, m, n );
  }

  dtable_data_write ( output, m, n, table );

  output.close ( );

  return;
}

//****************************************************************************80

void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
