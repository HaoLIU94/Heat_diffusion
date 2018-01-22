# include <iostream>
# include <fstream>

using namespace std;

void dtable_data_write ( ofstream &output, int m, int n )
{
  int i;
  int j;
  double value = 1.0/1000.0 ;

  for ( j = 0; j < m; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
        output << value*i << "\t";
        output << value*j << "\t";
        output << "\n";
    } 
  }

  return;
}
//****************************************************************************80

void dtable_write ( string output_filename, int m, int n )
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

  dtable_data_write ( output, m, n);

  output.close ( );

  return;
}


void triangulisation()
{
  int x_num = 1001;
  int y_num = 1001;
  int nodes[x_num][y_num];
  int triang[6][2*((x_num-1)/2)*((y_num-1)/2)]; //6 , (1001-1)^2 * 2
  for (int i = 0; i< x_num; i++)
    for (int j = 0; j< x_num; j++)
      nodes[i][j] = i+x_num*j+1;
  /*
  for (int i = 0;i<1001;i++)
    for (int j = 0; j< 1001; j++)
      cout<<nodes[i][j]<<endl;
  */
  for (i = 0; i< 2*((x_num-1)/2)*((y_num-1)/2); i++)
    for (j = 0; j< 6; j++)
      triang[j][i] = 


}


int main ( int argc, char *argv[] )
{
  /*
  int x_num = 1000;
  int y_num = 1000;
  string node_file = "elements.txt";
  dtable_write ( node_file, x_num, y_num );
  cout << " Elements data written to \"" << node_file << "\".\n";
  
  int node_num = 1001;
  string element_file = "nodes.txt";
  dtable_write ( "nodes.txt", node_num,node_num);
  cout << " Nodes data written to \"" << element_file << "\".\n";
  */
  triangulisation();
}