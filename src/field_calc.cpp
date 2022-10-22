#include "helix_magnet.h"
#include "token_parser.h"
#include <iostream>

using namespace std;

vector<Vector3d> parse_coord_csv(const vector<string>& coord_csv_vec);
void field_to_csv(vector< vector<Vector3d> >& data_vec, string filename);
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){
  string ftag="-f",ctag="-c";
  if(argv[2] == ctag){ // individual coordinate input
    // GET COORDINATES
    Vector3d location;
    location(0) = stod(argv[3]);
    location(1) = stod(argv[4]);
    location(2) = stod(argv[5]);

    // MAKE HELIX MAGNET OBJECT
    // first input is the magnet configuration
    string magnet_csv = argv[1];
    // instantiate the HELIX magnet
    helix my_helix(magnet_csv,-1);

    // CALCULATE+OUTPUT
    Vector3d calc = my_helix.B(location);
    cout << calc(0) << "\t" 
      << calc(1) << "\t" 
      << calc(2) << "\t" 
      << calc.norm() << endl;
  }
  else if(argv[2] == ftag){ // coordinate file input
    // GET COORDINATES
    // load coordinates filenames
    vector<string> coord_csv_vec;
    for(int i = 3; i < argc; i++){
      coord_csv_vec.push_back(argv[i]);
    }

    // make a vector of coordinates and field measurements
    vector<Vector3d> coord_vec = parse_coord_csv(coord_csv_vec);
    vector<Vector3d> field_vec(coord_vec.size(),Vector3d());

    // MAKE HELIX MAGNET OBJECT
    // first input is the magnet configuration
    string magnet_csv = argv[1];
    // instantiate the HELIX magnet
    helix my_helix(magnet_csv,-1);

    // CALCULATE FIELD
#pragma omp parallel for
    for(int i = 0; i < coord_vec.size(); i++){
      field_vec[i] = my_helix.B(coord_vec[i]);
      if(i%500 == 0){
        cout << coord_vec.size() - i << " points remaining..." << endl;
      }
    }

    // WRITE TO FILE
    vector< vector<Vector3d> > data_vec;
    data_vec.push_back(coord_vec);
    data_vec.push_back(field_vec);
    string calc_name = "calculated_field.csv";
    field_to_csv(data_vec,calc_name);

    cout << "Saved values to " << calc_name << endl;
  }
  else{
    throw invalid_argument("Invalid input...");
  }


  return 0;
}
//------------------------------------------------------------------------------
// Returns a vector of the form·
//  vector<Vector3d> coordinates from a vector of magnet coordinate files
vector<Vector3d> parse_coord_csv(const vector<string>& coord_csv_vec){
  vector<Vector3d> out;
  int count;
  for(int i = 0; i < coord_csv_vec.size(); i++){
    ifstream csv_stream(coord_csv_vec[i]);
    token_parser csv_row(',');
    count = 0;
    while(csv_stream >> csv_row){
      Vector3d coord;
      // first 3 columns are taken, all others are ignored
      coord(0) = stod(csv_row[0]);
      coord(1) = stod(csv_row[1]);
      coord(2) = stod(csv_row[2]);
      out.push_back(coord);
      count++;
    }
    cout << "Loaded " << count << " points from "
      << coord_csv_vec[i] << "..." << endl;
  }
  return out;
}
//------------------------------------------------------------------------------
// Saves magnetic field data to a csv file
void field_to_csv(vector< vector<Vector3d> >& data_vec, string filename){
  ofstream ofs(filename,ofstream::out);
  ofs.precision(16);
  for(int i = 0; i < data_vec[0].size(); i++){
    for(int j = 0; j < 3; j++){ ofs << data_vec[0][i](j) << ",";}
    for(int k = 0; k < 2; k++){ ofs << data_vec[1][i](k) << ",";}
    ofs << data_vec[1][i](2) << endl;
  }
  ofs.close();
}
