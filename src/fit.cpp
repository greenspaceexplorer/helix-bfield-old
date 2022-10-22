#include "helix_magnet.h"
#include <iostream>
#include <algorithm>
#include <gsl/gsl_multimin.h>
#include <omp.h>

using namespace std;

// struct for passing parameters to fit function
struct fit_params{
  helix* magnet;
  vector<AngleAxisd> og_rot; //original rotation 
  vector< vector<double> > og_ir; //original inner radius
  vector< vector<double> > og_or; //original outer radius
  vector< vector<Vector3d> >* data;
  unsigned int dilution;
};
// function to minimize to fit model to data
double fit_func(const gsl_vector* v, void* params);
// outputs vector of field coordinates and values given properly formatted 
//   raw scan data file
vector< vector<Vector3d> > parse_data_csv(const vector<string>& data_csv_vec);
//------------------------------------------------------------------------------
int main(int argc, char *argv[]){
  if(argc < 2){
    throw invalid_argument("FIT_DATA takes 3 or more arguments: \
        magnet_config.csv dilution data1.csv (data2.csv ...)");
  }
  // MAKE HELIX MAGNET OBJECT
  // first input is the magnet configuration
  string magnet_csv = argv[1];
  // instantiate the HELIX magnet
  int contraction_axis = 0;
  helix my_helix(magnet_csv,contraction_axis);
  // print data to check it is correctly loaded
  cout << "Loaded magnet data from " << magnet_csv << " ..." << endl;
  my_helix.print_magnet_info();

  // second input is the dilution of the scan data, to be used later
  int dilution = stoi(argv[2]);
  // LOAD MAGNET SCAN DATA
  // third and subsequent inputs are magnet scan files
  vector<string> data_csv_vec;
  for(int i = 3; i < argc; i++){
    data_csv_vec.push_back(argv[i]);
  }
  // make a vector of coordinates and field measurements
  vector< vector<Vector3d> > data_vec = parse_data_csv(data_csv_vec);

  // COIL FITTING ROUTINE
  // initialize gsl minimizer
  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_vector *step_size;
  gsl_multimin_function min_func;

  size_t iter = 0;
  int status;
  double size;

  // assign input parameters
  fit_params helix_fit_params; // struct containing magnet and data objects
  helix_fit_params.magnet = &my_helix;
  for(int i = 0; i < my_helix.coil_vec.size(); i++){
    helix_fit_params.og_rot.push_back(my_helix.coil_vec[i].get_rotation());
    helix_fit_params.og_ir.push_back(my_helix.coil_vec[i].get_inner_radius());
    helix_fit_params.og_or.push_back(my_helix.coil_vec[i].get_outer_radius());
  }
  helix_fit_params.data = &data_vec;
  helix_fit_params.dilution = dilution; // only evaluate every nth point

//  // assign free parameters and initial values
//  gsl_vector *v; // free parameter vector
//  v = gsl_vector_alloc(8);
//  gsl_vector_set_zero(v); // dtheta,dgamma,y,and z for each coil
//  gsl_vector_set(v,2,my_helix.coil_vec[0].get_origin()(1));
//  gsl_vector_set(v,3,my_helix.coil_vec[0].get_origin()(2));
//  gsl_vector_set(v,6,my_helix.coil_vec[1].get_origin()(1));
//  gsl_vector_set(v,7,my_helix.coil_vec[1].get_origin()(2));
  // assign free parameters and initial values
  gsl_vector *v; // free parameter vector
  v = gsl_vector_alloc(12);
  gsl_vector_set_zero(v); // dtheta,dgamma,x,y,z, and s for each coil
  gsl_vector_set(v,2,my_helix.coil_vec[0].get_origin()(0));
  gsl_vector_set(v,3,my_helix.coil_vec[0].get_origin()(1));
  gsl_vector_set(v,4,my_helix.coil_vec[0].get_origin()(2));
  gsl_vector_set(v,5,1.0);
  gsl_vector_set(v,8,my_helix.coil_vec[1].get_origin()(0));
  gsl_vector_set(v,9,my_helix.coil_vec[1].get_origin()(1));
  gsl_vector_set(v,10,my_helix.coil_vec[1].get_origin()(2));
  gsl_vector_set(v,11,1.0);

//  // assign initial step sizes
//  step_size = gsl_vector_alloc(8);
//  gsl_vector_set_all(step_size,0.001745); // start at about 0.1 deg step size
//  gsl_vector_set(step_size,2,0.001); // change offset in 1 mm increments
//  gsl_vector_set(step_size,3,0.001); 
//  gsl_vector_set(step_size,6,0.001); 
//  gsl_vector_set(step_size,7,0.001); 

  // assign initial step sizes
  step_size = gsl_vector_alloc(12);
  gsl_vector_set_all(step_size,0.001); // start w/ 1 mm shift, ~0.5 deg rotation, .1% scaling

  // initialize minimization method and iterate
  size_t max_iter = 10000; // maximum iterations of minimizer
  double min_size = 1e-4; // minimum simplex size
  //min_func.n = 8;
  min_func.n = 12;
  min_func.f = fit_func;
  min_func.params = &helix_fit_params;

  //s = gsl_multimin_fminimizer_alloc(T,8);
  s = gsl_multimin_fminimizer_alloc(T,12);
  gsl_multimin_fminimizer_set(s, &min_func, v, step_size);

  cout.precision(4);
  cout << scientific;
//  cout << "ITER" << "\t" 
//    << "TH1" << "\t"
//    << "GAM1" << "\t"
//    << "Y1" << "\t"
//    << "Z1" << "\t"
//    << "TH2" << "\t"
//    << "GAM2" << "\t"
//    << "Y2" << "\t"
//    << "Z2" << "\t"
//    << "FUNC" << "\t"
//    << "SIMPLEX" << endl;
  cout << "ITER" << "\t" 
      << "TH1" << "\t"
      << "GAM1" << "\t"
      << "X1" << "\t"
      << "Y1" << "\t"
      << "Z1" << "\t"
      << "S1" << "\t"
      << "TH2" << "\t"
      << "GAM2" << "\t"
      << "X2" << "\t"
      << "Y2" << "\t"
      << "Z2" << "\t"
      << "S2" << "\t"
      << "FUNC" << "\t"
      << "SIMPLEX" << endl;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);

    if(status){ break; }

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size,min_size);

    if(status == GSL_SUCCESS){
        cout << "converged to minimum at" << endl;
//        cout << iter << "\t" 
//            << gsl_vector_get(s->x,0) << "\t"
//            << gsl_vector_get(s->x,1) << "\t"
//            << gsl_vector_get(s->x,2) << "\t"
//            << gsl_vector_get(s->x,3) << "\t"
//            << gsl_vector_get(s->x,4) << "\t"
//            << gsl_vector_get(s->x,5) << "\t"
//            << gsl_vector_get(s->x,6) << "\t"
//            << gsl_vector_get(s->x,7) << "\t"
//            << s->fval << "\t"
//            << size << endl;
        cout << iter << "\t" 
            << gsl_vector_get(s->x,0) << "\t"
            << gsl_vector_get(s->x,1) << "\t"
            << gsl_vector_get(s->x,2) << "\t"
            << gsl_vector_get(s->x,3) << "\t"
            << gsl_vector_get(s->x,4) << "\t"
            << gsl_vector_get(s->x,5) << "\t"
            << gsl_vector_get(s->x,6) << "\t"
            << gsl_vector_get(s->x,7) << "\t"
            << gsl_vector_get(s->x,8) << "\t"
            << gsl_vector_get(s->x,9) << "\t"
            << gsl_vector_get(s->x,10) << "\t"
            << gsl_vector_get(s->x,11) << "\t"
            << s->fval << "\t"
            << size << endl;
    }
    if((iter-1)%10 == 0){
//    cout << iter << "\t" 
//      << gsl_vector_get(s->x,0) << "\t"
//      << gsl_vector_get(s->x,1) << "\t"
//      << gsl_vector_get(s->x,2) << "\t"
//      << gsl_vector_get(s->x,3) << "\t"
//      << gsl_vector_get(s->x,4) << "\t"
//      << gsl_vector_get(s->x,5) << "\t"
//      << gsl_vector_get(s->x,6) << "\t"
//      << gsl_vector_get(s->x,7) << "\t"
//      << s->fval << "\t"
//      << size << endl;
    cout << iter << "\t" 
        << gsl_vector_get(s->x,0) << "\t"
        << gsl_vector_get(s->x,1) << "\t"
        << gsl_vector_get(s->x,2) << "\t"
        << gsl_vector_get(s->x,3) << "\t"
        << gsl_vector_get(s->x,4) << "\t"
        << gsl_vector_get(s->x,5) << "\t"
        << gsl_vector_get(s->x,6) << "\t"
        << gsl_vector_get(s->x,7) << "\t"
        << gsl_vector_get(s->x,8) << "\t"
        << gsl_vector_get(s->x,9) << "\t"
        << gsl_vector_get(s->x,10) << "\t"
        << gsl_vector_get(s->x,11) << "\t"
        << s->fval << "\t"
        << size << endl;
    }
  }while(status == GSL_CONTINUE && iter < max_iter);
  cout << defaultfloat;
  cout.precision(9);

  cout << "HELIX magnet configuration after fitting: " << endl;
  my_helix.decontract();
  my_helix.print_magnet_info();

  // SAVE FITTED MAGNET CONFIGURATION TO FILE
  string config_name = "helix_config_newfit.csv";
  my_helix.save_magnet_info(config_name);
  cout << "Fitted coil output file: " << config_name << endl;

  // free memory
  gsl_vector_free(step_size);
  gsl_vector_free(v);
  gsl_multimin_fminimizer_free(s);

  return 0;
}

//------------------------------------------------------------------------------
// Minimizing this function will give the best fit b/w calculated and measured
//  magnetic fields.
//
//double fit_func(const gsl_vector* v, void* params){
//  // get parameters from void pointer
//  fit_params* par;
//  par = (fit_params*)params;
//  helix* magnet = par->magnet;
//  vector< vector<Vector3d> >* data = par->data;
//
//  // position magnet coils from values in v
//  magnet->coil_vec[0].set_rotation(par->og_rot[0]);
//  magnet->coil_vec[0].compose_rotation(
//      gsl_vector_get(v,0),
//      M_PI/2.,
//      gsl_vector_get(v,1));
//
//  magnet->coil_vec[1].set_rotation(par->og_rot[1]);
//  magnet->coil_vec[1].compose_rotation(
//      gsl_vector_get(v,4),
//      M_PI/2.,
//      gsl_vector_get(v,5));
//  double x0 = magnet->coil_vec[0].get_origin()(0);
//  magnet->coil_vec[0].set_origin(
//      x0,
//      gsl_vector_get(v,2),
//      gsl_vector_get(v,3));
//  double x1 = magnet->coil_vec[1].get_origin()(0);
//  magnet->coil_vec[1].set_origin(
//      x1,
//      gsl_vector_get(v,6),
//      gsl_vector_get(v,7));
//
//  // Sum over the modulus squared of the difference 
//  //  in measured and calculated fields.
//  double sum = 0.0;
//  if(par->dilution <= 1){
//#pragma omp parallel for
//      for(int i = 0; i < (*data)[0].size(); i++){
//          Vector3d diff = magnet->B((*data)[0][i]) - (*data)[1][i];
//          sum += diff.dot(diff);
//      }
//  }
//  else{
//#pragma omp parallel for
//      for(int i = 0; i < (*data)[0].size(); i++){
//          if(i%par->dilution == 0){
//              Vector3d diff = magnet->B((*data)[0][i]) - (*data)[1][i];
//              sum += diff.dot(diff);
//          }
//      }
//  }
//  return sum;
//
//}
double fit_func(const gsl_vector* v, void* params){
    // get parameters from void pointer
    fit_params* par;
    par = (fit_params*)params;
    helix* magnet = par->magnet;
    vector< vector<Vector3d> >* data = par->data;

    // position magnet coils from values in v
//    magnet->coil_vec[0].set_rotation(par->og_rot[0]);
//    magnet->coil_vec[0].compose_rotation(
//            gsl_vector_get(v,0),
//            M_PI/2.,
//            gsl_vector_get(v,1));
//    magnet->coil_vec[0].set_origin(
//            gsl_vector_get(v,2),
//            gsl_vector_get(v,3),
//            gsl_vector_get(v,4));
//
//    magnet->coil_vec[1].set_rotation(par->og_rot[1]);
//    magnet->coil_vec[1].compose_rotation(
//            gsl_vector_get(v,6),
//            M_PI/2.,
//            gsl_vector_get(v,7));
//    magnet->coil_vec[1].set_origin(
//            gsl_vector_get(v,8),
//            gsl_vector_get(v,9),
//            gsl_vector_get(v,10));
    int step = v->size/magnet->coil_vec.size();
    for(int i = 0; i < magnet->coil_vec.size(); i++){
        magnet->coil_vec[i].set_rotation(par->og_rot[i]);
        magnet->coil_vec[i].compose_rotation(
                gsl_vector_get(v,0+i*step),
                M_PI/2.,
                gsl_vector_get(v,1+i*step));
        magnet->coil_vec[i].set_origin(
                gsl_vector_get(v,2+i*step),
                gsl_vector_get(v,3+i*step),
                gsl_vector_get(v,4+i*step));    
        vector<double> irad(par->og_ir[i]),orad(par->og_or[i]);
        for(int j = 0; j < irad.size(); j++){
            irad[j]*=gsl_vector_get(v,5+i*step);
            orad[j]*=gsl_vector_get(v,5+i*step);
        }
        magnet->coil_vec[i].set_inner_radius(irad);
        magnet->coil_vec[i].set_outer_radius(orad);
    }
    // Sum over the modulus squared of the difference 
    //  in measured and calculated fields.
    double sum = 0.0;
#pragma omp parallel for
    for(int i = 0; i < (*data)[0].size(); i++){
        if(i%par->dilution == 0){
            Vector3d diff = magnet->B((*data)[0][i]) - (*data)[1][i];
            sum += diff.dot(diff);
        }
    }
    return sum;

}
//------------------------------------------------------------------------------
// Returns a vector of the form 
//  (vector<Vector3d> coordinates, vector<Vector3d> B_field)
//  from a vector of magnet scan csv files
vector< vector<Vector3d> > parse_data_csv(const vector<string>& data_csv_vec){
  vector<Vector3d> coordinates,bfield;
  vector<vector<Vector3d> > out;
  out.push_back(coordinates);
  out.push_back(bfield);
  int count;
  for(int i = 0; i < data_csv_vec.size(); i++){
    ifstream csv_stream(data_csv_vec[i]);
    if(!csv_stream.good()){
        cout << "Warning: " 
            << data_csv_vec[i] 
            << "is an invalid data file" << endl;
    }
    token_parser csv_row(',');
    count = 0;
    while(csv_stream >> csv_row){
      Vector3d coord;
      Vector3d field;
      coord(0) = stod(csv_row[0]);
      coord(1) = stod(csv_row[1]);
      coord(2) = stod(csv_row[2]);
      // field data appears to be rotated clockwise around z-axis in raw
      //  data files
      field(0) = stod(csv_row[3]);
      field(1) = stod(csv_row[4]);
      field(2) = stod(csv_row[5]);
      out[0].push_back(coord);
      out[1].push_back(field);
      count++;
    }
    cout << "Loaded " << count << " points from " 
      << data_csv_vec[i] << " ..." << endl;
  }
  return out;
}

