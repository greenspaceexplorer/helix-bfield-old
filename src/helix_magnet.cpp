#include "helix_magnet.h"

//------------------------------------------------------------------------------ 
helix::helix(vector<coil> helix_coils, int caxis)
{
    contraction = false;
    default_contraction = false;
    coil_vec = helix_coils;
    set_con_axis(caxis); // use caxis = -1 for axis with largest coil displacement
}
//------------------------------------------------------------------------------ 
helix::helix(const string magnet_config, int caxis)
{
    contraction = false;
    default_contraction = false;
    coil_vec = parse_magnet_csv(magnet_config);
    set_con_axis(caxis); // use caxis = -1 for axis with largest coil displacement
    if(default_contraction){
        contract();
    }
}
//------------------------------------------------------------------------------ 
helix::~helix(){};
//------------------------------------------------------------------------------ 
void helix::set_con_axis(int caxis){
    if(caxis == -1){
        if(coil_vec.size() <= 1){
            con_axis = 2;
        }
        else{
            vector<Vector3d> diff_vec;
            for(int i = 0; i < coil_vec.size(); i++){
                for(int j = i; j < coil_vec.size(); j++){
                    diff_vec.push_back(
                            coil_vec[i].get_origin()-coil_vec[j].get_origin()
                            );
                }
            }
            double max_delta = 0.;
            double new_delta;
            con_axis = 2;
            for(int i = 0; i < diff_vec.size(); i++){
                for(int j = 0; j < 3; j++){
                    new_delta = abs(diff_vec[i](j));
                    if(new_delta > max_delta){
                        max_delta = new_delta;
                        con_axis = j;
                    }
                }
            }
        }
    }
    else if(caxis < 0 || caxis > 2){
        cout << "Error: invalid contraction axis number" << endl;
        exit(1);
    }
    else{
        con_axis = caxis;
    }
}
//------------------------------------------------------------------------------ 
int helix::get_con_axis(){ return con_axis; }
//------------------------------------------------------------------------------ 
void helix::contract(){
    if(contraction){
        cout << "Thermal contraction is already applied to coils." << endl;
    }
    else{
        //cout << "-----Contracted coils-----" << endl;
        contraction = true;
        // clear z_shift vector in case previous contractions have been applied
        z_shift.clear();
        // get z coordinate of the base (most positive along contraction axis) coil
        double z_base = coil_vec[0].get_origin()(con_axis);
        double base_check;
        for(int i = 0; i < coil_vec.size(); i++){
            base_check = coil_vec[i].get_origin()(con_axis);
            if(base_check > z_base){
                z_base = base_check;
            }
        }

        // apply thermal contraction
        for(int i = 0; i < coil_vec.size(); i++){
            // apply copper contraction to coil width
            double coil_w = coil_vec[i].get_width();
            coil_vec[i].set_width(coil_w*cont_cu);

            // apply aluminum contraction to coil separation
            Vector3d coil_origin = coil_vec[i].get_origin();
            double zshift = cont_al*(z_base-coil_origin(con_axis));
            coil_origin(con_axis) = z_base-zshift;
            coil_vec[i].set_origin(coil_origin);

            // apply copper contraction to subcoil radii
            vector<double> coil_ir = coil_vec[i].get_inner_radius();
            vector<double> coil_or = coil_vec[i].get_outer_radius();
            for(int j = 0; j < coil_ir.size(); j++){
                coil_ir[j]*=cont_cu;
                coil_or[j]*=cont_cu;
            } 
            coil_vec[i].set_inner_radius(coil_ir);
            coil_vec[i].set_outer_radius(coil_or); 
        }  
    }
}
//------------------------------------------------------------------------------ 
void helix::decontract(){
    if(!contraction){
        cout << "Thermal contraction is not applied to coils." << endl;
    }
    else{
        //cout << "-----Decontracted coils-----" << endl;
        contraction = false;
        // get z coordinate of the base (most positive along contraction axis) coil
        double z_base = coil_vec[0].get_origin()(con_axis);
        double base_check;
        for(int i = 0; i < coil_vec.size(); i++){
            base_check = coil_vec[i].get_origin()(con_axis);
            if(base_check > z_base){
                z_base = base_check;
            }
        }
        for(int i = 0; i < coil_vec.size(); i++){
            // decontract coil width
            double coil_w = coil_vec[i].get_width();
            coil_vec[i].set_width(coil_w/cont_cu);

            // decontract coil separation
            Vector3d coil_origin = coil_vec[i].get_origin();
            double zshift = (z_base-coil_origin(con_axis))/cont_al;
            coil_origin(con_axis) = z_base - zshift;
            coil_vec[i].set_origin(coil_origin);

            // decontract subcoil radii
            vector<double> coil_ir = coil_vec[i].get_inner_radius();
            vector<double> coil_or = coil_vec[i].get_outer_radius();
            for(int j = 0; j < coil_ir.size(); j++){
                coil_ir[j]/=cont_cu;
                coil_or[j]/=cont_cu;
            } 
            coil_vec[i].set_inner_radius(coil_ir);
            coil_vec[i].set_outer_radius(coil_or); 
        }

    }
}
//------------------------------------------------------------------------------ 
Vector3d helix::B(Vector3d &position){
    Vector3d out(0.,0.,0.);
    for(int i=0; i < coil_vec.size(); i++){
        out += coil_vec[i].B(position);
    }
    return out;
}
//------------------------------------------------------------------------------ 
void helix::print_magnet_info(){
    if(contraction){
    cout << "Note: thermal contraction applied to coils." << endl;
    }
    for(int i = 0; i < coil_vec.size(); i++){
        cout << "------Coil " << i << "------" << endl;
        coil_vec[i].print_coil_info();
    }
}
//------------------------------------------------------------------------------ 
void helix::save_magnet_info(string filename){

    if(contraction){ decontract(); }

    ofstream ofs(filename,ofstream::out);
    ofs.precision(16);
    AngleAxisd coil_angleaxis;
    Vector3d coil_origin;
    vector<double> coil_ir,coil_or,coil_turns;
    vector<unsigned int> coil_rho_div,coil_z_div;

    ofs << "default_contraction," << default_contraction << endl;
    for(int i = 0; i < coil_vec.size(); i++){
        coil_angleaxis = coil_vec[i].get_rotation();
        coil_origin    = coil_vec[i].get_origin();
        coil_ir        = coil_vec[i].get_inner_radius();
        coil_or        = coil_vec[i].get_outer_radius();
        coil_turns     = coil_vec[i].get_subcoil_turns();
        coil_rho_div   = coil_vec[i].get_rho_elements();
        coil_z_div     = coil_vec[i].get_z_elements();

        ofs << "new_coil" << endl;
        ofs << "current," << coil_vec[i].get_current() << endl;
        ofs << "width,"   << coil_vec[i].get_width() << endl;
        ofs << "rotation_angleaxis," 
            << coil_angleaxis.angle() << ","
            << coil_angleaxis.axis()(0) << ","
            << coil_angleaxis.axis()(1) << ","
            << coil_angleaxis.axis()(2) << endl;
        ofs << "origin,"
            << coil_origin(0) << ","
            << coil_origin(1) << ","
            << coil_origin(2) << endl;
        ofs << "inner_radius,";
        for(int j = 0; j < coil_ir.size(); j++){
            if(j == coil_ir.size()-1){
                ofs << coil_ir[j] << endl;
            }
            else{
                ofs << coil_ir[j] << ",";
            }
        }
        ofs << "outer_radius,";
        for(int j = 0; j < coil_or.size(); j++){
            if(j == coil_or.size()-1){
                ofs << coil_or[j] << endl;
            }
            else{
                ofs << coil_or[j] << ",";
            }
        }
        ofs << "subcoil_turns,";
        for(int j = 0; j < coil_turns.size(); j++){
            if(j == coil_turns.size()-1){
                ofs << coil_turns[j] << endl;
            }
            else{
                ofs << coil_turns[j] << ",";
            }
        }

        ofs << "subcoil_rho_div,";
        for(int j = 0; j < coil_rho_div.size(); j++){
            if(j == coil_rho_div.size()-1){
                ofs << coil_rho_div[j] << endl;
            }
            else{
                ofs << coil_rho_div[j] << ",";
            }
        }
        ofs << "subcoil_z_div,";
        for(int j = 0; j < coil_z_div.size(); j++){
            if(j == coil_z_div.size()-1){
                ofs << coil_z_div[j] << endl;
            }
            else{
                ofs << coil_z_div[j] << ",";
            }
        }
    }
    ofs.close();
    cout << "HELIX magnet configuration saved to " << filename << endl;
}
//------------------------------------------------------------------------------
// Returns a vector of coil objects that can instantiate a helix object given
//  a properly formatted magnet configuration csv.
vector<coil> helix::parse_magnet_csv(const string csv_file){
    ifstream csv_stream(csv_file);
    if(!csv_stream.good()){
        cout << "Error: invalid magnet config file" << endl;
        exit(1);
    }
    vector<coil> out;
    token_parser csv_row(',');
    int count = -1;
    while(csv_stream >> csv_row){
        if(csv_row[0] == "default_contraction"){
            default_contraction = stoi(csv_row[1]);
        }
        else if(csv_row[0] == "new_coil"){
            count++;
            coil new_coil = coil();
            out.push_back(new_coil);
        }
        else if(csv_row[0] == "current"){
            double current = stod(csv_row[1]);
            out[count].set_current(current);
        }
        else if(csv_row[0] == "width"){
            double width = stod(csv_row[1]);
            out[count].set_width(width);
        }
        else if(csv_row[0] == "rotation_angleaxis"){
            double rotation_angle = stod(csv_row[1]);
            Vector3d rotation_axis;
            for(int i = 2; i < csv_row.size(); i++){
                rotation_axis(i-2) = stod(csv_row[i]);
            }
            out[count].set_rotation(rotation_angle,rotation_axis);
        }
        else if (csv_row[0] == "origin"){
            Vector3d origin;
            for(int i = 1; i < csv_row.size(); i++){
                origin(i-1) = stod(csv_row[i]);
            }
            out[count].set_origin(origin);
        }
        else if(csv_row[0] == "inner_radius"){
            vector<double> inner_radius;
            for(int i = 1; i < csv_row.size(); i++){
                double ir_d = stod(csv_row[i]);
                inner_radius.push_back(ir_d);
            }
            out[count].set_inner_radius(inner_radius);
        }
        else if(csv_row[0] == "outer_radius"){
            vector<double> outer_radius;
            for(int i = 1; i < csv_row.size(); i++){
                outer_radius.push_back(stod(csv_row[i]));
            }
            out[count].set_outer_radius(outer_radius);
        }
        else if(csv_row[0] == "subcoil_turns"){
            vector<double> subcoil_turns;
            for(int i = 1; i < csv_row.size(); i++){
                subcoil_turns.push_back(stod(csv_row[i]));
            }
            out[count].set_subcoil_turns(subcoil_turns);
        }
        else if(csv_row[0] == "subcoil_rho_div"){
            vector<unsigned int> subcoil_rho_div;
            for(int i = 1; i < csv_row.size(); i++){
                subcoil_rho_div.push_back(stoi(csv_row[i]));
            }
            out[count].set_rho_elements(subcoil_rho_div);
        }
        else if(csv_row[0] == "subcoil_z_div"){
            vector<unsigned int> subcoil_z_div;
            for(int i = 1; i < csv_row.size(); i++){
                subcoil_z_div.push_back(stoi(csv_row[i]));
            }
            out[count].set_z_elements(subcoil_z_div);
        }
        else{
            throw "Error: invalid row encountered while parsing magnet data.";
        }
    }
    return out; 
}
