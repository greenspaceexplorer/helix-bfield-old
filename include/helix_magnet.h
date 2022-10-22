#ifndef HELIX_MAGNET_H
#define HELIX_MAGNET_H

#include "coil.h"
#include "token_parser.h"
#include <fstream>

class helix{
    public:

        // constructor from coil objects
        helix(vector<coil> helix_coils, int caxis);

        // constructor from magnet configuration file
        helix(const string magnet_config, int caxis);

        // destructor
        ~helix();

        // sets thermal contraction axis
        //  the stationary coil is the last one in the helix_coils vector
        void set_con_axis(int caxis);

        // gets thermal contraction axis
        int get_con_axis();

        // applies thermal contraction to coil radii and positions
        void contract();

        // reverses thermal contraction
        void decontract();

        // returns B field at given position
        Vector3d B(Vector3d &position);

        // prints stored magnet information
        void print_magnet_info();

        // saves magnet information to csv file
        void save_magnet_info(string filename);

        // outputs vector of coil objects given properly formatted magnet
        //   configuration file
        vector<coil> parse_magnet_csv(const string csv_file);

        // vector containing helix coil objects
        vector<coil> coil_vec;

    private:

        // magnitude of thermal contraction for copper
        double cont_cu = 0.99674;

        // magnitude of thermal contraction for aluminum
        double cont_al = 0.99585;

        // axis along which thermal contraction takes place
        int con_axis;

        // boolean to keep track if contraction has been applied
        bool contraction;
        // boolean denoting if default contraction is desired
        bool default_contraction;

        // variables for reversing thermal contraction when saving data
        vector<double> z_shift;
};
#endif // HELIX_MAGNET_H
