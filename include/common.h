#ifndef COMMON_H
#define COMMON_H

extern "C"{
void heat_mag_(
        double *xin, 
        double *yin, 
        double *zin, 
        double *brho, 
        double *bz, 
        double *b
        );
void excoil_(
        double *rhoout,
        double *rhoin,
        double *w,
        double *i,
        double *nturns,
        int *nrho,
        int *nz,
        double *x,
        double *y,
        double *z,
        double *brhoec,
        double *bzec,
        double *modbec
        );
void coilb_(
        double *i,
        double *a,
        double *nturns,
        double *x,
        double *y,
        double *z,
        double *brho,
        double *bz,
        double *modb
        );
double cel(
        double *qqc,
        double *pp,
        double *aa,
        double *bb
        );
}

#endif
