#pragma once

#include "../sys/defines.h"

typedef struct {
    FEMBoolType use_tauber;
    FEMBoolType use_sliding_velocity;
    FEMScalarType NS_KDC_W, NS_KDC_A;
    FEMScalarType LSC_KDC;
    FEMScalarType TEM_KDC;
} VMSConfigType;

typedef struct {
    FEMIndexType NRES;
    FEMIndexType *max_iter, *min_iter;
    FEMScalarType *rtol, *atol;
} KSPConfigType;

typedef struct {
    FEMIndexType NRES;
    FEMIndexType *max_iter, *min_iter;
    FEMScalarType *rtol, *atol;
} NewtonRaphsonConfigType;

typedef struct {
    FEMBoolType iga;
    FEMBoolType fem_flag;
    FEMBoolType use_hessian;
    FEMBoolType calc_cfl;
    VMSConfigType *vms;
    KSPConfigType *ksp;
    NewtonRaphsonConfigType *newton_raphson;
} ConfigType;

ConfigType* Config_new(const char *config_file);

#include <stdio.h>
void Config_print(const ConfigType *config, FILE *fp);
void Config_free(ConfigType *config);

void print_config(const ConfigType *config, FILE *fp);