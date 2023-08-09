#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAX_LINE_LENGTH 1024

#define TOLOWERCASE(x) \
do { \
    for (char *p = x; *p != '\0'; p++) { \
        *p = tolower(*p); \
    } \
} while (0)

/*
 * Find the value of a key in a config file
 * The files are assumed to be of the form
    * # comment
    * key = value
    * ...
 */
static FEMErrorCode find_config_value_str(FILE *fp, const char *key, char *value) {
    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '#') continue;
        TOLOWERCASE(line);
        if (strstr(line, key) != NULL) {
            char *token = strtok(line, " =");
            token = strtok(NULL, " =");
            strcpy(value, token);
            return FEM_SUCCESS;
        }
    }
    *value = '\0';
    fprintf(stderr, "Error: could not find key %s in config file. "
            "Will use the default value: %s\n", key, value);
    return FEM_INVALID_INPUT;
}

static FEMErrorCode find_config_value_int(FILE *fp, const char *key, int *value) {
    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '#') continue;
        TOLOWERCASE(line);
        if (strstr(line, key) != NULL) {
            char *token = strtok(line, " =");
            token = strtok(NULL, " =");
            *value = atoi(token);
            return FEM_SUCCESS;
        }
    }
    *value = 0;
    fprintf(stderr, "Error: could not find key %s in config file. "
            "Will use the default value: %d\n", key, *value);
    return FEM_INVALID_INPUT;
}

static FEMErrorCode find_config_value_double(FILE *fp, const char *key, double* value) {
    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '#') continue;
        TOLOWERCASE(line);
        if (strstr(line, key) != NULL) {
            char *token = strtok(line, " =");
            token = strtok(NULL, " =");
            *value = atof(token);
            return FEM_SUCCESS;
        }
    }
    *value = 0.0;
    fprintf(stderr, "Error: could not find key %s in config file. "
            "Will use the default value: %g\n", key, *value);
    return FEM_INVALID_INPUT;
}

static FEMErrorCode find_config_value_bool(FILE *fp, const char *key, bool *value) {
    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (line[0] == '#') continue;
        if (strstr(line, key) != NULL) {
            char *token = strtok(line, " =");
            token = strtok(NULL, " =");
            if (strcmp(token, "true") == 0) {
                *value = true;
            } else if (strcmp(token, "false") == 0) {
                *value = false;
            } else {
                fprintf(stderr, "Error: could not parse boolean value %s\n", token);
                exit(FEM_INVALID_INPUT);
            }
            return FEM_SUCCESS;
        }
    }
    value = false;
    fprintf(stderr, "Error: could not find key %s in config file. "
            "Will use the default value: %s\n", key, *value ? "true" : "false");
    return FEM_INVALID_INPUT;
}

ConfigType* Config_new(const char *config_file) {

    FILE *fp = fopen(config_file, "r");
    printf("config_file = %s\n", config_file);
    if (fp == NULL) {
        fprintf(stderr, "Error: could not open config file %s\n", config_file);
        exit(FEM_INVALID_INPUT);
    }
    ConfigType *config = (ConfigType *)malloc(sizeof(ConfigType));

    if (config == NULL) {
        fprintf(stderr, "Error: could not allocate memory for config\n");
        exit(FEM_RUNTIME_ERROR);
    }
    FEMErrorCode err;
    FEM_CHECK_ERROR(find_config_value_bool(fp, "iga", &config->iga));
    FEM_CHECK_ERROR(find_config_value_bool(fp, "fem_flag", &config->fem_flag));
    FEM_CHECK_ERROR(find_config_value_bool(fp, "use_hessian", &config->use_hessian));
    FEM_CHECK_ERROR(find_config_value_bool(fp, "calc_cfl", &config->calc_cfl));


    int nRES;
    FEM_CHECK_ERROR(find_config_value_int(fp, "nres", &nRES));

    config->vms = (VMSConfigType *)malloc(sizeof(VMSConfigType));
    if (config->vms == NULL) {
        fprintf(stderr, "Error: could not allocate memory for vms config\n");
        exit(FEM_RUNTIME_ERROR);
    }
    FEM_CHECK_ERROR(find_config_value_bool(fp, "use_tauber", &config->vms->use_tauber));
    FEM_CHECK_ERROR(find_config_value_bool(fp, "use_sliding_velocity", &config->vms->use_sliding_velocity));
    FEM_CHECK_ERROR(find_config_value_double(fp, "ns_kdc_w", &config->vms->NS_KDC_W));
    FEM_CHECK_ERROR(find_config_value_double(fp, "ns_kdc_a", &config->vms->NS_KDC_A));
    FEM_CHECK_ERROR(find_config_value_double(fp, "lsc_kdc", &config->vms->LSC_KDC));
    FEM_CHECK_ERROR(find_config_value_double(fp, "tem_kdc", &config->vms->TEM_KDC));

    config->ksp = (KSPConfigType *)malloc(sizeof(KSPConfigType));
    if (config->ksp == NULL) {
        fprintf(stderr, "Error: could not allocate memory for ksp config\n");
        exit(FEM_RUNTIME_ERROR);
    }
    config->ksp->NRES = nRES;
    config->ksp->max_iter = (FEMIndexType *)malloc(nRES * sizeof(FEMIndexType));
    config->ksp->min_iter = (FEMIndexType *)malloc(nRES * sizeof(FEMIndexType));
    config->ksp->rtol = (FEMScalarType *)malloc(nRES * sizeof(FEMScalarType));
    config->ksp->atol = (FEMScalarType *)malloc(nRES * sizeof(FEMScalarType));
    if (config->ksp->max_iter == NULL || config->ksp->min_iter == NULL || config->ksp->rtol == NULL || config->ksp->atol == NULL) {
        fprintf(stderr, "Error: could not allocate memory for ksp config\n");
        exit(FEM_RUNTIME_ERROR);
    }
    for (int i = 0; i < nRES; i++) {
        char key[32];
        sprintf(key, "ksp_max_iter_%d", i);
        FEM_CHECK_ERROR(find_config_value_int(fp, key, &config->ksp->max_iter[i]));
        sprintf(key, "ksp_min_iter_%d", i);
        FEM_CHECK_ERROR(find_config_value_int(fp, key, &config->ksp->min_iter[i]));
        sprintf(key, "ksp_rtol_%d", i);
        FEM_CHECK_ERROR(find_config_value_double(fp, key, &config->ksp->rtol[i]));
        sprintf(key, "ksp_atol_%d", i);
        FEM_CHECK_ERROR(find_config_value_double(fp, key, &config->ksp->atol[i]));
    }

    config->newton_raphson = (NewtonRaphsonConfigType *)malloc(sizeof(NewtonRaphsonConfigType));
    if (config->newton_raphson == NULL) {
        fprintf(stderr, "Error: could not allocate memory for newton_raphson config\n");
        exit(FEM_RUNTIME_ERROR);
    }
    config->newton_raphson->NRES = nRES;
    config->newton_raphson->max_iter = (FEMIndexType *)malloc(1 * sizeof(FEMIndexType));
    config->newton_raphson->min_iter = (FEMIndexType *)malloc(1 * sizeof(FEMIndexType));
    config->newton_raphson->rtol = (FEMScalarType *)malloc(nRES * sizeof(FEMScalarType));
    config->newton_raphson->atol = (FEMScalarType *)malloc(nRES * sizeof(FEMScalarType));

    if (config->newton_raphson->max_iter == NULL || config->newton_raphson->min_iter == NULL || config->newton_raphson->rtol == NULL || config->newton_raphson->atol == NULL) {
        fprintf(stderr, "Error: could not allocate memory for newton_raphson config\n");
        exit(FEM_RUNTIME_ERROR);
    }
    FEM_CHECK_ERROR(find_config_value_int(fp, "newton_raphson_max_iter", &config->newton_raphson->max_iter[0]));
    FEM_CHECK_ERROR(find_config_value_int(fp, "newton_raphson_min_iter", &config->newton_raphson->min_iter[0]));
    for(int i = 0; i < nRES; i++) {
        char key[32];
        sprintf(key, "newton_raphson_rtol_%d", i);
        FEM_CHECK_ERROR(find_config_value_double(fp, key, &config->newton_raphson->rtol[i]));
        sprintf(key, "newton_raphson_atol_%d", i);
        FEM_CHECK_ERROR(find_config_value_double(fp, key, &config->newton_raphson->atol[i]));
    }
}

void Config_print(const ConfigType *config, FILE *fp) {
    if(fp == NULL) {
        fp = stdout;
    }
    fprintf(fp, "config.iga = %s\n", config->iga ? "true" : "false");
    fprintf(fp, "config.fem = %s\n", config->fem_flag ? "true" : "false");
    fprintf(fp, "config.use_hessian = %s\n", config->use_hessian ? "true" : "false");
    fprintf(fp, "config.calc_cfl = %s\n", config->calc_cfl ? "true" : "false");
    fprintf(fp, "config.vms.use_tauber = %s\n", config->vms->use_tauber ? "true" : "false");
    fprintf(fp, "config.vms.use_sliding_velocity = %s\n", config->vms->use_sliding_velocity ? "true" : "false");
    fprintf(fp, "config.vms.NS_KDC_W = %g\n", config->vms->NS_KDC_W);
    fprintf(fp, "config.vms.NS_KDC_A = %g\n", config->vms->NS_KDC_A);
    fprintf(fp, "config.vms.LSC_KDC = %g\n", config->vms->LSC_KDC);
    fprintf(fp, "config.vms.TEM_KDC = %g\n", config->vms->TEM_KDC);
    fprintf(fp, "config.ksp.NRES = %d\n", config->ksp->NRES);

    fprintf(fp, "config.ksp.max_iter[:] = ");
    for (int i = 0; i < config->ksp->NRES; i++) {
        fprintf(fp, "%d", config->ksp->max_iter[i]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "config.ksp.min_iter[:] = ");
    for (int i = 0; i < config->ksp->NRES; i++) {
        fprintf(fp, "%d", config->ksp->min_iter[i]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "config.ksp.rtol[:] = ");
    for (int i = 0; i < config->ksp->NRES; i++) {
        fprintf(fp, "%g", config->ksp->rtol[i]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "config.ksp.atol[:] = ");
    for (int i = 0; i < config->ksp->NRES; i++) {
        fprintf(fp, "%g", config->ksp->atol[i]);
    }
    fprintf(fp, "\n");

    fprintf(fp, "config.newton_raphson.NRES = %d\n", config->newton_raphson->NRES);
    fprintf(fp, "config.newton_raphson.max_iter[:] = %d\n",
            config->newton_raphson->max_iter[0]);
    fprintf(fp, "config.newton_raphson.min_iter[:] = %d\n",
            config->newton_raphson->min_iter[0]);
    fprintf(fp, "config.newton_raphson.rtol[:] = ");
    for (int i = 0; i < config->newton_raphson->NRES; i++) {
        fprintf(fp, "%g", config->newton_raphson->rtol[i]);
    }
    fprintf(fp, "\n");
    fprintf(fp, "config.newton_raphson.atol[:] = ");
    for (int i = 0; i < config->newton_raphson->NRES; i++) {
        fprintf(fp, "%g", config->newton_raphson->atol[i]);
    }
    fprintf(fp, "\n");
}

void Config_free(ConfigType *config) {
    if(config == NULL)
        return;
    if(config->vms) {
        free(config->vms);
        config->vms = NULL;
    }
    if(config->ksp) {
        free(config->ksp->max_iter);
        free(config->ksp->min_iter);
        free(config->ksp->rtol);
        free(config->ksp->atol);
        free(config->ksp);
        config->ksp = NULL;
    }
    if(config->newton_raphson) {
        free(config->newton_raphson->max_iter);
        free(config->newton_raphson->min_iter);
        free(config->newton_raphson->rtol);
        free(config->newton_raphson->atol);
        free(config->newton_raphson);
    }
    free(config);
}

void print_config(const ConfigType*config, FILE*fp) {}