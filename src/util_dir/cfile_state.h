#pragma once

#include <sys_dir/defines.h>
#include <stdio.h>


struct _p_FEMCFileState {
    FILE *fp;
    int line_number;
    FEMLogLevelType log_level;
};

typedef struct _p_FEMCFileState FEMCFileState;

FEMErrorCode FEMCFileStateCreate(FEMCFileState **state, const char *filename, FEMLogLevelType log_level);
FEMErrorCode FEMCFileStateCreateWithFile(FEMCFileState **state, FILE *fp, FEMLogLevelType log_level);
FEMErrorCode FEMCFileStateFree(FEMCFileState *state);
FEMErrorCode FEMCFileStateGetLine(FEMCFileState *state, int buffer_size, char *buffer);


