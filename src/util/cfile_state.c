#include "cfile_state.h"
#include <stdlib.h>
FEMErrorCode FEMFileStateCreateWithFile(FEMCFileState **ppFileState, FILE *pFile, FEMLogLevelType log_level) {
    FEMCFileState *pFileState = malloc(sizeof(FEMCFileState));
    if(pFileState == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "malloc failed in FEMFileStateCreateWithFile\n");
        return FEM_MEM_ERROR;
    }
    if(pFile == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "pFile is NULL in FEMFileStateCreateWithFile\n");
        return FEM_INVALID_INPUT;
    }
    pFileState->fp = pFile;
    pFileState->log_level = log_level;
    pFileState->line_number = 1;
    *ppFileState = pFileState;
    return FEM_SUCCESS;
}

FEMErrorCode FEMCFileStateCreate(FEMCFileState **ppFileState, const char *filename, FEMLogLevelType logLevel) {
    FILE *fp = fopen(filename, "r");
    if(fp == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "fopen failed in FEMCFileStateCreate\n");
        return FEM_FILE_ERROR;
    }
    FEMErrorCode err = FEMFileStateCreateWithFile(ppFileState, fp, logLevel);
    if(err != FEM_SUCCESS) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "FEMFileStateCreateWithFile failed in FEMCFileStateCreate\n");
    }
    return FEM_SUCCESS;
}

FEMErrorCode FEMCFileStateFree(FEMCFileState *pFileState) {
    if(pFileState == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "pFileState is NULL in FEMFileStateFree\n");
        return FEM_INVALID_INPUT;
    }
    if(pFileState->fp != NULL) {
        fclose(pFileState->fp);
    }
    free(pFileState);
    return FEM_SUCCESS;
}

FEMErrorCode FEMCFileStateGetLine(FEMCFileState *pFileState, int maxLineLength, char *line) {
    if(pFileState == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "pFileState is NULL in FEMFileStateGetLine\n");
        return FEM_INVALID_INPUT;
    }
    if(line == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "line is NULL in FEMFileStateGetLine\n");
        return FEM_INVALID_INPUT;
    }
    if(maxLineLength <= 0) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "maxLineLength is not positive in FEMFileStateGetLine\n");
        return FEM_INVALID_INPUT;
    }
    char *ret = fgets(line, maxLineLength, pFileState->fp);
    if(ret == NULL) {
        FEM_LOG(FEM_LOG_LEVEL_ERROR, "fgets gets an EOF or error in FEMFileStateGetLine\n");
        return FEM_FILE_ERROR;
    }
    pFileState->line_number++;
    return FEM_SUCCESS;
}
