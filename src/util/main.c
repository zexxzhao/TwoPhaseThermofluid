
#include "config.h"

int main(int argc, char **argv) {
    printf("Hello, world!\n");
    ConfigType *config = Config_new("config.dat");
    // Config_print(config, stdout);
    Config_free(config);
    return 0;
}