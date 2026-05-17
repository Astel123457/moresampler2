#include "wavtool.h"

int wavtool(wavtool_data* data) {
    
    printf("wavtool: out=%s in=%s offset=%.2f length=%.2f\n",
           data->output, data->input, data->offset, data->length);
    return 0;
}