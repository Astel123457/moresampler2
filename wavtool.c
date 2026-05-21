#include "wavtool.h"
#include "utils.h"
#include <stdio.h>
#include <libllsm/llsm.h>

int wavtool(wavtool_data* data) {
    
    int nfrm_in = 0, fs_in = 0, nbit_in = 0, nx_in = 0, nhop_in = 128;
    llsm_aoptions* opt_a_in = llsm_create_aoptions();
    llsm_soptions* opt_s_in = llsm_create_soptions(44100.0f);
    llsm_chunk* input = get_chunk_from_file(data->input, &nfrm_in, &fs_in, &nbit_in, &nx_in, &nhop_in, 1, opt_a_in, opt_s_in);
    if (!input) {
        fprintf(stderr, "Failed to load input file: %s\n", data->input);
        llsm_delete_aoptions(opt_a_in);
        llsm_delete_soptions(opt_s_in);
        return 1;
    }
    int nfrm_out = 0, fs_out = 0, nbit_out = 0, nx_out = 0, nhop_out = 128;
    llsm_aoptions* opt_a_out = llsm_create_aoptions();
    llsm_soptions* opt_s_out = llsm_create_soptions(44100.0f);
    llsm_chunk* output = get_chunk_from_file(data->output, &nfrm_out, &fs_out, &nbit_out, &nx_out, &nhop_out, 1, opt_a_out, opt_s_out);
    if (!output) {
        fprintf(stderr, "Failed to load output file: %s\n", data->output);
        llsm_delete_chunk(input);
        llsm_delete_aoptions(opt_a_in);
        llsm_delete_soptions(opt_s_in);
        llsm_delete_aoptions(opt_a_out);
        llsm_delete_soptions(opt_s_out);
        return 1;
    }

    
    
    return 0;
}