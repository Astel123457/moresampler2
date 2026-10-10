// Inspect phrase model/envelope state without adding public test-only API.
#include "../wavtool.c"

static void check(int ok, const char* message) {
    if (!ok) { fprintf(stderr, "%s\n", message); exit(1); }
}

static llsm_chunk* source(void) {
    llsm_aoptions* options = llsm_create_aoptions();
    llsm_container* conf = llsm_aoptions_toconf(options, 22050);
    llsm_container_attach(conf, LLSM_CONF_NFRM, llsm_create_int(20),
        llsm_delete_int, llsm_copy_int);
    llsm_chunk* chunk = llsm_create_chunk(conf, 1);
    llsm_delete_container(conf);
    llsm_delete_aoptions(options);
    for (int i = 0; i < 20; ++i) {
        llsm_container* frame = chunk->frames[i];
        llsm_container_attach(frame, LLSM_FRAME_F0, llsm_create_fp(200), llsm_delete_fp, llsm_copy_fp);
        llsm_container_attach(frame, LLSM_FRAME_RD, llsm_create_fp(1), llsm_delete_fp, llsm_copy_fp);
        FP_TYPE* vt = llsm_create_fparray(4);
        for (int j = 0; j < 4; ++j) vt[j] = -10;
        llsm_container_attach(frame, LLSM_FRAME_VTMAGN, vt, llsm_delete_fparray, llsm_copy_fparray);
        llsm_container_attach(frame, LLSM_FRAME_VSPHSE, llsm_create_fparray(1), llsm_delete_fparray, llsm_copy_fparray);
    }
    return chunk;
}

int main(void) {
    llsm_chunk* input = source();
    float hop = 128000.0f / 44100;
    // Repeated phrase construction also exercises frame ownership/teardown.
    for (int pass = 0; pass < 32; ++pass) {
        int level = pass % 2 ? 100 : 50;
        wavtool_phrase* phrase = wavtool_phrase_create();
        wavtool_data a = {0};
        a.length = 20 * hop;
        a.p3 = 10 * hop;
        a.v1 = a.v2 = a.v3 = level;
        check(!wavtool_phrase_add(phrase, input, 44100, 16, 0, &a), "first add failed");
        wavtool_data b = {0};
        b.length = 20 * hop;
        b.position = b.ovr = b.p2 = 10 * hop;
        b.v2 = b.v3 = b.v4 = level;
        check(!wavtool_phrase_add(phrase, input, 44100, 16, 0, &b), "second add failed");
        for (int i = 0; i < 30; ++i) {
            check(fabsf(phrase->envelope[i] - level) < .001f, "complementary fades must sum to the intended level");
            FP_TYPE* vt = llsm_container_get(phrase->chunk->frames[i], LLSM_FRAME_VTMAGN);
            check(fabs(vt[0] + 10) < .001, "stitch must not attenuate the raw spectrum");
        }
        wavtool_phrase_delete(phrase);
    }
    llsm_delete_chunk(input);
    return 0;
}
