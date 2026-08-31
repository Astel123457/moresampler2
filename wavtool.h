#include <libllsm/llsm.h>

// Envelope points: p = position in ms into the sample, v = volume (0-100).
// p1/p2 are measured forward from the start; p3/p4 backward from the end.
// p5 is an optional middle point measured forward (0 when unused).
typedef struct {
    char* output; // Path to the output wav file
    char* input;  // Path to the resampled input wav file
    float offset; // STP offset in milliseconds (start trim)
    float length; // Note length in milliseconds
    float p1;     // Position of envelope point 1 (ms, forward)
    float p2;     // Position of envelope point 2 (ms, forward)
    float p3;     // Position of envelope point 3 (ms, backward) — next-overlap start
    float p4;     // Position of envelope point 4 (ms, backward) — next-overlap end
    float p5;     // Position of optional middle envelope point (ms, forward; 0 if unused)
    float v1;     // Volume at p1
    float v2;     // Volume at p2
    float v3;     // Volume at p3
    float v4;     // Volume at p4
    float v5;     // Volume at p5 (100 if unused)
    float ovr;    // Overlap in milliseconds
    float delta;  // Trailing +/- correction parsed from ticks@tempo+delta
    // Phrase-renderer extensions. Ignored by the traditional CLI path.
    float position;       // Absolute phrase-relative position in milliseconds.
    int use_position;
    int first;
    int defer_synthesis;
} wavtool_data;

int wavtool(wavtool_data* data);

typedef struct wavtool_phrase wavtool_phrase;
wavtool_phrase* wavtool_phrase_create(void);
int wavtool_phrase_add(wavtool_phrase* phrase, const llsm_chunk* input,
    int fs, int nbit, float gain_db, const wavtool_data* data);
int wavtool_phrase_synthesize(wavtool_phrase* phrase, float** samples,
    int* sample_count, int* sample_rate);
void wavtool_phrase_delete(wavtool_phrase* phrase);
