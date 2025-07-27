/*
libgvps
===

Copyright (c) 2015, Kanru Hua
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software without
specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "gvps.h"
#include <stdlib.h>
#include <string.h>

gvps_obsrv* gvps_obsrv_create(int T)
{
    gvps_obsrv* ret = malloc(sizeof(gvps_obsrv));
    ret -> slice = calloc(T, sizeof(gvps_obsrv_slice*));
    ret -> T = T;
    memset(ret -> slice, 0, T * sizeof(gvps_obsrv_slice*));
    return ret;
}

gvps_obsrv_slice* gvps_obsrv_slice_create(int nstate)
{
    gvps_obsrv_slice* ret = malloc(sizeof(gvps_obsrv_slice));
    ret -> pair = calloc(nstate, sizeof(gvps_pair));
    ret -> N = nstate;
    return ret;
}

gvps_obsrv* gvps_obsrv_create_full(int T, int nstate)
{
    int t;
    gvps_obsrv* ret = gvps_obsrv_create(T);
    for(t = 0; t < T; t ++)
        ret -> slice[t] = gvps_obsrv_slice_create(nstate);
    return ret;
}

int gvps_obsrv_slice_free(gvps_obsrv_slice* slice)
{
    free(slice -> pair);
    free(slice);
    return 1;
}

int gvps_obsrv_free(gvps_obsrv* seq)
{
    int t;
    for(t = 0; t < seq -> T; t ++)
        if(seq -> slice[t] != 0)
            gvps_obsrv_slice_free(seq -> slice[t]);
    free(seq -> slice);
    free(seq);
    return 1;
}

