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

#ifndef GVPS
#define GVPS

typedef struct
{
    FP_TYPE p;
    int state;
    void* ptr;
} gvps_pair;

typedef struct
{
    gvps_pair* pair;
    int N;
} gvps_obsrv_slice;

typedef struct
{
    gvps_obsrv_slice** slice;
    int T;
} gvps_obsrv;

gvps_obsrv* gvps_obsrv_create(int T);
gvps_obsrv* gvps_obsrv_create_full(int T, int nstate);
gvps_obsrv_slice* gvps_obsrv_slice_create(int nstate);
int gvps_obsrv_slice_free(gvps_obsrv_slice* slice);
int gvps_obsrv_free(gvps_obsrv* seq);

typedef FP_TYPE (* gvps_fobsrv)(void* context, int t, int n);
typedef int (* gvps_fnobsrv)(void* context, int t);
typedef FP_TYPE (* gvps_ftran_sampled)(void* context, int ds, int t);
typedef FP_TYPE (* gvps_ftran)(void* context, int s1, int s2, int t2);
typedef int (* gvps_ftrexist)(void* context, int s1, int s2, int t2);
typedef int (* gvps_fntran)(void* context, int t);

#define GVPS_ARG_SPSAMPHIDDEN (int* dst, void* context, int nstate, gvps_obsrv* obsrv, \
    gvps_ftran_sampled ftransame, gvps_ftran_sampled ftrandiff, gvps_fntran fntran, int nhiddenprune)

#define GVPS_ARG_SPSAMP (int* dst, void* context, int nstate, gvps_obsrv* obsrv, \
    gvps_ftran_sampled ftran, gvps_fntran fntran)

#define GVPS_ARG_FULLSAMPHIDDEN (int* dst, void* context, gvps_obsrv* obsrv, \
    gvps_ftran_sampled ftransame, gvps_ftran_sampled ftrandiff, gvps_fntran fntran, int nhiddenprune)

#define GVPS_ARG_FULLSAMP (int* dst, void* context, gvps_obsrv* obsrv, \
    gvps_ftran_sampled ftran, gvps_fntran fntran)

FP_TYPE gvps_sparse_sampled_hidden GVPS_ARG_SPSAMPHIDDEN;
FP_TYPE gvps_sparse_sampled GVPS_ARG_SPSAMP;

FP_TYPE gvps_sparse_circular_hidden GVPS_ARG_SPSAMPHIDDEN;
FP_TYPE gvps_sparse_circular GVPS_ARG_SPSAMP;

FP_TYPE gvps_full_sampled_hidden GVPS_ARG_FULLSAMPHIDDEN;
FP_TYPE gvps_full_sampled GVPS_ARG_FULLSAMP;

FP_TYPE gvps_full_circular_hidden GVPS_ARG_FULLSAMPHIDDEN;
FP_TYPE gvps_full_circular GVPS_ARG_FULLSAMP;

FP_TYPE gvps_sparse_sampled_hidden_static GVPS_ARG_SPSAMPHIDDEN;
FP_TYPE gvps_sparse_sampled_static GVPS_ARG_SPSAMP;

FP_TYPE gvps_sparse_circular_hidden_static GVPS_ARG_SPSAMPHIDDEN;
FP_TYPE gvps_sparse_circular_static GVPS_ARG_SPSAMP;

FP_TYPE gvps_full_sampled_hidden_static GVPS_ARG_FULLSAMPHIDDEN;
FP_TYPE gvps_full_sampled_static GVPS_ARG_FULLSAMP;

FP_TYPE gvps_full_circular_hidden_static GVPS_ARG_FULLSAMPHIDDEN;
FP_TYPE gvps_full_circular_static GVPS_ARG_FULLSAMP;

/* not implemented due to tractability/performance issues
FP_TYPE gvps_full_hidden(int* dst, void* context, gvps_obsrv* obsrv,
    gvps_ftran ftransame, gvps_ftran ftrandiff, gvps_ftrexist ftrexist);
*/
FP_TYPE gvps_full(int* dst, void* context, gvps_obsrv* obsrv,
    gvps_ftran ftran, gvps_ftrexist ftrexist);
FP_TYPE gvps_full_static(int* dst, void* context, gvps_obsrv* obsrv,
    gvps_ftran ftran);

FP_TYPE gvps_variable(int* dst, void* context, gvps_obsrv* obsrv,
    gvps_ftran ftran, gvps_ftrexist ftrexist);

#endif

