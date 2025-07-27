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

#undef a_h
#undef p_h
#undef a
#undef p
#undef trandiff
#undef transame
#undef nstate_at
#undef state_at
#undef p_at
#undef tran

int T = obsrv -> T;
#ifdef FULL
int nstate = obsrv -> slice[0] -> N;
#define nstate_at(t) nstate
#define state_at(t, n) n
FP_TYPE*  a_ = calloc(T * nstate, sizeof(FP_TYPE));
int*      p_ = calloc(T * nstate, sizeof(int));
#define a(t, n) a_[(t) * nstate + (n)]
#define p(t, n) p_[(t) * nstate + (n)]
#else
int* nsparse = calloc(T, sizeof(int));
#define nstate_at(t) nsparse[t]
#define state_at(t, n) (obsrv -> slice[t] -> pair[n].state)
FP_TYPE** a_ = calloc(T, sizeof(FP_TYPE*));
int**     p_ = calloc(T, sizeof(int*));
#define a(t, n) a_[t][n]
#define p(t, n) p_[t][n]
#endif

#define p_at(t, n) (obsrv -> slice[t] -> pair[n].p)

#ifdef SAMPLED
#ifdef STATIC
FP_TYPE*  transame_ = calloc(nstate, sizeof(FP_TYPE));
#else
FP_TYPE*  transame_ = calloc(T * nprestore, sizeof(FP_TYPE));
#endif

#ifdef HIDDEN
FP_TYPE*  a_h_ = calloc(T * nstate, sizeof(FP_TYPE));
int*      p_h_ = calloc(T * nstate, sizeof(int));
#define a_h(t, n) a_h_[(t) * nstate + (n)]
#define p_h(t, n) p_h_[(t) * nstate + (n)]
#ifdef STATIC
FP_TYPE*  trandiff_ = calloc(nstate, sizeof(FP_TYPE));
#define trandiff(t, d) trandiff_[d]
#define transame(t, d) transame_[d]
#else
FP_TYPE*  trandiff_ = calloc(T * nprestore, sizeof(FP_TYPE));
#define trandiff(t, d) (((d) < nprestore) ? trandiff_[(t) * nprestore + (d)] : log((*ftrandiff)(context, d, t)))
#define transame(t, d) (((d) < nprestore) ? transame_[(t) * nprestore + (d)] : log((*ftransame)(context, d, t)))
#endif
#else
#ifdef STATIC
#define transame(t, d) transame_[d]
#else
#define transame(t, d) (((d) < nprestore) ? transame_[(t) * nprestore + (d)] : log((*ftran)(context, d, t)))
#endif
#endif

#else
#ifdef STATIC
FP_TYPE* tran_ = calloc(nstate * nstate, sizeof(FP_TYPE));
#define transame(t2, s1, s2) tran_[(s1) * nstate + (s2)]
#else
#define transame(t2, s1, s2) (*ftran)(context, s1, s2, t2)
#endif
#endif

int t, i, j;
t = i = j = 0;

#ifdef STATIC
#ifdef SAMPLED
int d;
for(d = 0; d < nstate; d ++)
{
#ifdef HIDDEN
    transame(0, d) = log((*ftransame)(context, d, 0));
    trandiff(0, d) = log((*ftrandiff)(context, d, 0));
#else
    transame(0, d) = log((*ftran)(context, d, 0));
#endif
}
#else
for(i = 0; i < nstate; i ++)
    for(j = 0; j < nstate; j ++)
        transame(0, i, j) = (*ftran)(context, i, j, 0);
#endif
#endif

#ifdef SAMPLED
for(t = 0; t < T; t ++)
{
#ifndef FULL
    nsparse[t] = obsrv -> slice[t] -> N;
    a_[t] = calloc(nsparse[t], sizeof(FP_TYPE));
    p_[t] = calloc(nsparse[t], sizeof(int));
#endif
#ifndef STATIC
    int d;
    for(d = 0; d < nprestore; d ++)
    {
#ifdef HIDDEN
        transame_[(t) * nprestore + d] = log((*ftransame)(context, d, t));
        trandiff_[(t) * nprestore + d] = log((*ftrandiff)(context, d, t));
#else
        transame_[(t) * nprestore + d] = log((*ftran)(context, d, t));
#endif
    }
#endif
}
#endif

#ifdef VARIABLE
for(t = 0; t < T; t ++) {
    nsparse[t] = obsrv -> slice[t] -> N;
    a_[t] = calloc(nsparse[t], sizeof(FP_TYPE));
    p_[t] = calloc(nsparse[t], sizeof(int));
}
#endif

// initialization
#ifdef HIDDEN
FP_TYPE pexposed = 0;
#endif
for(i = 0; i < nstate_at(0); i ++)
{
    FP_TYPE P = p_at(0, i);
    a(0, i) = log(P);
#ifdef HIDDEN
    pexposed += P;
#endif
}
#ifdef HIDDEN
FP_TYPE phidden = log(1 - pexposed < 0.00001 ? 0.00001 : 1 - pexposed);
for(i = 0; i < nstate; i ++)
    a_h(0, i) = phidden;
#endif

// DP searching
for(t = 1; t < T; t ++)
{
    int n, nprv;
    FP_TYPE maxp;
    int     maxj;
#ifdef SAMPLED
    int ntrans = (*fntran)(context, t);
#endif
#ifdef HIDDEN
    int ntransh = min(ntrans, nhiddenprune);
    pexposed = 0;
#endif

#undef transition_update
#ifdef CIRCULAR
    #define transition_update(f, a, newj) \
        int diff = abs(i - j); \
        FP_TYPE psum = a + f(t, diff > nstate / 2 ? nstate - diff : diff); \
        if(psum > maxp) \
        { \
            maxp = psum; \
            maxj = newj; \
        }
#else
#ifndef SAMPLED
    #define transition_update(f, a, newj) \
        FP_TYPE psum = a + f(t, j, i); \
        if(psum > maxp) \
        { \
            maxp = psum; \
            maxj = newj; \
        }
#else
    #define transition_update(f, a, newj) \
        FP_TYPE psum = a + f(t, abs(i - j)); \
        if(psum > maxp) \
        { \
            maxp = psum; \
            maxj = newj; \
        }
#endif
#endif
    
    // exposed/hidden -> exposed
    for(n = 0; n < nstate_at(t); n ++)
    {
        maxp = -30e10;
        maxj = 0;
        i = state_at(t, n);
#ifdef HIDDEN
        // hidden -> exposed
#ifdef CIRCULAR
        for(int cj = i - ntrans; cj < i + ntrans + 1; cj ++)
        {
            j = (cj + nstate) % nstate;
            transition_update(trandiff, a_h(t - 1, j), -j - 1);
        }
#else
        for(j = max(0, i - ntrans); j < min(nstate, i + ntrans + 1); j ++)
        {
            transition_update(trandiff, a_h(t - 1, j), -j - 1);
        }
#endif
#endif
        // exposed -> exposed
        for(nprv = 0; nprv < nstate_at(t - 1); nprv ++)
        {
            j = state_at(t - 1, nprv);
#ifdef SAMPLED
#ifndef CIRCULAR
            if(abs(i - j) > ntrans) continue;
#endif
#endif
#ifdef FTREXIST
            if(! (*ftrexist)(context, j, i, t)) continue;
#endif
            transition_update(transame, a(t - 1, nprv), nprv);
        }
#ifdef HIDDEN
        pexposed += p_at(t, n);
#endif
        a(t, n) = maxp + log(p_at(t, n));
        p(t, n) = maxj;
    }
    
#ifdef HIDDEN
    phidden = log(1 - pexposed < 0.00001 ? 0.00001 : 1 - pexposed);
    // exposed/hidden -> hidden
    for(i = 0; i < nstate; i ++)
    {
        maxp = -30e10;
        maxj = 0;
        // hidden -> hidden
#ifdef CIRCULAR
        for(int cj = i - ntransh; cj < i + ntransh + 1; cj ++)
        {
            j = (cj + nstate) % nstate;
            transition_update(transame, a_h(t - 1, j), -j - 1);
        }
#else
        for(j = max(0, i - ntransh); j < min(nstate, i + ntransh + 1); j ++)
        {
            transition_update(transame, a_h(t - 1, j), -j - 1);
        }
#endif
        
        // exposed -> hidden
        for(nprv = 0; nprv < nstate_at(t - 1); nprv ++)
        {
            j = state_at(t - 1, nprv);
#ifndef CIRCULAR
            if(abs(i - j) > ntrans) continue;
#endif
            transition_update(trandiff, a(t - 1, nprv), nprv);
        }
        
        a_h(t, i) = maxp + phidden;
        p_h(t, i) = maxj;
    }
#endif
#ifdef SAMPLED
    (void)ntrans; // suppress unused variable warning
#endif
}

// termination
FP_TYPE optim = -30e10;
int maxi = 0;
t = T - 1;
#ifdef HIDDEN
// hidden
for(i = 0; i < nstate; i ++)
    if(a_h(t, i) > optim)
    {
        optim = a_h(t, i);
        maxi = -i - 1;
    }
#endif
// exposed
for(i = 0; i < nstate_at(t); i ++)
    if(a(t, i) > optim)
    {
        optim = a(t, i);
        maxi = i;
    }

// back tracking
dst[t] = maxi;
while(t > 0)
{
    t --;
#ifdef HIDDEN
    dst[t] = dst[t + 1] < 0 ? p_h(t + 1, -dst[t + 1] - 1) : p(t + 1, dst[t + 1]);
#else
    dst[t] = p(t + 1, dst[t + 1]);
#endif
}
#ifdef HIDDEN
for(t = 0; t < T; t ++) // offset hidden states by nstate
    dst[t] = dst[t] < 0 ? nstate - dst[t] - 1 : state_at(t, dst[t]);
#else
for(t = 0; t < T; t ++) // pair index to state index
    dst[t] = state_at(t, dst[t]);
#endif

// free
#ifndef FULL
for(t = 0; t < T; t ++)
{
    free(a_[t]);
    free(p_[t]);
}
free(nsparse);
#endif
#if defined(STATIC) && ! defined(SAMPLED)
free(tran_);
#endif
free(a_); free(p_);
#ifdef SAMPLED
free(transame_);
#endif
#ifdef HIDDEN
free(a_h_); free(p_h_);
free(trandiff_);
#endif

return optim;

