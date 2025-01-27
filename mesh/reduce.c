/*
 *  Copyright (c) 2011-2019, Triad National Security, LLC.
 *  All rights Reserved.
 *
 *  CLAMR -- LA-CC-11-094
 *
 *  Copyright 2011-2019. Triad National Security, LLC. This software was produced 
 *  under U.S. Government contract 89233218CNA000001 for Los Alamos National 
 *  Laboratory (LANL), which is operated by Triad National Security, LLC 
 *  for the U.S. Department of Energy. The U.S. Government has rights to use, 
 *  reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR
 *  TRIAD NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR 
 *  ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified
 *  to produce derivative works, such modified software should be clearly marked,
 *  so as not to confuse it with the version available from LANL.
 *
 *  Additionally, redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the Triad National Security, LLC, Los Alamos 
 *       National Laboratory, LANL, the U.S. Government, nor the names of its 
 *       contributors may be used to endorse or promote products derived from 
 *       this software without specific prior written permission.
 *  
 *  THIS SOFTWARE IS PROVIDED BY THE TRIAD NATIONAL SECURITY, LLC AND 
 *  CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT 
 *  NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 *  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL TRIAD NATIONAL
 *  SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 *  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 *  OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 *  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *  
 *  CLAMR -- LA-CC-11-094
 *  This research code is being developed as part of the 
 *  2011 X Division Summer Workshop for the express purpose
 *  of a collaborative code for development of ideas in
 *  the implementation of AMR codes for Exascale platforms
 *  
 *  AMR implementation of the Wave code previously developed
 *  as a demonstration code for regular grids on Exascale platforms
 *  as part of the Supercomputing Challenge and Los Alamos 
 *  National Laboratory
 *  
 *  Authors: Bob Robey       XCP-2   brobey@lanl.gov
 *           Neal Davis              davis68@lanl.gov, davis68@illinois.edu
 *           David Nicholaeff        dnic@lanl.gov, mtrxknight@aol.com
 *           Dennis Trujillo         dptrujillo@lanl.gov, dptru10@gmail.com
 * 
 */
#include "reduce.h"
#ifdef HAVE_OPENCL
#include "ezcl/ezcl.h"
#endif

#ifdef HAVE_OPENCL
#include "reduce_kernel.inc"
#endif

#ifdef HAVE_OPENCL
/* Define the kernels declared in the header file */
cl_kernel   kernel_reduce_sum,
            kernel_reduce_sum_stage1of2,
            kernel_reduce_sum_stage2of2,
            kernel_reduce_sum_int_stage1of2,
            kernel_reduce_sum_int_stage2of2,
            kernel_reduce_product,
            kernel_reduce_max,
            kernel_reduce_max_stage1of2,
            kernel_reduce_max_stage2of2,
            kernel_reduce_min,
            kernel_reduce_min_stage1of2,
            kernel_reduce_min_stage2of2;
#endif

void init_kernels_reduce(void)
{
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_sum     = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_cl");
    kernel_reduce_sum_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_stage1of2_cl");
    kernel_reduce_sum_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_stage2of2_cl");
    kernel_reduce_sum_int_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_int_stage1of2_cl");
    kernel_reduce_sum_int_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_int_stage2of2_cl");
    kernel_reduce_product = ezcl_create_kernel_wsource(context, reduce_source, "reduce_product_cl");
    kernel_reduce_max     = ezcl_create_kernel_wsource(context, reduce_source, "reduce_max_cl");
    kernel_reduce_max_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_max_stage1of2_cl");
    kernel_reduce_max_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_max_stage2of2_cl");
    kernel_reduce_min     = ezcl_create_kernel_wsource(context, reduce_source, "reduce_min_cl");
    kernel_reduce_min_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_min_stage1of2_cl");
    kernel_reduce_min_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_min_stage2of2_cl");
#endif
}

void init_kernel_sum(void)
{
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_sum = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_cl");
#endif
}

void init_kernel_2stage_sum(void)
{
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_sum_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_stage1of2_cl");
    kernel_reduce_sum_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_stage2of2_cl");
#endif
}

void terminate_kernel_2stage_sum(void)
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_sum_stage1of2);
    ezcl_kernel_release(kernel_reduce_sum_stage2of2);
#endif
}

void init_kernel_2stage_sum_int(void)
{   
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_sum_int_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_int_stage1of2_cl");
    kernel_reduce_sum_int_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_sum_int_stage2of2_cl");
#endif
}

void terminate_kernel_2stage_sum_int(void)
{   
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_sum_int_stage1of2);
    ezcl_kernel_release(kernel_reduce_sum_int_stage2of2);
#endif
}

void init_kernel_product(void)
{
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_product = ezcl_create_kernel_wsource(context, reduce_source, "reduce_product_cl");
#endif
}

void init_kernel_max(void)
{   
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_max = ezcl_create_kernel_wsource(context, reduce_source, "reduce_max_cl");
#endif
}

void init_kernel_2stage_max(void)
{
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_max_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_max_stage1of2_cl");
    kernel_reduce_max_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_max_stage2of2_cl");
#endif
}

void init_kernel_min(void)
{   
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_min = ezcl_create_kernel_wsource(context, reduce_source, "reduce_min_cl");
#endif
}

void init_kernel_2stage_min(void)
{
#ifdef HAVE_OPENCL
    cl_context context = ezcl_get_context();
    kernel_reduce_min_stage1of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_min_stage1of2_cl");
    kernel_reduce_min_stage2of2 = ezcl_create_kernel_wsource(context, reduce_source, "reduce_min_stage2of2_cl");
#endif
}

void release_kernels_reduce()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_sum);
    ezcl_kernel_release(kernel_reduce_sum_stage1of2);
    ezcl_kernel_release(kernel_reduce_sum_stage2of2);
    ezcl_kernel_release(kernel_reduce_sum_int_stage1of2);
    ezcl_kernel_release(kernel_reduce_sum_int_stage2of2);
    ezcl_kernel_release(kernel_reduce_product);
    ezcl_kernel_release(kernel_reduce_max);
    ezcl_kernel_release(kernel_reduce_max_stage1of2);
    ezcl_kernel_release(kernel_reduce_max_stage2of2);
    ezcl_kernel_release(kernel_reduce_min);
    ezcl_kernel_release(kernel_reduce_min_stage1of2);
    ezcl_kernel_release(kernel_reduce_min_stage2of2);
#endif
}

void release_kernel_sum()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_sum);
#endif
}

void release_kernel_2stage_sum()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_sum_stage1of2);  
    ezcl_kernel_release(kernel_reduce_sum_stage2of2);
#endif
}

void release_kernel_2stage_sum_int()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_sum_int_stage1of2);  
    ezcl_kernel_release(kernel_reduce_sum_int_stage2of2);
#endif
}

void release_kernel_product()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_product);
#endif
}

void release_kernel_max()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_max);
#endif
}

void release_kernel_2stage_max()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_max_stage1of2);  
    ezcl_kernel_release(kernel_reduce_max_stage2of2);
#endif
}

void release_kernel_min()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_min);
#endif
}

void release_kernel_2stage_min()
{
#ifdef HAVE_OPENCL
    ezcl_kernel_release(kernel_reduce_min_stage1of2);  
    ezcl_kernel_release(kernel_reduce_min_stage2of2);
#endif
}

