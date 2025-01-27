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
 */  

#include <stdlib.h>
#include "l7.h"
#include "l7p.h"

#define L7_LOCATION "L7_PUSH_UPDATE"

int L7_Push_Update(
      const int               *array,
      int                     *return_array,
      const int               l7_push_id
      )
{
   /*
    * Purpose
    * =======
    * L7_Update collects into array data_buffer data located off-process,
    * appending it to owned (on-process) data data_buffer.
    * 
    * Arguments
    * =========
    * array              (input) data input
    *
    * return_array       (output) data input
    * 
    * l7_push_id         (input) const int
    *                    Handle to database containing conmmunication
    *                    requirements.
    * 
    * Notes:
    * =====
    * 1) Serial compilation creates a no-op
    * 
    */
#if defined HAVE_MPI
   
   /*
    * Local variables
    */

   int
     ierr;                        /* Error code for return                */

   l7_push_id_database
     *l7_push_id_db;
   
   /*
    * Executable Statements
    */
   
   if (! l7.mpi_initialized){
      return(0);
   }
    
   if (l7.initialized !=1){
      ierr = 1;
      L7_ASSERT(l7.initialized == 1, "L7 not initialized", ierr);
   }
   
   /*
    * Check input.
    */

   if (array == NULL){
      ierr = -1;
      L7_ASSERT( array != NULL, "array != NULL", ierr);
   }
   
   if (l7_push_id <= 0){
      ierr = -1;
      L7_ASSERT( l7_push_id > 0, "l7_push_id <= 0", ierr);
   }
   
   if (l7.numpes == 1){
      ierr = L7_OK;
      return(ierr);
   }
   
   /*
    * find database associated with input l7_id
    */
   
   l7_push_id_db = l7.first_push_db;
   
   while (l7_push_id_db){
      if (l7_push_id_db->l7_push_id == l7_push_id){
            break;
      } else {
         /* Move to next one */
         l7_push_id_db = l7_push_id_db->next_push_db;
      }
   }

   if (l7_push_id_db == NULL){
      ierr = -1;
      L7_ASSERT(l7_push_id_db != NULL, "Failed to find database.", ierr);
   }
   
   if (l7.numpes == 1){ /* No-op */
      ierr = L7_OK;
      return(ierr);
   }
   
   int sizeof_type = l7p_sizeof(L7_INT);
   struct l7_update_datatype *dt = &l7_push_id_db->nbr_state.update_datatypes[sizeof_type];
   MPI_Neighbor_alltoallw(
			  array, l7_push_id_db->nbr_state.mpi_send_counts, 
			  (MPI_Aint *)l7_push_id_db->nbr_state.mpi_send_offsets, dt->out_types, 
			  return_array, l7_push_id_db->nbr_state.mpi_recv_counts, 
			  (MPI_Aint *)l7_push_id_db->nbr_state.mpi_recv_offsets, dt->in_types, 
			  l7_push_id_db->nbr_state.comm);
   
#endif /* HAVE_MPI */
   
   return(L7_OK);
    
} /* End L7_Update */

