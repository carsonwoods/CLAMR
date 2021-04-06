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

typedef struct val_rank {
  double value;
  int rank;
};

typedef struct int_rank {
  int value;
  int rank;
};

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
   
   /*
    * Set some parameters base on input datatype.
    */
   double begin_setup = MPI_Wtime();
   for (int ip = 0; ip < l7_push_id_db->num_comm_partners; ip++){
      int count = l7_push_id_db->send_buffer_count[ip]; // for vectorization
      for (int ic = 0; ic < count; ic++){
         l7_push_id_db->send_buffer[ip][ic] = array[l7_push_id_db->send_database[ip][ic]];
      }    
   }    


// Send/Receives will be done in L7_Push_Update. Input will be send_buffer. Output will be in
// preallocated receive_buffer
   MPI_Request request[2*l7_push_id_db->num_comm_partners];
   MPI_Status  status[2*l7_push_id_db->num_comm_partners];
  
   int iloc = 0;
   // Barrier then begin
   double begin_Irecv = MPI_Wtime();

   for (int ip = 0; ip < l7_push_id_db->num_comm_partners; ip++){
      MPI_Irecv(&return_array[iloc], l7_push_id_db->recv_buffer_count[ip], MPI_INT,
                l7_push_id_db->comm_partner[ip], l7_push_id_db->comm_partner[ip], MPI_COMM_WORLD, &request[ip]);
      iloc += l7_push_id_db->recv_buffer_count[ip];
   }
   
   double begin_Isend = MPI_Wtime();

   for (int ip = 0; ip < l7_push_id_db->num_comm_partners; ip++){
      MPI_Isend(l7_push_id_db->send_buffer[ip], l7_push_id_db->send_buffer_count[ip], MPI_INT,
                l7_push_id_db->comm_partner[ip], l7.penum, MPI_COMM_WORLD, &request[l7_push_id_db->num_comm_partners+ip]);
   }

   double end_Isend = MPI_Wtime();

   MPI_Waitall(2*l7_push_id_db->num_comm_partners, request, status);

   double end_waitall = MPI_Wtime();

   double setup_time = begin_Irecv - begin_setup;
   double Irecv_time = begin_Isend - begin_Irecv;
   double Isend_time = end_Isend - begin_Isend;
   double waitall_time = end_waitall - end_Isend;

   int rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);


   struct val_rank Irecv_min_loc, Irecv_max_loc, Irecv_min_loc_out, Irecv_max_loc_out;
   struct val_rank Isend_min_loc, Isend_max_loc, Isend_min_loc_out, Isend_max_loc_out;
   struct val_rank waitall_min_loc, waitall_max_loc, waitall_min_loc_out, waitall_max_loc_out;
   struct val_rank setup_min_loc, setup_max_loc, setup_min_loc_out, setup_max_loc_out;

   Irecv_min_loc.value = Irecv_time;
   Irecv_min_loc.rank = rank;
   Irecv_max_loc.value = Irecv_time;
   Irecv_max_loc.rank = rank;
   Isend_min_loc.value = Isend_time;
   Isend_min_loc.rank = rank;
   Isend_max_loc.value = Isend_time;
   Isend_max_loc.rank = rank;
   waitall_min_loc.value = waitall_time;
   waitall_min_loc.rank = rank;
   waitall_max_loc.value = waitall_time;
   waitall_max_loc.rank = rank;
   setup_min_loc.value = setup_time;
   setup_min_loc.rank = rank;
   setup_max_loc.value = setup_time;
   setup_max_loc.rank = rank;




   double Irecv_sum,Irecv_average;
   double Isend_sum,Isend_average;
   double waitall_sum, waitall_average;
   double setup_sum, setup_average;

   // Structs to get number of partners  
   struct int_rank num_partners_min_loc, num_partners_max_loc, num_partners_min_loc_out, num_partners_max_loc_out;
   int num_partners_sum;
   double num_partners_average;

   num_partners_min_loc.value = l7_push_id_db->num_comm_partners;
   num_partners_min_loc.rank = rank;
   num_partners_max_loc.value = l7_push_id_db->num_comm_partners;
   num_partners_max_loc.rank = rank;



   // TODO: Define and implement struct and MPI_MINLOC, MAXLOC 
   // Min for Irecv_time;
   MPI_Reduce(&Irecv_min_loc, &Irecv_min_loc_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
   // Max for Irecv_time;
   MPI_Reduce(&Irecv_max_loc, &Irecv_max_loc_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   // SUM for Irecv_time;
   MPI_Reduce(&Irecv_time, &Irecv_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   // Min for Isend_time;
   MPI_Reduce(&Isend_min_loc, &Isend_min_loc_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
   // Max for Isend_time;
   MPI_Reduce(&Isend_max_loc, &Isend_max_loc_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   // SUM for Isend_time;
   MPI_Reduce(&Isend_time, &Isend_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   // Min for waitall_time;
   MPI_Reduce(&waitall_min_loc, &waitall_min_loc_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
   // Max for waitall_time;
   MPI_Reduce(&waitall_max_loc, &waitall_max_loc_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   // SUM for waitall_time;
   MPI_Reduce(&waitall_time, &waitall_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   // Min for setup_time;
   MPI_Reduce(&setup_min_loc, &setup_min_loc_out, 1, MPI_DOUBLE_INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
   // Max for setup_time;
   MPI_Reduce(&setup_max_loc, &setup_max_loc_out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   // SUM for setup_time;
   MPI_Reduce(&setup_time, &setup_sum, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

   // Min for num partners
   MPI_Reduce(&num_partners_min_loc, &num_partners_min_loc_out, 1, MPI_2INT, MPI_MINLOC, 0, MPI_COMM_WORLD);
   // Max for num partners
   MPI_Reduce(&num_partners_max_loc, &num_partners_max_loc_out, 1, MPI_2INT, MPI_MAXLOC, 0, MPI_COMM_WORLD);
   // SUM of num partners
   MPI_Reduce(&l7_push_id_db->num_comm_partners, &num_partners_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  
   // Grabbing comm size to calculate average
   int comm_size;
   MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

   Irecv_average = Irecv_sum / comm_size;
   Isend_average = Isend_sum / comm_size;
   waitall_average = waitall_sum / comm_size;
   setup_average = setup_sum / comm_size;
   num_partners_average = num_partners_sum / (comm_size * 1.0f);

   // TODO: Print to std error
   if (rank == 0) {
      fprintf(stderr, "\n=====================\n");
      fprintf(stderr, "\nsetup min: %f (rank: %d)\n", setup_min_loc_out.value, setup_min_loc_out.rank);
      fprintf(stderr, "\nsetup max: %f (rank: %d)\n", setup_max_loc_out.value, setup_max_loc_out.rank);
      fprintf(stderr, "\nsetup avg: %f\n", setup_average);
      fprintf(stderr, "\nIrecv min: %f (rank: %d)\n", Irecv_min_loc_out.value, Irecv_min_loc_out.rank);
      fprintf(stderr, "\nIrecv max: %f (rank: %d)\n", Irecv_max_loc_out.value, Irecv_max_loc_out.rank);
      fprintf(stderr, "\nIrecv avg: %f\n", Irecv_average);
      fprintf(stderr, "\nIsend min: %f (rank: %d)\n", Isend_min_loc_out.value, Isend_min_loc_out.rank);
      fprintf(stderr, "\nIsend max: %f (rank: %d)\n", Isend_max_loc_out.value, Isend_max_loc_out.rank);
      fprintf(stderr, "\nIsend avg: %f\n", Isend_average);
      fprintf(stderr, "\nwaitall min: %f (rank: %d)\n", waitall_min_loc_out.value, waitall_min_loc_out.rank);
      fprintf(stderr, "\nwaitall max: %f (rank: %d)\n", waitall_max_loc_out.value, waitall_max_loc_out.rank);
      fprintf(stderr, "\nwaitall avg: %f\n", waitall_average);
      fprintf(stderr, "\nNum partners min: %d (rank: %d)\n", num_partners_min_loc_out.value, num_partners_min_loc_out.rank);
      fprintf(stderr, "\nNum partners max: %d (rank: %d)\n", num_partners_max_loc_out.value, num_partners_max_loc_out.rank);
      fprintf(stderr, "\nNum partners avg: %f\n", num_partners_average);
  }
  

/*
   if (ncycle >= 1) { 
      for (int ib = 0; ib<receive_count_total; ib++){
         fprintf(fp,"DEBUG receive %d is %d\n",ib,border_data_receive[ib]);
      }    
   }    
*/
   
#endif /* HAVE_MPI */
   
   return(L7_OK);
    
} /* End L7_Update */

