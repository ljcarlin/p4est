/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/** \file gmt2.c
 * 
 * Search based refinement. There are 3 models: synthetic, sphere and latlong.
 * 
 * Usage of the sphere model follows the following pipeline.
 *  -# Prepare a csv file of geodesics, following the convention described
 *     in \ref sphere_preprocessing.c . Note that world-map datasets frequently
 *     come in .shp files. The raw line geodesic data can be extracted from
 *     these easily using the python geopandas library. However, there is 
 *     currently only limited library support for .shp files in C.
 *  -# Run the preprocessing script as described in \ref sphere_preprocessing.c
 *  -# Run the sphere model: p4est_gmt --sphere --resolution <resolution>
*/

#include <p4est_extended.h>
#include <p4est_search.h>
#include <p4est_vtk.h>
#include <sc_options.h>
#include "gmt_models.h"

static const double irootlen = 1. / (double) P4EST_ROOT_LEN;

typedef struct global
{
  int                 minlevel;
  int                 maxlevel;
  int                 balance;
  int                 resolution;
  int                 synthetic;
  int                 latlongno;
  int                 sphere; /* globe sphere model */
  int                 distributed; /* are we running the distributed version */
  sc_MPI_Comm         mpicomm;
  p4est_t            *p4est;
  p4est_gmt_model_t  *model;
}
global_t;

static int
setup_model (global_t * g)
{
  /* this function populates the model on successful initialization */
  P4EST_ASSERT (g->model == NULL);

  /* initialize model depending on command line options */
  if (g->synthetic >= 0) {
    switch (g->synthetic) {
    case 0:
      g->model = p4est_gmt_model_synth_new (g->synthetic, g->resolution);
      break;
    default:
      P4EST_GLOBAL_LERROR ("Synthetic model number exceeded\n");
    }
  }
  else if (g->latlongno >= 0) {
    p4est_gmt_model_latlong_params_t ap;

    switch (g->latlongno) {
    case 0:
      ap.latitude[0] = -50.;
      ap.latitude[1] = 0.;
      ap.longitude[0] = 0.;
      ap.longitude[1] = 60.;
      ap.resolution = g->resolution;
      ap.load_filename = "africa.gmt.data";
      ap.output_prefix = "africa";

      /* load data (possibly GMT, or file, or synthetic) */
      g->model = p4est_gmt_model_latlong_new (&ap);
      break;
    default:
      P4EST_GLOBAL_LERROR ("Latitute-longitude model number exceeded\n");
    }
  }
  else if (g->sphere == 1) {
    g->model = p4est_gmt_model_sphere_new(g->resolution, g->distributed, g->mpicomm);
  }

  /* on successful initalization the global model is set */
  return g->model == NULL ? -1 : 0;
}

static void
quad_init (p4est_t * p4est,
           p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  quadrant->p.user_int = 0;
}

static int
quad_refine (p4est_t * p4est,
             p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  return quadrant->p.user_int;
}

static int
quad_point (p4est_t * p4est,
            p4est_topidx_t which_tree, p4est_quadrant_t * quadrant,
            p4est_locidx_t local_num, void *point)
{
  int                 result;
  size_t              m;
  double              coord[4];
  p4est_qcoord_t      qh;
  p4est_gmt_model_t  *model;
  global_t           *g = (global_t *) p4est->user_pointer;

  /* sanity checks */
  P4EST_ASSERT (g != NULL);
  P4EST_ASSERT (g->p4est == p4est);
  model = g->model;
  P4EST_ASSERT (model != NULL);
  P4EST_ASSERT (model->intersect != NULL);

  /* retrieve object index and model */
  P4EST_ASSERT (point != NULL);
  m = *(size_t *) point;
  P4EST_ASSERT (m < model->M);

  /* provide rectangle coordinates */
  qh = P4EST_QUADRANT_LEN (quadrant->level);
  coord[0] = irootlen * quadrant->x;
  coord[1] = irootlen * quadrant->y;
  coord[2] = irootlen * (quadrant->x + qh);
  coord[3] = irootlen * (quadrant->y + qh);

  /* execute intersection test */
  if ((result = g->model->intersect (which_tree, coord, m, model)) &&
      local_num >= 0 && quadrant->level < g->maxlevel) {
    /* set refinement indicator for a leaf quadrant */
    quadrant->p.user_int = 1;
  }
  return result;
}

void
run_program (global_t * g)
{
  int                 refiter;
  size_t              zz;
  char                filename[BUFSIZ];
  sc_array_t         *points;
  p4est_gloidx_t      gnq_before;
  const size_t        quad_data_size = 0;

  /* create mesh */
  P4EST_GLOBAL_PRODUCTION ("Create initial mesh\n");
  g->p4est = p4est_new_ext (g->mpicomm, g->model->conn, 0, g->minlevel, 1,
                            quad_data_size, quad_init, g);

  /* run mesh refinement based on data */
  P4EST_GLOBAL_PRODUCTIONF ("Setting up %lld search objects\n",
                            (long long) g->model->M);
  points = sc_array_new_count (sizeof (size_t), g->model->M);

  for (zz = 0; zz < g->model->M; ++zz) {
    *(size_t *) sc_array_index (points, zz) = zz;
  }
  for (refiter = 0;; ++refiter) {
    P4EST_GLOBAL_PRODUCTIONF ("Into refinement iteration %d\n", refiter);
    snprintf (filename, BUFSIZ, "p4est_gmt_%s_%02d",
              g->model->output_prefix, refiter);
    P4EST_ASSERT(g->model != NULL);
    P4EST_ASSERT(g->model->model_geom != NULL);
    p4est_vtk_write_file (g->p4est, g->model->model_geom, filename);
    gnq_before = g->p4est->global_num_quadrants;

    P4EST_GLOBAL_PRODUCTION ("Run object search\n");
    p4est_search_reorder (g->p4est, 1, NULL, NULL, NULL, quad_point, points);

    P4EST_GLOBAL_PRODUCTION ("Run mesh refinement\n");
    p4est_refine (g->p4est, 0, quad_refine, quad_init);

    if (g->balance) {
      P4EST_GLOBAL_PRODUCTION ("Run 2:1 mesh balance\n");
      p4est_balance (g->p4est, P4EST_CONNECT_FULL, quad_init);
    }

    if (gnq_before < g->p4est->global_num_quadrants) {
      P4EST_GLOBAL_PRODUCTION ("Run mesh repartition\n");
      p4est_partition (g->p4est, 0, NULL);
    }
    else {
      P4EST_GLOBAL_PRODUCTION ("Done refinement iterations\n");
      break;
    }
  }
  sc_array_destroy (points);

  /* cleanup */
  p4est_destroy (g->p4est);
}

/* TODO: This is a WIP version of run_program where each process only knows
 *  some of the points being searched. 
 *  
 * While refinement makes progress (globally):
 *    For each process q:
 *        Initialise empty array p_sending_to_q
 *    For each point x process p knows:
 *        for each process q in owners(x):
 *            Add point x to array p_sending_to_q
 * 
 *    Send each array and aggregate incoming arrays.
 *    The incoming arrays are now the points process p knows.
 * 
 *    Remove duplicates (is there some clever way to do this?)
 *    (We could give each point x a unique owner who is responsible for
 *    passing messages, and every other process knowing x only uses it
 *    to refine.)
 */
void
run_program_distributed (global_t * g) 
{
  int                 mpiret;
  int                 refiter;
  size_t              zz;
  char                filename[BUFSIZ];
  sc_array_t         *points;
  p4est_gloidx_t      gnq_before;
  const size_t        quad_data_size = 0;
  int                 num_procs, rank;
  void               *send_buffer;

  /**TODO: this is hardcoded at the moment. Do it in a more universal way **/
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)g->model->model_data;

  /** communication data for each iteration. p is the local process **/
  /* q -> {points q should own in the next iteration, and be responsible for} */
  sc_array_t        **owned_resp; 
  /* q -> {points q should own in the next iteration, but not be responsible for } */
  sc_array_t        **owned;
  /* number of points that p should own in the next iteration (responsible)*/
  int                 num_owned_resp; 
  /* number of points that p should own in the next iteration (not responsible)*/
  int                 num_owned;
  /* q -> number of bytes we are receiving from q  */
  int                *recv_counts;
  /* q -> byte offset to receive at */
  int                *displs;
  
  /** per point data **/
  sc_array_t         *owners; /* who to send the point to */
  int                 responsible; /* who is in charge of propagating the point */


  /* Get rank and total process count */
  mpiret = sc_MPI_Comm_size (g->mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (g->mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* create mesh */
  P4EST_GLOBAL_PRODUCTION ("Create initial mesh\n");
  g->p4est = p4est_new_ext (g->mpicomm, g->model->conn, 0, g->minlevel, 1,
                            quad_data_size, quad_init, g);

  /* Initialise index of outgoing message arrays */
  owned = P4EST_ALLOC(sc_array_t*, num_procs);
  owned_resp = P4EST_ALLOC(sc_array_t*, num_procs);

  /* Initialise recv_count and displacement arrays */
  recv_counts = P4EST_ALLOC(int, num_procs);
  displs = P4EST_ALLOC(int, num_procs);

  /* Initially all points are read once and the process that reads them
  is responsible for their propagation*/
  num_owned = 0;
  num_owned_resp = g->model->M;
  printf("Proc %d owns %d points\n", rank, num_owned_resp);
 
  for (refiter = 0;; ++refiter) {
    P4EST_GLOBAL_PRODUCTIONF ("Into refinement iteration %d\n", refiter);
    snprintf (filename, BUFSIZ, "p4est_gmt_%s_%02d",
              g->model->output_prefix, refiter);
    P4EST_ASSERT(g->model != NULL);
    P4EST_ASSERT(g->model->model_geom != NULL);
    p4est_vtk_write_file (g->p4est, g->model->model_geom, filename);
    gnq_before = g->p4est->global_num_quadrants;

    printf("Proc %d: Initialise array\n", rank);
    /* Initialise arrays for communicating points */
    for (int q = 0; q < num_procs; q++) {
      owned[q] = sc_array_new(g->model->point_size);
      owned_resp[q] = sc_array_new(g->model->point_size);
    }

    printf("Proc %d: Prepare outgoing\n", rank);
    /* Prepare outgoing arrays of points we are responsible for */
    for (int i = 0; i < num_owned_resp; i++) {
      printf("Proc %d point %d\n", rank, i);
      /* Determine which processes should receive this point */
      printf("Proc %d owners\n", rank);
      owners = g->model->owners_fn(&(sdata->points[i]), g->p4est);
      /* Determine the receiver responsible for propagating this point in next
         iteration */
      printf("Proc %d responsible\n", rank);
      responsible = g->model->resp_fn(&(sdata->points[i]), owners);

      printf("Proc %d add to receivers\n", rank);
      /* Add this point to the arrays of the processes who receive it */
      for (int j = 0; j < owners->elem_count; j++) {
        int q = *(int*)sc_array_index_int(owners, j);
        if (responsible == q)
        {
          /* Process q should own point i and be responsible for its propagation */
          sc_array_push(owned_resp[q]);
          *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_index_int(owned_resp[q], owned_resp[q]->elem_count-1) = sdata->points[i];
        }
        else 
        {
          /* Process q should own point i but not be responsible for its propagation */
          sc_array_push(owned[q]);
          *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_index_int(owned[q], owned[q]->elem_count-1) = sdata->points[i];
        }
      }
      sc_array_destroy(owners);
    } 

    printf("Proc %d: Send points counts\n", rank);
    /* Tell each process q how many points it will own */
    for (int q = 0; q < num_procs; q++) { 

      mpiret = sc_MPI_Reduce(&(owned[q]->elem_count), &num_owned, 1, sc_MPI_INT, 
                        sc_MPI_SUM, q, g->mpicomm);
      SC_CHECK_MPI (mpiret);
      mpiret = sc_MPI_Reduce(&(owned_resp[q]->elem_count), &num_owned_resp, 1, sc_MPI_INT, 
                        sc_MPI_SUM, q, g->mpicomm);
      SC_CHECK_MPI (mpiret);
    }

    printf("Proc %d: Free old points\n", rank);
    /* Free the points we received last iteration */
    P4EST_FREE(sdata->points);
    // /* Create array for incoming points */
    sdata->points = P4EST_ALLOC(p4est_gmt_sphere_geodesic_seg_t, num_owned_resp + num_owned);
    /* Update the amount of points known */
    g->model->M = num_owned_resp + num_owned;

    printf("Proc %d: Send points\n", rank);
    /* Transmit points */
    for (int q = 0; q < num_procs; q++) 
      /* Tell process q how many (resp) bytes it will receive from each process */
      mpiret = sc_MPI_Gather(&((owned_resp[q]->elem_count) * sizeof(p4est_gmt_sphere_geodesic_seg_t)),
                                sizeof(int),
                                sc_MPI_BYTE,
                                recv_counts,
                                sizeof(int),
                                sc_MPI_BYTE,
                                q,
                                g->mpicomm
                            );
      SC_CHECK_MPI (mpiret);

      /* Compute receiver offsets (only relevant for process q) */
      if (rank == q) {
        displs[0] = 0;
        for (int i=1; i < num_procs; i++) {
          displs[i] = displs[i-1] + recv_counts[i-1];
        }
      }

      /* Set send buffer */
      if (owned_resp[q]->elem_count != 0) {
        send_buffer = sc_array_index_int(owned_resp[q], 0);
      }
      else {
        send_buffer = NULL;
      }

      /* Gather points (resp) to process q */
      mpiret = sc_MPI_Gatherv(send_buffer, 
                                (owned_resp[q]->elem_count) * (g->model->point_size),
                                sc_MPI_BYTE,
                                sdata->points,
                                recv_counts,
                                displs,
                                sc_MPI_BYTE,
                                q,
                                g->mpicomm);
      SC_CHECK_MPI (mpiret);

      /* Tell process q how many (resp) bytes it will receive from each process */
      mpiret = sc_MPI_Gather(&((owned[q]->elem_count) * sizeof(p4est_gmt_sphere_geodesic_seg_t)),
                                sizeof(int),
                                sc_MPI_BYTE,
                                recv_counts,
                                sizeof(int),
                                sc_MPI_BYTE,
                                q,
                                g->mpicomm
                            );
      SC_CHECK_MPI (mpiret);

      /* Compute receiver offsets (only relevant for process q) */
      if (rank == q) {
        displs[0] = 0;
        for (int i=1; i < num_procs; i++) {
          displs[i] = displs[i-1] + recv_counts[i-1];
        }
      }
      

      /* Set send buffer */
      if (owned[q]->elem_count != 0) {
        send_buffer = sc_array_index_int(owned[q], 0);
      }
      else {
        send_buffer = NULL;
      }

      /* Gather points to process q TODO */
      mpiret = sc_MPI_Gatherv(send_buffer, 
                                (owned[q]->elem_count) * (g->model->point_size),
                                sc_MPI_BYTE,
                                &(sdata->points[num_owned_resp]), /* rec buffer is offset */
                                (owned[q]->elem_count) * (g->model->point_size),
                                sc_MPI_BYTE,
                                q,
                                g->mpicomm);
      SC_CHECK_MPI (mpiret);
    }

    printf("Proc %d: cleanup sending arrays \n", rank);
    /* Clean up arrays used for sending points */
    for (int q = 0; q < num_procs; q++) {
      sc_array_destroy(owned[q]);
      sc_array_destroy(owned_resp[q]);
    }

    /* Set up search objects */
    P4EST_GLOBAL_PRODUCTIONF ("Setting up %lld search objects\n",
                            (long long) g->model->M);
    points = sc_array_new_count (sizeof (size_t), g->model->M); //TODO

    for (zz = 0; zz < g->model->M; ++zz) {
      *(size_t *) sc_array_index (points, zz) = zz;
    }

    /* Set refinement flags for leaves intersecting points this process owns */
    P4EST_GLOBAL_PRODUCTION ("Run object search\n");
    p4est_search_reorder (g->p4est, 1, NULL, NULL, NULL, quad_point, points);

    /* Refine based on refinement flags */
    P4EST_GLOBAL_PRODUCTION ("Run mesh refinement\n");
    p4est_refine (g->p4est, 0, quad_refine, quad_init);

    /* Cleanup points */
    sc_array_destroy (points);

    if (g->balance) {
      P4EST_GLOBAL_PRODUCTION ("Run 2:1 mesh balance\n");
      p4est_balance (g->p4est, P4EST_CONNECT_FULL, quad_init);
    }

    if (gnq_before < g->p4est->global_num_quadrants) {
      P4EST_GLOBAL_PRODUCTION ("Run mesh repartition\n");
      p4est_partition (g->p4est, 0, NULL);
    }
    else {
      P4EST_GLOBAL_PRODUCTION ("Done refinement iterations\n");
      break;
    }
  }
  

  /* cleanup */
  p4est_destroy (g->p4est);
}

static int
usagerrf (sc_options_t * opt, const char *fmt, ...)
{
  va_list             ap;
  char                msg[BUFSIZ];

  va_start (ap, fmt);
  vsnprintf (msg, BUFSIZ, fmt, ap);
  va_end (ap);

  P4EST_GLOBAL_LERROR ("ERROR/\n");
  P4EST_GLOBAL_LERRORF ("ERROR: %s\n", msg);
  P4EST_GLOBAL_LERROR ("ERROR\\\n");
  return 1;
}

static int
usagerr (sc_options_t * opt, const char *msg)
{
  return usagerrf (opt, "%s", msg);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 ue, fa;
  sc_options_t       *opt;
  global_t            sg, *g = &sg;

  /* initialize MPI */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* initialize global context */
  memset (g, 0, sizeof (*g));
  g->mpicomm = sc_MPI_COMM_WORLD;

  /* set global logging options for p4est */
  sc_init (g->mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initialize global application state */
  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'l', "minlevel", &g->minlevel, 0,
                      "Minimum refinement level");
  sc_options_add_int (opt, 'L', "maxlevel", &g->maxlevel, P4EST_QMAXLEVEL,
                      "Maximum refinement level");
  sc_options_add_bool (opt, 'b', "balance", &g->balance, 0,
                       "Execute 2:1 balance algorithm");
  sc_options_add_int (opt, 'r', "resolution", &g->resolution, 0,
                      "Level of resolution (model specific)");
  sc_options_add_int (opt, 'S', "synthetic", &g->synthetic, -1,
                      "Choose specific synthetic model");
  sc_options_add_int (opt, 'M', "latlongno", &g->latlongno, -1,
                      "Choose specific latitude-longitude model");
  sc_options_add_bool (opt, 'W', "sphere", &g->sphere, 0,
                       "Use sphere model");
  sc_options_add_bool (opt, 'd', "distributed", &g->distributed, 0,
                       "Distributed read mode");

  /* proceed in run-once loop for cleaner error checking */
  ue = 0;
  do {
    /* parse command line and assign configuration variables */
    fa = sc_options_parse (p4est_package_id, SC_LP_DEFAULT, opt, argc, argv);
    if (fa < 0 || fa != argc) {
      ue = usagerr (opt, "invalid option format or non-option argument");
      break;
    }
    P4EST_GLOBAL_PRODUCTIONF ("Manifold dimension is %d\n", P4EST_DIM);
    sc_options_print_summary (p4est_package_id, SC_LP_PRODUCTION, opt);

    /* check consistency of parameters */
    if (g->minlevel < 0 || g->minlevel > P4EST_QMAXLEVEL) {
      ue = usagerrf (opt, "minlevel not between 0 and %d", P4EST_QMAXLEVEL);
    }
    if (g->maxlevel < g->minlevel || g->maxlevel > P4EST_QMAXLEVEL) {
      ue = usagerrf (opt, "maxlevel not between minlevel and %d",
                     P4EST_QMAXLEVEL);
    }
    if (g->synthetic >= 0 ? 
          (g->latlongno >= 0 || g->sphere == 1) 
        : (g->latlongno >= 0 && g->sphere == 1)
       ) {
      ue = usagerrf (opt, "set only one of the synthetic, sphere and latlong models");
    }
    if (g->synthetic < 0 && g->latlongno < 0 && g->sphere == 0) {
      ue = usagerrf (opt, "set one of the synthetic, sphere, and latlong models");
    }
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }

  /* setup appplication model */
  if (!ue && setup_model (g)) {
    P4EST_ASSERT (g->model == NULL);
    ue = usagerr (opt, "model-specific initialization error");
  }

  /* execute application model */
  if (!ue && g->distributed == 0) {
    P4EST_ASSERT (g->model != NULL);
    run_program (g);
  }

  /* execute application model */
  if (!ue && g->distributed == 1) {
    P4EST_ASSERT (g->model != NULL);
    run_program_distributed (g);
  }

  /* cleanup application model */
  if (g->model != NULL) {
    p4est_connectivity_destroy (g->model->conn);
    p4est_gmt_model_destroy (g->model);
  }

  /* deinit main program */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return ue ? EXIT_FAILURE : EXIT_SUCCESS;
}
