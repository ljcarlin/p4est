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

#include "gmt_models.h"
#include <p4est_communication.h>
#include<p4est_base.h>

static void
model_set_geom(p4est_gmt_model_t *model,
               const char *name, p4est_geometry_X_t X)
{
  model->sgeom.name = name;
  // TODO: in p4est_geometry_connectivity_X user is assumed to point to
  // the relevant connectivity. The following line will cause geometries
  //  to violate that assumption.
  model->sgeom.user = model;
  model->sgeom.X = X;
  model->sgeom.destroy = NULL;
  model->model_geom = &model->sgeom;
}

/** Return the lowest rank owner. 
 * 
 * TODO this will result in an unbalanced distribution of responsibilities
 * and should eventually be replaced with hashing.
*/
static int model_resp_first (void *point, sc_array_t *owners) {
  int ret = *(int*)sc_array_index_int(owners, 0);
  return ret;
}

typedef struct p4est_gmt_model_synth
{
  int synthno;
  int resolution;
  size_t num_points;
  double *points;
} p4est_gmt_model_synth_t;

static void
model_synth_destroy_data(void *vmodel_data)
{
  p4est_gmt_model_synth_t *sdata = (p4est_gmt_model_synth_t *)vmodel_data;
  P4EST_FREE(sdata->points);
  P4EST_FREE(sdata);
}

static int
model_synth_intersect(p4est_topidx_t which_tree, const double coord[4],
                      size_t m, void *vmodel)
{
  p4est_gmt_model_t *model = (p4est_gmt_model_t *)vmodel;
  p4est_gmt_model_synth_t *sdata;
  const double *pco;
  double hx, hy;

  P4EST_ASSERT(model != NULL);
  P4EST_ASSERT(m < model->M);
  sdata = (p4est_gmt_model_synth_t *)model->model_data;
  P4EST_ASSERT(sdata != NULL && sdata->points != NULL);
  pco = sdata->points + 2 * m;
  P4EST_ASSERT(sdata->resolution >= 0);

  /* In this model we have only one tree, the unit square. */
  P4EST_ASSERT(which_tree == 0);

  /* Rectangle coordinates are in [0, 1] for the numbered reference tree and
   * stored as { lower left x, lower left y, upper right x, upper right y }. */

  /* We do not refine if target resolution is reached. */
  hx = coord[2] - coord[0];
  hy = coord[3] - coord[1];
  if (SC_MAX(hx, hy) <= pow(.5, sdata->resolution))
  {
    return 0;
  }

  /* In this synthetic example the point IS the object.  There are no lines. */
  if ((coord[0] <= pco[0] && pco[0] <= coord[2]) &&
      (coord[1] <= pco[1] && pco[1] <= coord[3]))
  {
    return 1;
  }

  /* We have exhausted the refinement criteria. */
  return 0;
}

static void
model_synth_geom_X(p4est_geometry_t *geom, p4est_topidx_t which_tree,
                   const double abc[3], double xyz[3])
{
  /* In this model we have only one tree, the unit square. */
  P4EST_ASSERT(which_tree == 0);

  /* We work with the unit square as physical space. */
  memcpy(xyz, abc, 3 * sizeof(double));
}

p4est_gmt_model_t *
p4est_gmt_model_synth_new(int synthno, int resolution)
{
  p4est_gmt_model_t *model = P4EST_ALLOC_ZERO(p4est_gmt_model_t, 1);
  p4est_gmt_model_synth_t *sdata = NULL;
  double *p;

  /* initalize model */
  switch (synthno)
  {
  case 0:
    model->output_prefix = "triangle";
    model->conn = p4est_connectivity_new_unitsquare();
    model->model_data = sdata = P4EST_ALLOC(p4est_gmt_model_synth_t, 1);
    sdata->synthno = synthno;
    sdata->resolution = resolution;
    sdata->num_points = model->M = 3;
    p = sdata->points = P4EST_ALLOC(double, 6);
    p[0] = 0.2;
    p[1] = 0.1;
    p[2] = 0.7;
    p[3] = 0.4;
    p[4] = 0.5;
    p[5] = 0.8;
    model->destroy_data = model_synth_destroy_data;
    model->intersect = model_synth_intersect;
    model_set_geom(model, model->output_prefix, model_synth_geom_X);
    break;
    /* possibly add more cases that work with polygon segments */
  default:
    SC_ABORT_NOT_REACHED();
  }

  /* return initialized model */
  P4EST_ASSERT(sdata != NULL);
  sdata->synthno = synthno;
  return model;
}

static int
model_latlong_intersect(p4est_topidx_t which_tree, const double coord[4],
                        size_t m, void *vmodel)
{
  p4est_gmt_model_t *model = (p4est_gmt_model_t *)vmodel;

  P4EST_ASSERT(model != NULL);
  P4EST_ASSERT(m < model->M);

  /* Rectangle coordinates are in [0, 1] for the numbered reference tree and
   * stored as { lower left x, lower left y, upper right x, upper right y }. */

  return 0;
}

static void
model_latlong_geom_X(p4est_geometry_t *geom, p4est_topidx_t which_tree,
                     const double abc[3], double xyz[3])
{
#if LATLONG_DATA_HAS_BEEN_PROGRAMMED
  p4est_gmt_model_t *model = (p4est_gmt_model_t *)geom->user;

  /* put the parameters latitude, longitude into the model data */
  longitude =
      ((typecast into gmt lanlong model data *)model->model_data)->longitude;
  latitude = ...;

  xyz[0] = longitude[0] + (longitude->[1] - longitude[0]) * abc[0];
  xyz[1] = latitude[0] + (latitude->[1] - latitude[0]) * abc[1];
#else
  xyz[0] = abc[0];
  xyz[1] = abc[1];
#endif
  xyz[2] = 0.;
}

p4est_gmt_model_t *
p4est_gmt_model_latlong_new(p4est_gmt_model_latlong_params_t *params)
{
  p4est_gmt_model_t *model = P4EST_ALLOC_ZERO(p4est_gmt_model_t, 1);

  /* the latlong models live on the unit square as reference domain */
  model->conn = p4est_connectivity_new_unitsquare();

  /* load model properties */
  model->model_data = NULL; /* <- Load something from params->load_filename,
                               also deep copy the parameters into it. */

  /* set virtual functions */
  model->intersect = model_latlong_intersect;
  model->destroy_data = NULL; /* <- needs to free whatever is in model_data */

  /* setup input/output parameters */
  model->output_prefix = params->output_prefix;
  model_set_geom(model, params->output_prefix, model_latlong_geom_X);

  /* the model is ready */
  model->M = 17; /* <- update to actual value */
  return model;
}

/** Represents the intersection of a geodesic with one of the cubed
 * connectivity faces.
 *
 * The endpoints p1 and p2 are given in tree-local coordinates and
 * in local coordinates the geodesic is just the line segment between
 * them.
 */

static void
model_sphere_destroy_data(void *vmodel_data)
{
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)vmodel_data;
  P4EST_FREE(sdata->points);
  P4EST_FREE(sdata);
}

/** Returns 1 if the line segments (p0 to p1) and (p2 to p3) intersect, otherwise 0 */
static int lines_intersect(double p0_x, double p0_y, double p1_x, double p1_y,
                           double p2_x, double p2_y, double p3_x, double p3_y)
{
  /* We solve the matrix equation (p1-p0, p2-p3) (s, t)^T = (p2-p0),
   * by inverting the matrix (p1-p0, p2-p3). */

  /* Precompute reused values for efficiency */
  double s1_x, s1_y, s2_x, s2_y;
  s1_x = p1_x - p0_x;
  s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;
  s2_y = p3_y - p2_y;

  /* Compute line intersection */
  double s, t;
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) / (-s2_x * s1_y + s1_x * s2_y);
  t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) / (-s2_x * s1_y + s1_x * s2_y);

  /* Check intersection lies on relevant segment */
  if (s >= 0 && s <= 1 && t >= 0 && t <= 1)
  {
    return 1;
  }

  return 0;
}

/** Returns 1 if the given geodesic intersects the given rectangle and 0 otherwise.
 *
 * \param[in] which_tree  tree id inside forest
 * \param[in] coord       rectangle for intersection checking. Rectangle coordinates
 *                        are in [0, 1] for the numbered reference tree and stored as
 *                        { lower left x, lower left y, upper right x, upper right y }.
 * \param[in] m           index of the geodesic we are checking
 * \param[in] vmodel      spherical model
 *
 */
static int
model_sphere_intersect(p4est_topidx_t which_tree, const double coord[4],
                       size_t m, void *vmodel)
{
  p4est_gmt_model_t *model = (p4est_gmt_model_t *)vmodel;
  p4est_gmt_model_sphere_t *sdata;
  const p4est_gmt_sphere_geodesic_seg_t *pco; /* mth geodesic segment */
  double hx, hy;                        /* width, height */

  P4EST_ASSERT(model != NULL);
  P4EST_ASSERT(m < model->M);
  sdata = (p4est_gmt_model_sphere_t *)model->model_data;
  P4EST_ASSERT(sdata != NULL && sdata->points != NULL);
  pco = sdata->points + m;
  P4EST_ASSERT(sdata->resolution >= 0);

  /* In this model we have 6 trees */
  P4EST_ASSERT(which_tree >= 0 && which_tree <= 5);

  /* Check the segment is on the relevant tree */
  if (pco->which_tree != which_tree)
  {
    return 0;
  }

  /* We do not refine if target resolution is reached. */
  hx = coord[2] - coord[0];
  hy = coord[3] - coord[1];
  if (SC_MAX(hx, hy) <= pow(.5, sdata->resolution))
  {
    return 0;
  }

  /* Check if the line segment L between p1 and p2 intersects the edges of
   * the rectangle.
   */

  /* Check if L intersects the bottom edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[1], coord[2], coord[1]))
  {
    return 1;
  }
  /* Check if L intersects the top edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[3], coord[2], coord[3]))
  {
    return 1;
  }
  /* Check if L intersects the left edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[0], coord[1], coord[0], coord[3]))
  {
    return 1;
  }
  /* Check if L intersects the right edge of rectangle */
  if (lines_intersect(pco->p1x, pco->p1y, pco->p2x, pco->p2y, coord[2], coord[1], coord[2], coord[3]))
  {
    return 1;
  }

  /* Check if L is contained in the interior of rectangle.
   * Since we have already ruled out intersections it suffices
   * to check if one of the endpoints of L is in the interior.
   */
  if (pco->p1x >= coord[0] && pco->p1x <= coord[2] && pco->p1y >= coord[1] && pco->p1y <= coord[3])
  {
    return 1;
  }

  /* We have exhausted the refinement criteria. */
  return 0;
}

/** Compute the processes whose domain potentially intersects a geodesic segment.
 * 
 *  We do this by taking the smallest quadrant containing the segment, and then
 *  checking which processes own the first and last atom of this quadrant.
 */
static sc_array_t* 
model_sphere_owners(void *point, p4est_t * p4est) 
{
  int first, last;
  p4est_gmt_sphere_geodesic_seg_t *seg = (p4est_gmt_sphere_geodesic_seg_t*)point;
  p4est_quadrant_t q;
  sc_array_t *owners;

  // /* DEBUG */

  // owners = sc_array_new_count(sizeof(int), p4est->mpisize);

  // for (int i = 0; i < p4est->mpisize; i++) {
  //   *(int*)sc_array_index_int(owners, i) = i;
  // }

  // return owners;

  // /* END DEBUG */

  /* First atom */
  q.p.which_tree = seg->which_tree;
  q.level = P4EST_MAXLEVEL;
  q.x = seg->bb1x;
  q.y = seg->bb1y;
  first = p4est_comm_find_owner(p4est, seg->which_tree, &q, 0); //TODO: Better guess?
  printf("First owner: %d\n\n", first);

  /* Last atom */
  q.x = seg->bb2x;
  q.y = seg->bb2y;
  last = p4est_comm_find_owner(p4est, seg->which_tree, &q, 0); //TODO: Better guess?
  printf("Last owner: %d\n\n", first);

  owners = sc_array_new_count(sizeof(int), last-first+1);

  printf("OWNERS: ");
  for (int i = 0; i <= last-first; i++) {
    *(int*)sc_array_index_int(owners, i) = first+i;
    printf("%d, ", first+i);
  }
  printf("\n");
  
  return owners;
}

p4est_gmt_model_t *
p4est_gmt_model_sphere_new(int resolution, int dist, sc_MPI_Comm mpicomm)
{
  //FILE *geodesic_file;
  sc_MPI_File file_handle;
  p4est_gmt_model_t *model = P4EST_ALLOC_ZERO(p4est_gmt_model_t, 1);
  p4est_gmt_model_sphere_t *sdata = NULL;
  int rank, num_procs;
  int mpiret;
  int global_num_points, local_num_points;
  int count;
  p4est_gloidx_t offset_mine, offset_next;
  sc_MPI_Offset mpi_offset;

  model->model_data = sdata = P4EST_ALLOC(p4est_gmt_model_sphere_t, 1);

  /* Get rank and number of processes */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* Open file of precomputed geodesic segments */
  mpiret = sc_io_open (mpicomm, "geodesics", SC_IO_READ, sc_MPI_INFO_NULL,
                       &file_handle);
  SC_CHECK_MPI (mpiret);

  if (rank == 0) {
    /* read the global number of points from file */
    mpiret = sc_io_read_at (file_handle, 0, &global_num_points,
                            sizeof (int), sc_MPI_BYTE, &count);
    SC_CHECK_MPI (mpiret);
    SC_CHECK_ABORT (count == (int) sizeof (int),
                    "Read number of global points: count mismatch");
  }

  /* broadcast the global number of points */
  mpiret = sc_MPI_Bcast (&global_num_points, sizeof (int),
                          sc_MPI_BYTE, 0, mpicomm);
  SC_CHECK_MPI (mpiret);

  /* set read offsets depending on whether we are running distributed */
  if (dist==0) {
    mpi_offset = 0;
    local_num_points = global_num_points;
  }
  else {
    /* offset to first point of current MPI process */
    offset_mine = p4est_partition_cut_gloidx (global_num_points,
                                            rank, num_procs);

    /* offset to first point of successor MPI process */
    offset_next = p4est_partition_cut_gloidx (global_num_points,
                                            rank + 1, num_procs);

    /* set file offset (in bytes) for this calling process */
    mpi_offset = (sc_MPI_Offset) offset_mine * sizeof(p4est_gmt_sphere_geodesic_seg_t);

    local_num_points = offset_next - offset_mine;
  }

  /* allocate buffer for geodesic segments */
  sdata->points = P4EST_ALLOC(p4est_gmt_sphere_geodesic_seg_t, local_num_points);

  /* TODO: is this right in non-distributed mode? */
  /* each mpi process reads its data for its own offset */
  mpiret = sc_io_read_at_all (file_handle, mpi_offset + sizeof (size_t),
                            &(sdata->points[0]), 
                            local_num_points * sizeof (p4est_gmt_sphere_geodesic_seg_t),
                            sc_MPI_BYTE, &count);
  SC_CHECK_MPI (mpiret);
  SC_CHECK_ABORT (count == (int) (local_num_points 
                                    * sizeof (p4est_gmt_sphere_geodesic_seg_t)),
                    "Read points: count mismatch");

  /* close the file collectively */
  mpiret = sc_io_close (&file_handle);
  SC_CHECK_MPI (mpiret);

  /* Set final geodesic count */
  sdata->num_geodesics = model->M = local_num_points; 

  /* Set point size */
  model->point_size = sizeof(p4est_gmt_sphere_geodesic_seg_t);

  /* the sphere model lives on the cube surface reference */
  model->conn = p4est_connectivity_new_cubed();
  model->output_prefix = "sphere";

  /* Assign resolution, intersector, and destructor */
  sdata->resolution = resolution;
  model->destroy_data = model_sphere_destroy_data;
  model->intersect = model_sphere_intersect;

  /* Assign owner and responsibility functions */
  model->owners_fn = model_sphere_owners;
  model->resp_fn = model_resp_first;

  /* Assign geometry */
  model->model_geom = p4est_geometry_new_sphere2d(model->conn, 1.0);

  /* the model is ready */
  return model;
}

void p4est_gmt_model_destroy(p4est_gmt_model_t *model)
{
  if (model->destroy_data != NULL)
  {
    /* only needed for non-trivial free code */
    model->destroy_data(model->model_data);
  }
  else
  {
    /* the default clears a standard allocation or respects NULL */
    P4EST_FREE(model->model_data);
  }
  p4est_geometry_destroy(model->model_geom);
  P4EST_FREE(model);
}