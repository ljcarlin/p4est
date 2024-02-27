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
#include<sc_notify.h>
#include "gmt_global.h"
#include <p4est_search.h>

static const double irootlen = 1. / (double) P4EST_ROOT_LEN;

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
  double s1_x, s1_y, s2_x, s2_y, det_inv;
  s1_x = p1_x - p0_x;
  s1_y = p1_y - p0_y;
  s2_x = p3_x - p2_x;
  s2_y = p3_y - p2_y;
  det_inv = -s2_x * s1_y + s1_x * s2_y;

  /* Compute line intersection */
  double s, t;
  s = (-s1_y * (p0_x - p2_x) + s1_x * (p0_y - p2_y)) * det_inv;
  t = (s2_x * (p0_y - p2_y) - s2_y * (p0_x - p2_x)) * det_inv;

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

  /* we do not refine if target resolution is reached */
  hx = coord[2] - coord[0];
  hy = coord[3] - coord[1];
  if (SC_MAX(hx, hy) <= pow(.5, sdata->resolution)) {
    return 0;
  }

  /* cull obvious non-intersections before running more expensive line
   * intersection computations
   */
  /* check if segment is left of quadrant */
  if (pco->p1x < coord[0] && pco->p2x < coord[0]) {
    return 0;
  }
  /* check if segment is below quadrant */
  if (pco->p1y < coord[1] && pco->p2y < coord[1]) {
    return 0;
  }
  /* check if segment is right of quadrant */
  if (pco->p1x > coord[2] && pco->p2x > coord[2]) {
    return 0;
  }
  /* check if segment is above quadrant */
  if (pco->p1y > coord[3] && pco->p2y > coord[3]) {
    return 0;
  }

  /* Check if segment is contained in the interior of quadrant */
  if (pco->p1x >= coord[0] && pco->p1x <= coord[2] 
      && pco->p1y >= coord[1] && pco->p1y <= coord[3]
      && pco->p2x >= coord[0] && pco->p2x <= coord[2] 
      && pco->p2y >= coord[1] && pco->p2y <= coord[3]) {
    return 1;
  }

  /* Check if the segment intersects the edges of quadrant */
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

  /* First atom */
  q.p.which_tree = seg->which_tree;
  q.level = P4EST_QMAXLEVEL;
  q.x = seg->bb1x;
  q.y = seg->bb1y;
  first = p4est_comm_find_owner(p4est, seg->which_tree, &q, p4est->mpirank);

  /* Last atom */
  q.x = seg->bb2x;
  q.y = seg->bb2y;
  last = p4est_comm_find_owner(p4est, seg->which_tree, &q, p4est->mpirank);

  owners = sc_array_new_count(sizeof(int), last-first+1);
  for (int i = 0; i <= last-first; i++) {
    *(int*)sc_array_index_int(owners, i) = first+i;
  }

  return owners;
}

/** Prepare outgoing buffers of points to propagate.
 * 
 * \param[out] resp Comm data for points whose receiver will be responsible
 *                  for propagating them in the next iteration
 * \param[out] own  Comm data for points whose receiver will *not* be
 *                  responsible for propagating in the next iteration
 * \param[in] p4est The forest
 * \param[in] model Sphere model
 * \param[in] num_procs Number of MPI processes
 */
static void 
sphere_compute_outgoing_points(p4est_gmt_sphere_comm_t *resp,
                                p4est_gmt_sphere_comm_t *own,
                                p4est_t *p4est,
                                p4est_gmt_model_t *model,
                                int num_procs
                              )
{
  /* Sphere model */
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)model->model_data;
  /** per point data **/
  sc_array_t         *owners; /* who to send the point to */
  int                 responsible; /* who is in charge of propagating the point */

  /* Initialise index of outgoing message buffers */
  own->to_send = P4EST_ALLOC(sc_array_t*, num_procs);
  resp->to_send = P4EST_ALLOC(sc_array_t*, num_procs);

  /* Initialise outgoing message buffers */
  /* TODO: It might make sense to only initialize when we are actually sending */
  for (int q = 0; q < num_procs; q++) {
    own->to_send[q] = sc_array_new(model->point_size);
    resp->to_send[q] = sc_array_new(model->point_size);
  }

  /* Iterate through the points we are responsible for propagating */
  for (int i = 0; i < (int)sdata->num_owned_resp; i++) {

    /* Determine which processes should receive this point */
    owners = model->owners_fn(&(sdata->points[i]), p4est);

    /* Determine the receiver responsible for propagating this point in next
        iteration */
    responsible = model->resp_fn(&(sdata->points[i]), owners);
    P4EST_ASSERT(0 <= responsible && responsible < num_procs);

    /* Add this point to the relevant outgoing message buffers */
    for (size_t j = 0; j < owners->elem_count; j++) {
      int q = *(int*)sc_array_index_int(owners, j);
      if (responsible == q)
      {
        /* Process q should own point i and be responsible for its propagation */
        *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_push(resp->to_send[q]) = sdata->points[i];
      }
      else 
      {
        /* Process q should own point i but not be responsible for its propagation */
        *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_push(own->to_send[q]) = sdata->points[i];
      }
    }

    /* clean up to prepare for next point*/
    sc_array_destroy(owners);
  }
}

/* quadrant callback for search_partition */
static int
sphere_partition_quadrant (p4est_t * p4est, p4est_topidx_t which_tree,
                            p4est_quadrant_t * quadrant, int pfirst,
                            int plast, void *point)
{

  P4EST_ASSERT (0 <= pfirst && pfirst <= plast);
  P4EST_ASSERT (point == NULL);

  /* always continue, as recursion is controlled by point callbacks */
  return 1;
}

/* point callback for search_partition */
static int
sphere_partition_point (p4est_t * p4est, p4est_topidx_t which_tree,
                        p4est_quadrant_t * quadrant, int pfirst, int plast,
                        void *point)
{
  /* global context */
  global_t           *g = (global_t *) p4est->user_pointer;
  /* sphere model */
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)g->model->model_data;
  /* last process's send buffer this point was added to */
  int                last_proc;
  /* quadrant coords*/
  double              coord[4];
  p4est_qcoord_t      qh;
  size_t pi = *(size_t *) point;

  /* sanity checks */
  P4EST_ASSERT (g != NULL && g->p4est == p4est);
  P4EST_ASSERT (0 <= pfirst && pfirst <= plast);
  P4EST_ASSERT (point != NULL);
  P4EST_ASSERT (pi < sdata->num_owned_resp);

  /* quadrant coordinates */
  qh = P4EST_QUADRANT_LEN (quadrant->level);
  coord[0] = irootlen * quadrant->x;
  coord[1] = irootlen * quadrant->y;
  coord[2] = irootlen * (quadrant->x + qh);
  coord[3] = irootlen * (quadrant->y + qh);

  /* if current quadrant has multiple owners */
  if (pfirst < plast) {
    /* point follows recursion when it intersects the quadrant */
    /* TODO: bounding box heuristic */
    return model_sphere_intersect (quadrant->p.which_tree, coord, pi, g->model);
  }

  /* current quadrant has a single owner */
  P4EST_ASSERT (pfirst == plast);

  if (!model_sphere_intersect (quadrant->p.which_tree, coord, pi, g->model)) {
    /* point does not intersect this quadrant */
    return 0;
  }

  /* get last process whose domain we have already recorded as intersecting 
   * this point
   */
  last_proc = *(int *) sc_array_index (sdata->geoseg_procs, pi);
  
  /* since we traverse in order we expect not to have seen this point in
   * in higher process domains yet
   */
  P4EST_ASSERT (last_proc <= pfirst);

  if (last_proc == pfirst) {
    /* we have found an already recorded process */
    return 0;
  }
  /* otherwise we have found a new process intersecting the point */

  /* record this new process */
  *(int *)sc_array_index (sdata->geoseg_procs, pi) = pfirst;

  /* add point to corresponding send buffer */
  if (last_proc == -1) {
    /* process should own point and be responsible for its propagation */
    *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_push (sdata->resp.to_send[pfirst]) = sdata->points[pi];
  }
  else {
    /* process should own point but not be responsible for its propagation */
    *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_push (sdata->own.to_send[pfirst]) = sdata->points[pi];
  }

  return 0;
}

/** an alternative version of sphere_compute_outgoing_points using 
 * search_partition 
 */
static void 
sphere_compute_outgoing_points_alt(p4est_gmt_sphere_comm_t *resp,
                                    p4est_gmt_sphere_comm_t *own,
                                    p4est_t *p4est,
                                    p4est_gmt_model_t *model,
                                    int num_procs
                                  )
{

  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)model->model_data;
  sc_array_t               *points;

  /* initialise to -1 to signify no points have been added to send buffers */
  sdata->geoseg_procs = sc_array_new_count(sizeof (int), sdata->num_owned_resp);
  sc_array_memset (sdata->geoseg_procs, -1);

  /* initialise index of outgoing message buffers */
  own->to_send = P4EST_ALLOC(sc_array_t*, num_procs);
  resp->to_send = P4EST_ALLOC(sc_array_t*, num_procs);

  /* initialise outgoing message buffers */
  /* TODO: It might make sense to only initialize when we are actually sending */
  for (int q = 0; q < num_procs; q++) {
    own->to_send[q] = sc_array_new(model->point_size);
    resp->to_send[q] = sc_array_new(model->point_size);
  }

  /* set up search objects for partition search */
  points = sc_array_new_count (sizeof (size_t), sdata->num_owned_resp);
  for (size_t zz = 0; zz < sdata->num_owned_resp; ++zz) {
    *(size_t *) sc_array_index (points, zz) = zz;
  }

  /* add points to the relevant send buffers (by partition search) */
  p4est_search_partition (p4est, 0, sphere_partition_quadrant,
                          sphere_partition_point, points);

  /* clean up search objects */
  sc_array_destroy_null (&points);
  sc_array_destroy_null (&sdata->geoseg_procs);
}

/** Update communication data with which processes p is sending points to and
 *  how many points are being sent to each of these. 
 * 
 *  The output is stored in the fields comm->receivers and comm->recvs_counts.
 *  We assume comm->to_send is already populated.
 * 
 * \param[in,out] comm communication data.
 * \param[in] num_procs number of mpi processes
 */
static void sphere_compute_receivers(p4est_gmt_sphere_comm_t *comm, int num_procs) 
{
  /* initialize receivers and counts*/
  comm->receivers = sc_array_new(sizeof(int));
  comm->recvs_counts = sc_array_new(sizeof(size_t));

  /* compute receivers and counts */
  for (int q = 0; q < num_procs; q++) {
    if (comm->to_send[q]->elem_count != 0) {
      /* add q to receivers */
      *(int*)sc_array_push(comm->receivers) = q;
      /* record how many points p is sending to q */
      *(size_t*)sc_array_push(comm->recvs_counts) = comm->to_send[q]->elem_count;
    }
  }
  P4EST_ASSERT(comm->receivers->elem_count == comm->recvs_counts->elem_count);
}

/** Update communication data with total number of incoming points, and 
 *  offsets to receive incoming points at. 
 * 
 *  The outputs are stored in the fields comm->num_incoming and comm->offsets.
 *  We assume that comm->senders and comm->senders_counts are already
 *  populated.
 * 
 * \param[in,out] comm communication data.
 */
static void sphere_compute_offsets(p4est_gmt_sphere_comm_t *comm) 
{
  /* initialize offset array */
  comm->num_incoming = 0;
  comm->offsets = P4EST_ALLOC(size_t, comm->senders->elem_count);

  /* compute offsets */
  for (int i = 0; i < (int)comm->senders->elem_count; i++) {
    comm->offsets[i] = comm->num_incoming;
    comm->num_incoming += *(size_t*)sc_array_index_int(comm->senders_counts, i);
  } 
}

/** Post non-blocking sends for points in the given communication data.
 * 
 * To each rank q in comm->receivers we send the points stored at 
 * comm->to_send[q]
 * 
 * \param[in]   comm        communication data
 * \param[in]   mpicomm     MPI communicator
 * \param[out]  req         request storage of same length as comm->receivers
 * \param[in]   point_size  size of a single point
 */
static void sphere_post_sends(p4est_gmt_sphere_comm_t *comm,
                              sc_MPI_Comm mpicomm,
                              sc_MPI_Request * req,
                              size_t point_size
                             ) 
{
  int mpiret;

  for (int i = 0; i < (int)comm->receivers->elem_count; i++) {
    int q = *(int*)sc_array_index_int(comm->receivers, i);
    mpiret = sc_MPI_Isend(sc_array_index_int(comm->to_send[q], 0), /* out buffer */
                 comm->to_send[q]->elem_count * point_size, /* num. of bytes */
                 sc_MPI_BYTE,
                 q, /* receiving process */
                 0, /* arbitrary tag */
                 mpicomm,
                 req + i /* request */
                 );
    SC_CHECK_MPI (mpiret);
  }
}

/** Post non-blocking receives for senders in the given communication data.
 * 
 *  We expect to receive points from each sender in comm->senders. The number
 *  of points each sender is sending is stored in comm->senders_counts (with
 *  corresponding indexing). We receive each message at the offset stored in
 *  comm->offsets (again with corresponding indexing).
 * 
 *  \param[in] comm communication data
 *  \param[in,out] buffer points to array where received points are stored
 *  \param[out] req request storage of same length as comm->senders
 *  \param[in] mpicomm MPI communicator
 *  \param[in] point_size size of a single point
*/
static void sphere_post_receives(p4est_gmt_sphere_comm_t *comm,
                                 p4est_gmt_sphere_geodesic_seg_t *buffer,
                                 sc_MPI_Request *req,
                                 sc_MPI_Comm mpicomm,
                                 size_t point_size
                                )
{
  int mpiret;

  for (int i = 0; i < (int)comm->senders->elem_count; i++) {
    int q = *(int*)sc_array_index_int(comm->senders, i);
    mpiret = sc_MPI_Irecv(buffer + comm->offsets[i],
                (*(size_t*)sc_array_index_int(comm->senders_counts, i)) * point_size,
                sc_MPI_BYTE,
                q,  /* sending process */
                0, /* arbitrary tag */
                mpicomm,
                req + i /* request */
    );
    SC_CHECK_MPI (mpiret);
  }
}

/** Communicate geodesics to the processes whose domain they *may* overlap.
 * 
 * \param[in] mpicomm   MPI communicator
 * \param[in] p4est     The forest
 * \param[in] model     Sphere model
 */
static void
model_sphere_communicate_points(sc_MPI_Comm mpicomm,
                                    p4est_t *p4est,
                                    p4est_gmt_model_t *model) 
{
  /* Sphere model */
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)model->model_data;
  /* MPI data */
  int                 mpiret;
  int                 num_procs;
  p4est_gmt_sphere_comm_t *own = &sdata->own; /* not responsible*/
  p4est_gmt_sphere_comm_t *resp = &sdata->resp; /* responsible */

  /* TODO: This naming convention is kind of confusing. Why not outgoing/incoming?*/
  sc_MPI_Request *recv_req; /* Requests for sends to receivers */
  sc_MPI_Request *sender_req; /* Requests for receives from senders */

  /* get total process count */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);

  /* prepare outgoing buffers of points to send */
  sphere_compute_outgoing_points_alt(resp, own, p4est, model, num_procs);

  /* free the points we received last iteration */
  P4EST_FREE(sdata->points);

  /* record which processes p is sending points to and how many points each
     process receives */
  sphere_compute_receivers(own, num_procs);
  sphere_compute_receivers(resp, num_procs);

  /* initialise outgoing request arrays */
  recv_req = P4EST_ALLOC(sc_MPI_Request, own->receivers->elem_count + resp->receivers->elem_count);

  /* post non-blocking sends */
  sphere_post_sends(resp, mpicomm, recv_req, model->point_size);
  sphere_post_sends(own, mpicomm, recv_req + resp->receivers->elem_count, model->point_size);

  /* initialize buffers for receiving communication data with sc_notify_ext */
  own->senders = sc_array_new(sizeof(int));
  resp->senders = sc_array_new(sizeof(int));
  own->senders_counts = sc_array_new(sizeof(size_t));
  resp->senders_counts = sc_array_new(sizeof(size_t));

  /* notify processes receiving points from p and determine processes sending
   to p. Also exchange counts of points being sent. */
  sc_notify_ext(own->receivers, own->senders, own->recvs_counts, own->senders_counts, mpicomm);
  sc_notify_ext(resp->receivers, resp->senders, resp->recvs_counts, resp->senders_counts, mpicomm);
  P4EST_ASSERT(own->senders->elem_count == own->senders_counts->elem_count);
  P4EST_ASSERT(resp->senders->elem_count == resp->senders_counts->elem_count);

  /* compute offsets for storing incoming points */
  sphere_compute_offsets(own);
  sphere_compute_offsets(resp);

  /* update counts of known points */
  model->M = resp->num_incoming + own->num_incoming;
  sdata->num_owned_resp = resp->num_incoming;
  sdata->num_owned = own->num_incoming;

  /* allocate memory for incoming points */
  sdata->points = P4EST_ALLOC(p4est_gmt_sphere_geodesic_seg_t, model->M);

  /* initialise incoming request arrays */
  sender_req = P4EST_ALLOC(sc_MPI_Request, own->senders->elem_count + resp->senders->elem_count);

  /* post non-blocking receives */
  sphere_post_receives(resp, sdata->points,
                        sender_req,
                        mpicomm, model->point_size);
  sphere_post_receives(own, sdata->points + resp->num_incoming, 
                        sender_req + resp->senders->elem_count, 
                        mpicomm, model->point_size);

  /* wait for messages to send */
  mpiret = sc_MPI_Waitall (own->receivers->elem_count + resp->receivers->elem_count, 
      recv_req,
      sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up data used for sending */
  sc_array_destroy(own->receivers);
  sc_array_destroy(resp->receivers);
  sc_array_destroy(own->recvs_counts);
  sc_array_destroy(resp->recvs_counts);
  for (int q = 0; q < num_procs; q++) {
    sc_array_destroy(own->to_send[q]);
    sc_array_destroy(resp->to_send[q]);
  }
  P4EST_FREE(own->to_send);
  P4EST_FREE(resp->to_send); 
  P4EST_FREE(own->offsets);
  P4EST_FREE(resp->offsets);
  P4EST_FREE(recv_req);

  /* Wait to receive messages */
  mpiret = sc_MPI_Waitall (own->senders->elem_count + resp->senders->elem_count, 
      sender_req,
      sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up data used for receiving */
  sc_array_destroy(own->senders); 
  sc_array_destroy(resp->senders); 
  sc_array_destroy(own->senders_counts);
  sc_array_destroy(resp->senders_counts);
  P4EST_FREE(sender_req);
}

p4est_gmt_model_t *
p4est_gmt_model_sphere_new(int resolution, int dist, sc_MPI_Comm mpicomm)
{
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
  if (!dist) {
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
  sdata->num_owned_resp = model->M = local_num_points;
  /* We are responsible for propagating all points that we intially read */
  sdata->num_owned = 0; 

  /* Set point size */
  model->point_size = sizeof(p4est_gmt_sphere_geodesic_seg_t);

  /* the sphere model lives on the cube surface reference */
  model->conn = p4est_connectivity_new_cubed();
  model->output_prefix = "sphere";

  /* Assign functions */
  sdata->resolution = resolution;
  model->destroy_data = model_sphere_destroy_data;
  model->intersect = model_sphere_intersect;
  if (dist) {
    model->communicate_points = model_sphere_communicate_points;  
  }

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
