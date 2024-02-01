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
  q.level = P4EST_QMAXLEVEL;
  q.x = seg->bb1x;
  q.y = seg->bb1y;
  first = p4est_comm_find_owner(p4est, seg->which_tree, &q, 0); //TODO: Better guess?
  //printf("First owner: %d\n\n", first);

  /* Last atom */
  q.x = seg->bb2x;
  q.y = seg->bb2y;
  last = p4est_comm_find_owner(p4est, seg->which_tree, &q, 0); //TODO: Better guess?
  //printf("Last owner: %d\n\n", first);

  owners = sc_array_new_count(sizeof(int), last-first+1);
  for (int i = 0; i <= last-first; i++) {
    *(int*)sc_array_index_int(owners, i) = first+i;
  }

  return owners;
}

/*TODO: Most of the code gets duplicated. We could fix this by putting into functions?*/
/*Another option would be to make a struct for the duplicated stuff*/
/** Communicate geodesics to the processes whose domain they *may* overlap.
 * 
 * We iterate over all locally known geodesics, work out which domains they
 * may overlap, and then send them to these processes.
 * 
 * TODO: This is a WIP version using Gatherv. It will be much more efficient
 * with point-to-point communication. 
 * 
 * \param[in] mpicomm   MPI communicator
 * \param[in] p4est     The forest
 * \param[in] model     Sphere model
 */
static void
p4est_gmt_sphere_communicate_points(sc_MPI_Comm mpicomm,
                                    p4est_t *p4est,
                                    p4est_gmt_model_t *model) 
{
  /* Sphere model */
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)model->model_data;

  /* MPI data */
  int                 mpiret;
  int                 num_procs, rank;

  /** communication data. The local process is referred to as p **/
  /* number of points that p receives in this iteration (responsible)*/
  int                 num_incoming_resp; 
  /* number of points that p receives in this iteration (not responsible)*/
  int                 num_incoming;
  /* q -> {points q should own in the next iteration, and be responsible for} */
  sc_array_t        **to_send_resp; 
  /* q -> {points q should own in the next iteration, but not be responsible for } */
  sc_array_t        **to_send;

  sc_array_t *receivers; /* Ranks receiving points from p */
  sc_array_t *receivers_resp; /* Ranks receiving points from p */
  sc_array_t *senders; /* Ranks sending points to p */
  sc_array_t *senders_resp; /* Ranks sending points to p */
  sc_array_t *recvs_counts; /* Number of points each receiver gets from p */
  sc_array_t *recvs_counts_resp; /* Number of points each receiver gets from p */
  sc_array_t *senders_counts; /* Number of points p gets from each sender */
  sc_array_t *senders_counts_resp; /* Number of points p gets from each sender */
  size_t * offsets; /* Offsets for incoming messages */
  size_t * offsets_resp; /* Offsets for incoming messages */

  /* TODO: This naming convention is kind of confusing. Why not outgoing/incoming?*/
  sc_array_t *recv_req; /* Requests for sends to receivers */
  sc_array_t *sender_req; /* Requests for receives from senders */

  /** per point data **/
  sc_array_t         *owners; /* who to send the point to */
  int                 responsible; /* who is in charge of propagating the point */

  /* Get rank and total process count */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* Initialise index of outgoing message buffers */
  to_send = P4EST_ALLOC(sc_array_t*, num_procs);
  to_send_resp = P4EST_ALLOC(sc_array_t*, num_procs);

  /* Initialise outgoing message buffers */
  /* TODO: It might make sense to only initialize when we are actually sending */
  for (int q = 0; q < num_procs; q++) {
    to_send[q] = sc_array_new(model->point_size);
    to_send_resp[q] = sc_array_new(model->point_size);
  }

  /* Iterate through points to add them to outgoing messages */
  for (int i = 0; i < (int)sdata->num_owned_resp; i++) {
    /* Determine which processes should receive this point */
    owners = model->owners_fn(&(sdata->points[i]), p4est);
    /* Determine the receiver responsible for propagating this point in next
        iteration */
    //printf("Proc %d responsible\n", rank);
    responsible = model->resp_fn(&(sdata->points[i]), owners);
    P4EST_ASSERT(0 <= responsible && responsible < num_procs);

    //printf("Proc %d add to receivers\n", rank);
    /* Add this point to the arrays of the processes who receive it */
    for (size_t j = 0; j < owners->elem_count; j++) {
      int q = *(int*)sc_array_index_int(owners, j);
      if (responsible == q)
      {
        /* Process q should own point i and be responsible for its propagation */
        sc_array_push(to_send_resp[q]);
        *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_index_int(to_send_resp[q], to_send_resp[q]->elem_count-1) = sdata->points[i];
      }
      else 
      {
        /* Process q should own point i but not be responsible for its propagation */
        sc_array_push(to_send[q]);
        *(p4est_gmt_sphere_geodesic_seg_t*)sc_array_index_int(to_send[q], to_send[q]->elem_count-1) = sdata->points[i];
      }
    }
    sc_array_destroy(owners);
  }

  /* Initialize receivers and payload arrays */
  receivers = sc_array_new(sizeof(int));
  receivers_resp = sc_array_new(sizeof(int));
  recvs_counts = sc_array_new(sizeof(size_t));
  recvs_counts_resp = sc_array_new(sizeof(size_t));
  /* Compute receivers and payloads */
  for (int q = 0; q < num_procs; q++) {
    if (to_send[q]->elem_count != 0) {
      /* add q to receivers */
      sc_array_push(receivers);
      *(int*)sc_array_index_int(receivers, receivers->elem_count-1) = q;
      /* record how many points p is sending to q */
      sc_array_push(recvs_counts);
      *(size_t*)sc_array_index_int(recvs_counts, recvs_counts->elem_count-1) 
          = to_send[q]->elem_count;
    }
    if (to_send_resp[q]->elem_count != 0) {
      /* add q to receivers */
      sc_array_push(receivers_resp);
      *(int*) sc_array_index_int(receivers_resp, receivers_resp->elem_count-1) = q;
      /* record how many points (resp) p is sending to q */
      sc_array_push(recvs_counts_resp);
      *(size_t*)sc_array_index_int(recvs_counts_resp, recvs_counts_resp->elem_count-1) 
          = to_send_resp[q]->elem_count;
    }
  }
  P4EST_ASSERT(receivers->elem_count == recvs_counts->elem_count);
  P4EST_ASSERT(receivers_resp->elem_count == recvs_counts_resp->elem_count);

  /* initialize incoming message buffers */
  senders = sc_array_new(sizeof(int));
  senders_resp = sc_array_new(sizeof(int));
  senders_counts = sc_array_new(sizeof(size_t));
  senders_counts_resp = sc_array_new(sizeof(size_t));

  /* notify processes receiving points from p and determine processes sending
   to p */
  sc_notify_ext(receivers, senders, recvs_counts, senders_counts, mpicomm);
  sc_notify_ext(receivers_resp, senders_resp, recvs_counts_resp, senders_counts_resp, mpicomm);

  P4EST_ASSERT(senders->elem_count == senders_counts->elem_count);
  P4EST_ASSERT(senders_resp->elem_count == senders_counts_resp->elem_count);

  /* determine the number of points p will receive */
  /* TODO: also determine the offsets */
  num_incoming = 0;
  num_incoming_resp = 0;
  offsets = P4EST_ALLOC(size_t, senders->elem_count);
  offsets_resp = P4EST_ALLOC(size_t, senders_resp->elem_count);
  for (int i = 0; i < (int)senders_counts->elem_count; i++) {
    offsets[i] = num_incoming;
    num_incoming += *(size_t*)sc_array_index_int(senders_counts, i);
  } 
  for (int i = 0; i < (int)senders_counts_resp->elem_count; i++) {
    offsets_resp[i] = num_incoming_resp;
    num_incoming_resp += *(size_t*)sc_array_index_int(senders_counts_resp, i);
  }

  /* Free the points we received last iteration */
  P4EST_FREE(sdata->points);
  /* Create array for incoming points and update counts */
  model->M = num_incoming_resp + num_incoming;
  sdata->points = P4EST_ALLOC(p4est_gmt_sphere_geodesic_seg_t, model->M);
  sdata->num_owned_resp = num_incoming_resp;
  sdata->num_owned = num_incoming;

  /* Initialise request arrays */
  recv_req = sc_array_new_count (sizeof (sc_MPI_Request), receivers->elem_count + receivers_resp->elem_count);

  /* Post non-blocking sends */
  for (int i = 0; i < (int)receivers_resp->elem_count; i++) {
    int q = *(int*)sc_array_index_int(receivers_resp, i);
    mpiret = sc_MPI_Isend(sc_array_index_int(to_send_resp[q], 0), /* out buffer */
                 to_send_resp[q]->elem_count * model->point_size, /* num. of bytes */
                 sc_MPI_BYTE,
                 q, /* receiving process */
                 0, /* arbitrary tag */
                 mpicomm,
                 (sc_MPI_Request *)sc_array_index_int(recv_req, i)
                 );
    SC_CHECK_MPI (mpiret);
  }
  for (int i = 0; i < (int)receivers->elem_count; i++) {
    int q = *(int*)sc_array_index_int(receivers, i);
    mpiret = sc_MPI_Isend(sc_array_index_int(to_send[q], 0), /* out buffer */
                 to_send[q]->elem_count * model->point_size, /* num. of bytes */
                 sc_MPI_BYTE,
                 q, /* receiving process */
                 0, /* arbitrary tag */
                 mpicomm,
                 (sc_MPI_Request *)sc_array_index_int(recv_req, i + receivers_resp->elem_count)
                 );
    SC_CHECK_MPI (mpiret);
  }

  /* Initialise request arrays */
  sender_req = sc_array_new_count (sizeof (sc_MPI_Request), senders->elem_count + senders_resp->elem_count);

  /* Post non-blocking receives */
  for (int i = 0; i < (int)senders_resp->elem_count; i++) {
    int q = *(int*)sc_array_index_int(senders_resp, i);
    mpiret = sc_MPI_Irecv(sdata->points + offsets_resp[i], /* TODO buffer offset */
                (*(size_t*)sc_array_index_int(senders_counts_resp, i)) * model->point_size,
                sc_MPI_BYTE,
                q,  /* sending process */
                0, /* arbitrary tag */
                mpicomm,
                (sc_MPI_Request *)sc_array_index_int(sender_req, i)
    );
    SC_CHECK_MPI (mpiret);
  }
  for (int i = 0; i < (int)senders->elem_count; i++) {
    int q = *(int*)sc_array_index_int(senders, i);
    mpiret = sc_MPI_Irecv(sdata->points + offsets[i] + num_incoming_resp,  /* TODO buffer offset */
                (*(size_t*)sc_array_index_int(senders_counts, i)) * model->point_size,
                sc_MPI_BYTE,
                q,  /* sending process */
                0, /* arbitrary tag */
                mpicomm,
                (sc_MPI_Request *)sc_array_index_int(sender_req, i + senders_resp->elem_count)
    );
    SC_CHECK_MPI (mpiret);
  }

  /* Wait for messages to send */
  mpiret = sc_MPI_Waitall (recv_req->elem_count, 
      (sc_MPI_Request *)sc_array_index_int (recv_req, 0),
      sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up data used for sending */
  sc_array_destroy(receivers);
  sc_array_destroy(receivers_resp);
  sc_array_destroy(recv_req);
  sc_array_destroy(recvs_counts);
  sc_array_destroy(recvs_counts_resp);
  for (int q = 0; q < num_procs; q++) {
    sc_array_destroy(to_send[q]);
    sc_array_destroy(to_send_resp[q]);
  }
  P4EST_FREE(to_send);
  P4EST_FREE(to_send_resp); 
  P4EST_FREE(offsets);
  P4EST_FREE(offsets_resp);

  /* Wait to receive messages */
  mpiret = sc_MPI_Waitall (sender_req->elem_count, 
      (sc_MPI_Request *)sc_array_index_int (sender_req, 0),
      sc_MPI_STATUSES_IGNORE);
  SC_CHECK_MPI (mpiret);

  /* clean up data used for receiving */
  sc_array_destroy(senders); 
  sc_array_destroy(senders_resp); 
  sc_array_destroy(senders_counts);
  sc_array_destroy(senders_counts_resp);
  sc_array_destroy(sender_req);
}

/** Communicate geodesics to the processes whose domain they *may* overlap.
 * 
 * We iterate over all locally known geodesics, work out which domains they
 * may overlap, and then send them to these processes.
 * 
 * TODO: This is a WIP version using Gatherv. It will be much more efficient
 * with point-to-point communication. 
 * 
 * \param[in] mpicomm   MPI communicator
 * \param[in] p4est     The forest
 * \param[in] model     Sphere model
 */
static void
p4est_gmt_sphere_communicate_points_old(sc_MPI_Comm mpicomm,
                                    p4est_t *p4est,
                                    p4est_gmt_model_t *model)
{
  /* Sphere model */
  p4est_gmt_model_sphere_t *sdata = (p4est_gmt_model_sphere_t *)model->model_data;

  /* MPI data */
  int                 mpiret;
  int                 num_procs, rank;

  /** communication data. The local process is referred to as p **/
  /* number of points that p will own in the next iteration (responsible)*/
  int                 num_owned_resp = sdata->num_owned_resp; 
  /* number of points that p will own in the next iteration (not responsible)*/
  int                 num_owned = sdata->num_owned;
  /* q -> {points q should own in the next iteration, and be responsible for} */
  sc_array_t        **owned_resp; 
  /* q -> {points q should own in the next iteration, but not be responsible for } */
  sc_array_t        **owned;
  /* the number of bytes we are sending to q */
  int                 send_count;
  /* q -> number of bytes we are receiving from q  */
  int                *recv_counts;
  /* q -> byte offset to receive at */
  int                *displs;
  /** per point data **/
  sc_array_t         *owners; /* who to send the point to */
  int                 responsible; /* who is in charge of propagating the point */
  void               *send_buffer; /* TODO: We don't need this if we're doing P2P*/

  /* Get rank and total process count */
  mpiret = sc_MPI_Comm_size (mpicomm, &num_procs);
  SC_CHECK_MPI (mpiret);
  mpiret = sc_MPI_Comm_rank (mpicomm, &rank);
  SC_CHECK_MPI (mpiret);

  /* Initialise index of outgoing message buffers */
  owned = P4EST_ALLOC(sc_array_t*, num_procs);
  owned_resp = P4EST_ALLOC(sc_array_t*, num_procs);

  /* Initialise outgoing message buffers */
  for (int q = 0; q < num_procs; q++) {
    owned[q] = sc_array_new(model->point_size);
    owned_resp[q] = sc_array_new(model->point_size);
  }

  //printf("Proc %d: Prepare outgoing\n", rank);
  /* Iterate through points to add them to outgoing messages */
  for (int i = 0; i < num_owned_resp; i++) {
    /* Determine which processes should receive this point */
    owners = model->owners_fn(&(sdata->points[i]), p4est);
    /* Determine the receiver responsible for propagating this point in next
        iteration */
    //printf("Proc %d responsible\n", rank);
    responsible = model->resp_fn(&(sdata->points[i]), owners);
    P4EST_ASSERT(0 <= responsible && responsible < num_procs);

    //printf("Proc %d add to receivers\n", rank);
    /* Add this point to the arrays of the processes who receive it */
    for (size_t j = 0; j < owners->elem_count; j++) {
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
                      sc_MPI_SUM, q, mpicomm);
    SC_CHECK_MPI (mpiret);
    mpiret = sc_MPI_Reduce(&(owned_resp[q]->elem_count), &num_owned_resp, 1, sc_MPI_INT, 
                      sc_MPI_SUM, q, mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  printf("Proc %d: Free old points\n", rank);
  /* Free the points we received last iteration */
  P4EST_FREE(sdata->points);
  /* Create array for incoming points */
  sdata->points = P4EST_ALLOC(p4est_gmt_sphere_geodesic_seg_t, num_owned_resp + num_owned);
  /* Update the amount of points known */
  model->M = num_owned_resp + num_owned;
  sdata->num_owned_resp = num_owned_resp;
  sdata->num_owned = num_owned;

  /* Initialise recv_count and displacement arrays */
  recv_counts = P4EST_ALLOC(int, num_procs);
  displs = P4EST_ALLOC(int, num_procs);

  printf("Proc %d: Send points\n", rank);
  /* Transmit points */
  for (int q = 0; q < num_procs; q++) {
    /* Tell process q how many (resp) bytes it will receive from each process */
    send_count = (owned_resp[q]->elem_count) * sizeof(p4est_gmt_sphere_geodesic_seg_t);
    mpiret = sc_MPI_Gather(&send_count,
                              sizeof(int),
                              sc_MPI_BYTE,
                              recv_counts,
                              sizeof(int),
                              sc_MPI_BYTE,
                              q,
                              mpicomm
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
                              (owned_resp[q]->elem_count) * (model->point_size),
                              sc_MPI_BYTE,
                              sdata->points, /* rec buffer */
                              recv_counts,
                              displs,
                              sc_MPI_BYTE,
                              q,
                              mpicomm);
    SC_CHECK_MPI (mpiret);

    /* Tell process q how many (not resp) bytes it will receive from each process */
    send_count = (owned[q]->elem_count) * sizeof(p4est_gmt_sphere_geodesic_seg_t);
    mpiret = sc_MPI_Gather(&send_count,
                              sizeof(int),
                              sc_MPI_BYTE,
                              recv_counts,
                              sizeof(int),
                              sc_MPI_BYTE,
                              q,
                              mpicomm
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
                              (owned[q]->elem_count) * (model->point_size),
                              sc_MPI_BYTE,
                              &(sdata->points[num_owned_resp]), /* rec buffer is offset */
                              recv_counts,
                              displs,
                              sc_MPI_BYTE,
                              q,
                              mpicomm);
    SC_CHECK_MPI (mpiret);
  }

  printf("Proc %d: cleanup sending arrays \n", rank);
  /* Cleanup*/
  /* Destroy outgoing message buffers */
  for (int q = 0; q < num_procs; q++) {
    sc_array_destroy(owned[q]);
    sc_array_destroy(owned_resp[q]);
  }
  /* Free message buffer indexing*/
  P4EST_FREE(owned);
  P4EST_FREE(owned_resp);
  /* Free displacement and receiver counts*/
  P4EST_FREE(displs);
  P4EST_FREE(recv_counts);
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
    model->communicate_points = p4est_gmt_sphere_communicate_points;  
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
