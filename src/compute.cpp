/* --------------------------------------------------------------------------
 * Much of this file is copied (and modified) from LAMMPS - Large-scale
 * Atomic/Molecular Massively Parallel Simulator
 * http://lammps.sandia.gov, Sandia National Laboratories
 * Steve Plimpton, sjplimp@sandia.gov
 *
 * Copyright (2003) Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software.  This software is distributed under
 * the GNU General Public License.
 *
 * See the README file in the top-level LAMMPS directory.
 * ------------------------------------------------------------------------- */

#include "compute.h"
#include <cstring>
#include <cctype>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "stdbool.h"
#include <vector>
#include <string>

// allocate space for static class instance variable and initialize it

int Compute::instance_total = 0;

/* ---------------------------------------------------------------------- */

Compute::Compute(char arg[30][30]) :
  id(NULL), style(NULL),
  vector(NULL), array(NULL), vector_atom(NULL),
  array_atom(NULL), vector_local(NULL), array_local(NULL), extlist(NULL),
  tlist(NULL)
{
  int ii,jj;
  
  instance_me = instance_total++;

  std::vector<std::string> strs;
  std::string str = "";
  ii = 0;
  while ( arg[ii] && strncmp(arg[ii],"",30) ) {
    strs.push_back(std::string(arg[ii]));
    str = str + " " + strs[ii];
    ii++;
  }

  // compute ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = atoi(arg[1]);

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;
  size_vector_variable = size_array_rows_variable = 0;

  tempflag = pressflag = peflag = 0;
  pressatomflag = peatomflag = 0;
  create_attribute = 0;
  tempbias = 0;

  timeflag = 0;

  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_peratom = invoked_local = -1;
  invoked_flag = 0;

  // set modify defaults

  // setup list of timesteps

  ntime = maxtime = 0;

  // data masks

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Compute::Compute(char **arg) :
  id(NULL), style(NULL),
  vector(NULL), array(NULL), vector_atom(NULL),
  array_atom(NULL), vector_local(NULL), array_local(NULL), extlist(NULL),
  tlist(NULL)
{
  int ii,jj;
  
  instance_me = instance_total++;

  ii = 0;
  while ( arg[ii] && strncmp(arg[ii],"",30) ) {
    strs.push_back(std::string(arg[ii]));
    str = str + strs[ii] + " ";
    ii++;
  }
  
  // compute ID, group, and style
  // ID must be all alphanumeric chars or underscores

  int n = strlen(arg[0]) + 1;
  id = new char[n];
  strcpy(id,arg[0]);

  igroup = atoi(arg[1]);

  n = strlen(arg[2]) + 1;
  style = new char[n];
  strcpy(style,arg[2]);

  // set child class defaults

  scalar_flag = vector_flag = array_flag = 0;
  peratom_flag = local_flag = 0;
  size_vector_variable = size_array_rows_variable = 0;

  tempflag = pressflag = peflag = 0;
  pressatomflag = peatomflag = 0;
  create_attribute = 0;
  tempbias = 0;

  timeflag = 0;

  invoked_scalar = invoked_vector = invoked_array = -1;
  invoked_peratom = invoked_local = -1;
  invoked_flag = 0;

  // set modify defaults

  // setup list of timesteps

  ntime = maxtime = 0;

  // data masks

  copymode = 0;
}

/* ---------------------------------------------------------------------- */

Compute::~Compute()
{
  if (copymode) return;

  delete [] id;
  delete [] style;
  delete[] tlist;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

