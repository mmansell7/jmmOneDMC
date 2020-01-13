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


#ifndef COMPUTE_H
#define COMPUTE_H

#include <stdbool.h>
#include <vector>
#include <string>

class Compute {
  public:
   static int instance_total;     // # of Compute classes ever instantiated

   char *id,*style;
   int igroup,groupbit;

   double scalar;            // computed global scalar
   double *vector;           // computed global vector
   double **array;           // computed global array
   double *vector_atom;      // computed per-atom vector
   double **array_atom;      // computed per-atom array
   double *vector_local;     // computed local vector
   double **array_local;     // computed local array

   int scalar_flag;          // 0/1 if compute_scalar() function exists
   int vector_flag;          // 0/1 if compute_vector() function exists
   int array_flag;           // 0/1 if compute_array() function exists
   int size_vector;          // length of global vector
   int size_array_rows;      // rows in global array
   int size_array_cols;      // columns in global array
   int size_vector_variable;      // 1 if vec length is unknown in advance
   int size_array_rows_variable;  // 1 if array rows is unknown in advance

   int peratom_flag;         // 0/1 if compute_peratom() function exists
   int size_peratom_cols;    // 0 = vector, N = columns in peratom array

   int local_flag;           // 0/1 if compute_local() function exists
   int size_local_rows;      // rows in local vector or array
   int size_local_cols;      // 0 = vector, N = columns in local array

   int extscalar;            // 0/1 if global scalar is intensive/extensive
   int extvector;            // 0/1/-1 if global vector is all int/ext/extlist
   int *extlist;             // list of 0/1 int/ext for each vec component
   int extarray;             // 0/1 if global array is all intensive/extensive

   int tempflag;       // 1 if Compute can be used as temperature
                       // must have both compute_scalar, compute_vector
   int pressflag;      // 1 if Compute can be used as pressure (uses virial)
                       //  // must have both compute_scalar, compute_vector
   int pressatomflag;  // 1 if Compute calculates per-atom virial
   int peflag;         // 1 if Compute calculates PE (uses Force energies)
   int peatomflag;     // 1 if Compute calculates per-atom PE
   int create_attribute;    // 1 if compute stores attributes that need
                            // setting when a new atom is created
 
   int tempbias;       // 0/1 if Compute temp includes self/extra bias
 
   int timeflag;       // 1 if Compute stores list of timesteps it's called on
   int ntime;          // # of entries in time list
   int maxtime;        // max # of entries time list can hold
   unsigned long int *tlist;      // list of timesteps the Compute is called on

   std::string str;
   std::vector<std::string> strs;

  int invoked_flag;       // non-zero if invoked or accessed this step, 0 if not
  unsigned long int invoked_scalar;  // last timestep on which compute_scalar() was invoked
  unsigned long int invoked_vector;  // ditto for compute_vector()
  unsigned long int invoked_array;   // ditto for compute_array()
  unsigned long int invoked_peratom; // ditto for compute_peratom()
  unsigned long int invoked_local;   // ditto for compute_local()

  double dof;         // degrees-of-freedom for temperature

  int comm_forward;         // size of forward communication (0 if none)
  int comm_reverse;         // size of reverse communication (0 if none)
  int dynamic_group_allow;  // 1 if can be used with dynamic group, else 0

  int copymode;

  Compute(char[30][30]);
  Compute(char**);
  virtual ~Compute();
  void modify_params(int, char **);
  
  
  // virtual void init() = 0;
  virtual void setup() {}
  virtual double compute_scalar() {return 0.0;}
  virtual void compute_vector() {}
  virtual void compute_array() {}
  virtual void compute_peratom() {}
  virtual void compute_local() {}

 protected:
  int instance_me;             // which Compute class instantiation I am

};

#endif

