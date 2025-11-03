#pragma once
#include <string>
#include <vector>

#include "exsqs/lattice.hpp"
#include "exsqs/linalg.hpp"

namespace exsqs {

struct Structure {
  Mat3 cell;                       // rows = basis vectors (Angstrom)
  std::vector<Vec3> frac;          // fractional coordinates in [0,1)
  std::vector<int> species;        // index into names
  std::vector<std::string> names;  // one entry per species
  int natoms() const { return static_cast<int>(frac.size()); }
  int nspecies() const { return static_cast<int>(names.size()); }
  std::vector<int> counts() const;
  std::vector<double> concentrations() const;  // achieved x-tilde [A5]
};

// General integer supercell: new basis rows A' = H * A0. det(H) must be > 0.
// Site order is deterministic: lexicographic in the old-lattice translation,
// basis site index innermost.
Structure make_supercell(const Prototype& p, const Mat3i& H,
                         const std::vector<std::string>& names, int species_fill = 0);
Structure make_supercell_diag(const Prototype& p, int n1, int n2, int n3,
                              const std::vector<std::string>& names, int species_fill = 0);

// Replace decoration of an existing geometry (sizes must match).
Structure decorate(const Structure& geom, const std::vector<int>& sigma,
                   const std::vector<std::string>& names);

// Stable reorder of atoms by species index (POSCAR block order). This pins the
// atom-enumeration convention shared between the in-memory cell and its POSCAR
// round-trip, which the T-D1 oracle comparison relies on.
Structure grouped_by_species(const Structure& s);


// Minimal VASP5 POSCAR reader (Direct or Cartesian; positive scale factor).
// Used for `lattice.file` prototypes; species labels are read but a prototype
// treats all sites as one geometric set.
Structure read_poscar(const std::string& path);

// VASP5 POSCAR (Direct coordinates, species grouped, 16 decimals).
void write_poscar(const Structure& s, const std::string& path, const std::string& comment = "");

}  // namespace exsqs
