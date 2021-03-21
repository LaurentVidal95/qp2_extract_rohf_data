=================
extract_rohf_data
=================

A small plugin that generates four files that are to be read by the julia package
ROHFToolkit [ see (insert github here)] to start ROHF computations.

The files are:
- `overlap_matrix.dat` that contains the overlap matrix.
- `H_core.dat`  that contains the core hamiltonian matrix.
- `nums_orbitals.dat` that containes the number of basis functions, doubly occupied and
  singly occupied orbitals.
- `Four_index_tensor.dat` that containes the four-index tensor in a condensed writing
  by taking symmetries into account. 
