/* Sets the two private fields that Mandyoc fills in createSwarm_2d
 * (dmcell->ops->locatepoints and ->getneighbors). The callbacks are in Julia. */
#include <petsc/private/dmimpl.h>

PetscErrorCode mandyoc_set_cell_dm_ops(DM cell_dm,
  PetscErrorCode (*locate_points)(DM, Vec, DMPointLocationType, PetscSF),
  PetscErrorCode (*get_neighbors)(DM, PetscInt *, const PetscMPIInt **))
{
  cell_dm->ops->locatepoints = locate_points;
  cell_dm->ops->getneighbors = get_neighbors;
  return PETSC_SUCCESS;
}
