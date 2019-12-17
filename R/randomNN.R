#####################################
# HELPER FUNCTIONS:
# randomNNs
#####################################

# randomNN -------------------------------------------------------------------------
# Find the random nearest neighbors.
#
# Gives the indices of nearest neighbors of points that are not unique in the data set.
# For each point we know that to which cluster of points it belongs (repeated points),
# we chose one of those indices that is not equal to the same index of our given point
# at random as its nearest neighbor.
#
# @param ids: a vector of ids of groups that each index is a member of.
#
# @return a vector of indices.
.randomNN <- function(ids) {
  m <- length(ids)

  x <- sample(x = (m - 1), m, replace = TRUE)
  x <- x + (x >= (1:m))

  return(ids[x])
}
