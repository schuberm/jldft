# Overlap operator (acting on 3d data sets)
#
# Usage: out=O(in)
#
# in: input 3d data set
# out: output 3d data set
#
# Uses GLOBAL variable(s) ---
# gbl_R: Lattice vectors

function O(in)
  global gbl_R; ## Must declare all globals with such statements to access them
  out= det(gbl_R)*in
  return out
end
