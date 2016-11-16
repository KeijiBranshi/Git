# Rotates coordinates in a single frame xyz file
require 'matrix'
require_relative 'Ruby-lib/xyz_file.rb'
require_relative 'Ruby-lib/matrix-ops.rb'
include Math
include XYZ

def allZero(array)
  array.each do |a|
    if a != 0
      return false
    end
  end
  true
end
# Get rotation and translation data matrix
thetaX = ARGV[1].to_f
thetaY = ARGV[2].to_f
thetaZ = ARGV[3].to_f

dx = ARGV[4].to_f
dy = ARGV[5].to_f
dz = ARGV[6].to_f

# Get coordinate data from file. File returns array with [header line, comment line, matrix of coord data]
header, comment, mat = getData(ARGV[0])

#Apply Transormations
mat = translate(rotate(mat.transpose, thetaX, thetaY, thetaZ), dx, dy, dz).transpose.round(2)

#Determine what to rename file
transformation = ""
if !allZero([dx,dy,dz])
  transformation += "_TRANS"
end
if !allZero([thetaX,thetaY,thetaZ])
  transformation += "_ROT"
end

# Write back to file
filename = ARGV[0].dup
filename.insert(filename.index('.'), transformation)
lineEnding = nil
putData([header, comment, mat], filename, lineEnding)
