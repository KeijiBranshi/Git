require 'matrix'

def getRotation(thetaX, thetaY, thetaZ)
  rotationX = Matrix[[1,0,0],[0,cos(thetaX),-sin(thetaX)],[0,sin(thetaX),cos(thetaX)]]
  rotationY = Matrix[[cos(thetaY),0,sin(thetaY)],[0,1,0],[-sin(thetaY),0,cos(thetaY)]]
  rotationZ = Matrix[[cos(thetaZ),-sin(thetaZ),0],[sin(thetaZ),cos(thetaZ),0],[0,0,1]]
  rotation = rotationX * rotationY * rotationZ
end

def rotate(matrix, thetaX, thetaY, thetaZ)
  getRotation(thetaX, thetaY, thetaZ) * matrix
end

def getTranslation(dx, dy, dz)
  translation = Matrix[[1,0,0,dx],[0,1,0,dy],[0,0,1,dz],[0,0,0,1]]
end

def translate(matrix, dx, dy, dz)
  tempMatrix = Matrix.rows(matrix.to_a << Array.new(matrix.column_count(), 1))  #Add row of ones for translation
  newMatrixArray = (getTranslation(dx, dy, dz) * tempMatrix).to_a
  newMatrixArray.pop    #pop the last row of ones

  Matrix.rows(newMatrixArray)   #return the final matrix
end
