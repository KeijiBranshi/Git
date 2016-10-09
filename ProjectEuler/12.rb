require './lib/math/TriangleNumbers'
include TriangleNumbers

divisorTarget = 500
nDivisors = 0
count = 0

while nDivisors < divisorTarget
  triNum = TriangleNumbers::getFactorized(count)
  key = triNum.keys
  nDivisors = triNum[key[0]].size
  count += 1
end

p triNum
