#Triangle number calculator from Project Euler problem 12
module TriangleNumbers
  triNumArray = []
  
  #Returns the Nth triangle number
  def get(n=0)
    sum = 0
    (0...n).each do |i|
      sum += i
    end
    return sum
  end

  def getFactorized(n=0)
    require './lib/math/Fixnum'
    sum = self::get(n)
    return {sum => sum.factors}
  end

  #Returns an array of all triangle numbers up to N
  def gather(n=0)
    array = []
    sum = 0
    (0...n).each do |i|
      sum += i
      array << sum
    end
    return array
  end

  def gatherFactorized(n=0)
    require './lib/math/Fixnum'
    hash = {}
    sum = 0
    for i in (0...n)
      sum += i
      hash[sum] = sum.factors
    end
    return hash
  end
end
