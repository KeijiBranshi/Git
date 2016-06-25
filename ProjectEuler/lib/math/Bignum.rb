class Bignum
  # Only aquires digits in base 10 right now
  def digits(base: 10)
    arr = Array.new()
    self.to_s.each_char do |c|
      arr << c.to_i
    end
    return arr
  end

  # DONT USE BROKEN: num_digits obtained from http://stackoverflow.com/questions/14005524/how-do-i-determine-the-length-of-a-fixnum-in-ruby
  def num_digits_10
    Math.log10(self).to_i + 1
  end
end
