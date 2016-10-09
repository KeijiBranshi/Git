class Fixnum
  # digits method obtiained from http://stackoverflow.com/questions/13091558/how-do-i-iterate-through-the-digits-of-an-integer
  def digits(base: 10)
    quotient, remainder = divmod(base)
    quotient == 0 ? [remainder] : [*quotient.digits(base: base), remainder]
  end

  # returns an array with factors of a given number
  def factors() (1..self).select { |n| (self % n).zero? } end
end
