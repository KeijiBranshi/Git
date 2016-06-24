class Fixnum
  # digits method obtiained from http://stackoverflow.com/questions/13091558/how-do-i-iterate-through-the-digits-of-an-integer
  def digits(base: 10)
    quotient, remainder = divmod(base)
    quotient == 0 ? [remainder] : [*quotient.digits(base: base), remainder]
  end

  # DONT USE BROKEN: num_digits obtained from http://stackoverflow.com/questions/14005524/how-do-i-determine-the-length-of-a-fixnum-in-ruby
  def num_digits_10
    Math.log10(self).to_i + 1
  end
end
