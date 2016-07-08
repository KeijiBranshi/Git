=begin
To help setup the bounds of my loops that will find a, b, and c, I need to figure out reasonable bounds to iterate through.

Given that:
1. a < b < c
2. a^2 + b^2 = c^2
3. a+b+c=1000

By (1) and (2), we know that a,b,c > 0 (otherwise two within a,b,c could be equal). Therefore, the minimum they any of them could be is 1. If a, b, xor c is 1, then then 2 is the smallest value either of the others could be.
-> By (3), the biggest a,b,c could possibly be is 997.


=end

@a = 0
@b = 0
@c = 0
found = false

(1..996).each do |a|
  (a...996).each do |b|
    c = Math.sqrt(a**2 + b**2)

    if (c < a || c < b)
      next
    elsif (a+b+c == 1000)
      @a = a
      @b = b
      @c = c
      found = true
      break
    end
  end
  if found then break end
end

puts ("a=#{@a} b=#{@b} c=#{@c}\n")
solution = @a * @b * @c
p solution
=begin
print ("#{@a}^2 = ")
p @a**2
print ("#{@b}^2 = ")
p @b**2
print ("#{@c}^2 = ")
p @c**2
=end
