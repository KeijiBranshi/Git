require "./lib/math/factorial.rb"

#Easy way (build smallest number via pen and paper)
answer = (2**4) * (3**2) * 5 * 7 * 11 * 13 * 17 * 19
puts "Quicker Answer: #{answer}"

#Long way (to help me get used to Ruby)
=begin
last = factorial(20)
@found = false
@array = [6,7,8,9,11,13,16,17,19,20]
answer = 0

(20..last).step(2) do |num|

  @array.each do |i|
    temp = (num % i)
    puts "#{num} % #{i} = #{temp}"

    if temp != 0
      @found = false
      break
    else
      @found = true
    end
  end
  
  if @found
    @answer = num
    break
  else
    puts "Not this one\n\n"
  end

end

puts "Answer: #{@answer}"
=end

