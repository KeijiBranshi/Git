sum = 100 * (101) / 2
sumofsq = 0
sqofsum = 0
for i in 1..100
  sumofsq += i**2
end

sqofsum = sum**2

answer = sqofsum - sumofsq
puts "Answer: #{answer}"
