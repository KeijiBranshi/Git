fibo1 = 0
fibo2 = 1
sum = 0

until fibo2 >= 4000000
    temp = fibo2
    fibo2 += fibo1
    fibo1 = temp

    puts "(#{fibo1}, #{fibo2})"

    sum += fibo1 if fibo1%2 == 0
end

puts sum
