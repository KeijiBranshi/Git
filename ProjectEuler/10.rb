require 'prime'

max = 2000000
sum = 0

(2...max).each do |option|
  if (option < max && Prime.prime?(option))
    sum += option
  else next
  end
end

p sum
