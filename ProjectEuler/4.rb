#Largest palindromic number that is the multiple of two 3-digit numbers
require './lib/palindrome.rb'
include Palindrome

$answer = 0
$x = 0
$y = 0

999.downto(100).each do |i|
	i.downto(100).each do |j|
		product = i*j
=begin
		if product.to_s[0] == product.to_s[product.to_s.length-1]
			print (product)
			puts " = #{i} * #{j}"
		end
=end
		next unless Palindrome::is_palindrome?(product)

		if product > @answer
			@answer = product
			@x = i
			@y = j
		end
	end
end

puts "ANSWER: #{@answer} = #{@x} * #{@y}"
