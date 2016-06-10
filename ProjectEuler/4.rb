#Largest palindromic number that is the multiple of two 3-digit numbers
require '../Ruby/lib/palindrome.rb'
include Palindrome

999.downto(100).each do |i|
	@answer = 0
	@found = false

	i.downto(100).each do |j|
		print (i*j)
		puts " = #{i} * #{j}"
		next unless Palindrome::is_palindrome?(i * j)

		@found = true
		@answer = i*j
		puts "ANSWER: #{@answer} = #{i} * #{j}"
		break
	end

	if @found 
		break
	end
end

2895 before

24 per room water sewer trash 

elect internet, gas, A,b,c 100 per month

pay them 24 per month for electricity and plumbing
we do internet and gas
