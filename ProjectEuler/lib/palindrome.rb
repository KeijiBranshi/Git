module Palindrome
	def is_palindrome?(input)
		string = ""
		if !input.instance_of? String
			string = input.to_s
=begin
			if string[string.length-1] == "0"
				string.insert 0, "0"
			end
=end
		else
			string = input
		end

		if string == string.reverse
			return true
		else
			return false
		end
	end

	#make a palindrome out of the input to the right
	def forward(input)
		string
		if !input.instance_of? String
			string = input.to_s
		else
			string = input
		end

		string << string.reverse

		if input.instance_of? Integer
			return string.to_i
		else
			return string
		end
	end

	#same as above but to the left
	def backward(input)
		string
		if !input.instance_of? String
			string = input.to_s
		else
			string = input
		end

		string.insert 0, string.reverse

		if input.instance_of? Integer
			return string.to_i
		else
			return string
		end
	end

end