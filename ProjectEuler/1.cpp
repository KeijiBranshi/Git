#include <iostream>

int main()
{
	int sum = 0;
	for (int i; i<1000; ++i)
	{
		if ( !(i%3) || !(i%5) )
			sum += i;
		else
			continue;
	}

	std::cout << "Answer is " << sum << std::endl;
}