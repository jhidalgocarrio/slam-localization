#include <iostream>
#include <asguard_localization/sckf.hpp>
#include <asguard_localization/lesma.hpp>

int main(int argc, char** argv)
{
	localization::sckf mysckf;
	mysckf.welcome();

	return 0;
}
