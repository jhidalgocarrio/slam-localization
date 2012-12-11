#include <iostream>
#include <rover_localization/sckf.hpp>
#include <rover_localization/lesma.hpp>

int main(int argc, char** argv)
{
	localization::sckf mysckf;
	mysckf.welcome();

	return 0;
}
