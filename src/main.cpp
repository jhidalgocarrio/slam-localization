#include <iostream>
#include <rover_localization/Sckf.hpp>
#include <rover_localization/Lesma.hpp>

int main(int argc, char** argv)
{
	localization::Sckf mysckf;
	mysckf.welcome();

	return 0;
}
