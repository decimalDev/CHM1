// lab1.cpp: определяет точку входа для приложения.
//

#include "lab1.h"
#include "Lagrange.cpp"
#include "Fourier.cpp"
#include "Task5.cpp"
#include "Task6.cpp"

using namespace std;

int main()
{
	//cout << "Hello CMake." << endl;
	//Lagrange_task();
	//NM_task4();
	//task5();
	task6();
	std::cout << "write any letter" << std::endl;
	char c;
	std::cin >> c;
	return 0;
}
