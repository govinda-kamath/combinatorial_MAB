#include <iostream>
#include "Test.h"
using namespace std;

Test::Test(string name):name(name){}

void Test::display(){
	cout << "A student with name " << this->name << endl;
}

