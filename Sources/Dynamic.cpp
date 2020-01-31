#define _USE_MATH_DEFINES

#include "Dynamic.h"

#include <cmath>
#include <algorithm>

using namespace std;

Dynamic::~Dynamic () {
	clear ();
}
