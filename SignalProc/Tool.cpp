#include "Tool.h"


Tool::Tool()
{
}


Tool::~Tool()
{
}

long Tool::multiTo(int n)
{
	if (n == 0) return 1;

	long p = 1;
	for (int i = 1; i <= n; i++){
		p *= i;
	}
	return p;
}
