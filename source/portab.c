#include "portab.h"
#ifdef _GCC_SUN_
#include <stdlib.h>
void *memmove(void *dest, const void *src, size_t n)
{
	bcopy (src, dest, n);
	return dest;
}
#endif

