#include <signal.h>
#include "ma.h"

typedef void (*fptr)(int);

char CtrlBreakHit = 0;
fptr saveInt1B;

void Int1BProg (int dummy)
{
	signal (SIGINT, (fptr)Int1BProg);
	CtrlBreakHit = !0;
}

void InitBreakControl (void)
{
	saveInt1B = signal (SIGINT, Int1BProg);
}

void DoneBreakControl (void)
{
	if (saveInt1B!=SIG_ERR) signal (SIGINT, saveInt1B);
}
