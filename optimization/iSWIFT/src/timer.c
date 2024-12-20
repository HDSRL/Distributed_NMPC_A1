#include "timer.h"

#if (defined WIN32 || _WIN64)

void tic(timer* t)
{
	QueryPerformanceFrequency(&t->freq);
	QueryPerformanceCounter(&t->tic);
}

realqp toc(timer* t)
{
	QueryPerformanceCounter(&t->toc);
	return ((t->toc.QuadPart - t->tic.QuadPart) / (realqp)t->freq.QuadPart);
}


#elif (defined __APPLE__)

void tic(timer* t)
{
	/* read current clock cycles */
	t->tic = mach_absolute_time();
}

realqp toc(timer* t)
{

	uint64_t duration; /* elapsed time in clock cycles*/

	t->toc = mach_absolute_time();
	duration = t->toc - t->tic;

	/*conversion from clock cycles to nanoseconds*/
	mach_timebase_info(&(t->tinfo));
	duration *= t->tinfo.numer;
	duration /= t->tinfo.denom;

	return (realqp)duration / 1000000000;
}



#else

/* read current time */
void tic(timer* t)
{
	clock_gettime(CLOCK_MONOTONIC, &t->tic);
}


/* return time passed since last call to tic on this timer */
double toc(timer* t)
{
	struct timespec temp;

	clock_gettime(CLOCK_MONOTONIC, &t->toc);

	if ((t->toc.tv_nsec - t->tic.tv_nsec)<0) {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec - 1;
		temp.tv_nsec = 1000000000 + t->toc.tv_nsec - t->tic.tv_nsec;
	}
	else {
		temp.tv_sec = t->toc.tv_sec - t->tic.tv_sec;
		temp.tv_nsec = t->toc.tv_nsec - t->tic.tv_nsec;
	}
	return (realqp)temp.tv_sec + (realqp)temp.tv_nsec / 1000000000;
}

#endif
