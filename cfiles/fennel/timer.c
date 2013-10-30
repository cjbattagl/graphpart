#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "timer.h"

/* =================================================== */
/* 
 * Timing functions 
 */
#if !defined(HAVE_TIMER)
#  define TIMER_DESC "gettimeofday"

#define USE_STD_CREATE
#define USE_STD_DESTROY

#include <sys/time.h>

struct stopwatch_t
{
  struct timeval t_start_;
  struct timeval t_stop_;
  int is_running_;
};

static
long double
elapsed (struct timeval start, struct timeval stop)
{
  return (long double)(stop.tv_sec - start.tv_sec)
    + (long double)(stop.tv_usec - start.tv_usec)*1e-6;
}

long double
stopwatch_elapsed (struct stopwatch_t* T)
{
  long double dt = 0;
  if (T) {
    if (T->is_running_) {
      struct timeval stop;
      gettimeofday (&stop, 0);
      dt = elapsed (T->t_start_, stop);
    } else {
      dt = elapsed (T->t_start_, T->t_stop_);
    }
  }
  return dt;
}

void
stopwatch_init (void)
{
  fprintf (stderr, "Timer: %s\n", TIMER_DESC);
  fprintf (stderr, "Timer resolution: ~ 1 us (?)\n");
  fflush (stderr);
}

void
stopwatch_start (struct stopwatch_t* T)
{
  assert (T);
  T->is_running_ = 1;
  gettimeofday (&(T->t_start_), 0);
}

long double
stopwatch_stop (struct stopwatch_t* T)
{
  long double dt = 0;
  if (T) {
    if (T->is_running_) {
      gettimeofday (&(T->t_stop_), 0);
      T->is_running_ = 0;
    }
    dt = stopwatch_elapsed (T);
  }
  return dt;
}

#  define HAVE_TIMER 1
#endif

#if defined(USE_STD_CREATE)
struct stopwatch_t *
stopwatch_create (void)
{
  struct stopwatch_t* new_timer =
    (struct stopwatch_t *)malloc (sizeof (struct stopwatch_t));
  if (new_timer)
    memset (new_timer, 0, sizeof (struct stopwatch_t));
  return new_timer;
}
#endif

#if defined(USE_STD_DESTROY)
void
stopwatch_destroy (struct stopwatch_t* T)
{
  if (T) {
    stopwatch_stop (T);
    free (T);
  }
}
#endif
/* =================================================== */

