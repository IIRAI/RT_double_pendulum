#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <allegro.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <time.h>
#include "graphlib.h"
#include "tasklib.h"

int main()
{
int 	i;

	read_data();
	display_init();
	task_init();

	listen_keyboard();

	for (i = 0; i < N_TASK; i++) wait_for_task_end(i);

	allegro_exit();
	return 0;
}
