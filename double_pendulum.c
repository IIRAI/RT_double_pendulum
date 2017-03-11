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
char 	scan = 0;	

	read_data();
	display_init();
	task_init();

	do {
		scan = get_scancode();
		execute_scan(scan);
	} while (scan != KEY_ESC);

	end = 1;
	for (i = 0; i <= npend; i++) wait_for_task_end(i);

	allegro_exit();	
	return 0;
}
