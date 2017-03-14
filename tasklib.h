#ifndef TASKLIB_H
#define TASKLIB_H

// -----------------------------------------------------------------------------
// PHISICAL CONSTANT
// -----------------------------------------------------------------------------
#define G			9.81		// gravity
#define TSCALE		1			// time scale factor
// -----------------------------------------------------------------------------
// TASK CONSTANT
// -----------------------------------------------------------------------------
#define N_TASK		10			// maximum number of task
// -----------------------------------------------------------------------------
// STRUCTURE
// -----------------------------------------------------------------------------
// Structure containing the task parameters
struct task_par {
	int arg;
	int period;
	int deadline;
	int priority;
	int dmiss;
	struct timespec at;
	struct timespec dl;
};
// -----------------------------------------------------------------------------
// EXTERN VARIABLES
// -----------------------------------------------------------------------------
extern int 		npend;
extern struct 	sched_param *mypar;
extern struct 	task_par 	tp[MAX_DP + 1];
extern int 		end;
// -----------------------------------------------------------------------------
// TASK FUNCTIONS
// -----------------------------------------------------------------------------
void 	*pend_task(void *arg);
void	*dline_wintask(void *arg);
// -----------------------------------------------------------------------------
// THREAD MANAGEMENT FUNCIONS
// -----------------------------------------------------------------------------
void 	task_init();
void 	task_create(int i, void *task_fun, int period, int dline, int pri);
void 	wait_for_task_end(int i);

#endif