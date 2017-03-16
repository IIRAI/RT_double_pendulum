/*
COMMENT HERE
*/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <pthread.h>
#include <sched.h>
#include <time.h>
#include "graphlib.h"
#include "tasklib.h"

struct 			sched_param *mypar;
struct 			task_par 	tp[N_TASK];
pthread_attr_t 	myatt[N_TASK];
pthread_t 		pend[N_TASK];
int 			end = 0;						// flag to terminate the task
float 			h;								// increment interval

// -----------------------------------------------------------------------------
// MATRIXES TO STORE THE VALUES FOR THE PENDULUM'S MOTION
// -----------------------------------------------------------------------------
double 			k1[MAX_DP][4];
double 			k2[MAX_DP][4];
double 			k3[MAX_DP][4];
double 			k4[MAX_DP][4];
double 			m[MAX_DP][4];			// memorizza l'incremento del parametro
// -----------------------------------------------------------------------------
// PENDULUM'S MOTION FUNCTIONS
// -----------------------------------------------------------------------------
static void 	RK4(int i);
static void 	M(int i);
static void 	K1(int i);
static void 	K2(int i);
static void 	K3(int i);
static void 	K4(int i);
static double 	f1(double tht1_dot);
static double 	f2(double tht2_dot);
static double 	f3(double tht1, double tht2, double tht1_dot, double tht2_dot, int i);
static double 	f4(double tht1, double tht2, double tht1_dot, double tht2_dot, int i);
// -----------------------------------------------------------------------------
// THREAD MANAGEMENT FUNCIONS
// -----------------------------------------------------------------------------
static void 	set_priority(int i, struct sched_param *mypar);
static int 		task_argument(void *arg);
static int 		task_period(int i);
static void 	set_period(int index);
static void 	wait_for_period(int i);
static int 		deadline_miss(int i);
static void 	time_add_ms(struct timespec *t, int ms);
static void 	time_copy(struct timespec *td, struct timespec ts);
static int 		time_cmp(struct timespec t1, struct timespec t2);

// -----------------------------------------------------------------------------
// PEND_TASK: thread function of the double pendulum. It evaluates the 
// coordinate values and draw them on the backgrund bitmap
// -----------------------------------------------------------------------------
void *pend_task(void *arg)
{
int 	i;
int 	dmiss;

	i = task_argument(arg);
	set_period(i);

	h = TSCALE * (float)task_period(i) / 1000.0;

	while (!end) {

			RK4(i);				// calculates the new value of the pendulum	
			draw_pend(i);

			deadline_miss(i);
			wait_for_period(i);
	}
}

// -----------------------------------------------------------------------------
// DLINE_WINTASK: thread function of the deadline status window. 
// -----------------------------------------------------------------------------
void *dline_wintask(void *arg)
{
int 	a;

	a = task_argument(arg);
	set_period(a);

	while(!end) {

		draw_dline();

		deadline_miss(a);
		wait_for_period(a);
	}	
}

// -----------------------------------------------------------------------------
// RK4: it evaluates the new values of the pendulum
// -----------------------------------------------------------------------------
void RK4(int i)
{
	// aggiorna i valori di incremento
	M(i);

	pnd[i].tht1 		+= m[i][0];
	pnd[i].tht2 		+= m[i][1];
	pnd[i].tht1_dot 	+= m[i][2];
	pnd[i].tht2_dot 	+= m[i][3];
	// printf("theta_dot_1: %1.15f\n", pnd[i].tht1_dot);
	// printf("theta_dot_2: %1.15f\n", pnd[i].tht2_dot);
}

// -----------------------------------------------------------------------------
// M: it evaluates the increments of the new values of the pendulum
// -----------------------------------------------------------------------------
void M(int i)
{
int 	a;
	// aggiorna i valori di K
	K1(i);
	K2(i);
	K3(i);
	K4(i);

	for (a = 0; a < 4; a++)
		m[i][a] = (k1[i][a] + (2 * k2[i][a]) + (2 * k3[i][a]) + k4[i][a])/6;
}

void K1(int i)
{
double 	tht1_k;
double 	tht2_k;
double 	tht1_dot_k;
double 	tht2_dot_k;

	tht1_k 		= pnd[i].tht1;
	tht2_k 		= pnd[i].tht2;
	tht1_dot_k 	= pnd[i].tht1_dot;
	tht2_dot_k 	= pnd[i].tht2_dot;

	k1[i][0] = h * f1(tht1_dot_k);
	k1[i][1] = h * f2(tht2_dot_k);
	k1[i][2] = h * f3(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
	k1[i][3] = h * f4(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
}

void K2(int i)
{
double 	tht1_k;
double 	tht2_k;
double 	tht1_dot_k;
double 	tht2_dot_k;

	tht1_k 		= pnd[i].tht1 		+ (k1[i][0] / 2.0);
	tht2_k 		= pnd[i].tht2 		+ (k1[i][1] / 2.0);
	tht1_dot_k 	= pnd[i].tht1_dot 	+ (k1[i][2] / 2.0);
	tht2_dot_k 	= pnd[i].tht2_dot 	+ (k1[i][3] / 2.0);

	k2[i][0] = h * f1(tht1_dot_k);
	k2[i][1] = h * f2(tht2_dot_k);
	k2[i][2] = h * f3(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
	k2[i][3] = h * f4(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
}

void K3(int i)
{
double 	tht1_k;
double 	tht2_k;
double 	tht1_dot_k;
double 	tht2_dot_k;

 	tht1_k 		= pnd[i].tht1 		+ (k2[i][0] / 2.0);
 	tht2_k 		= pnd[i].tht2 		+ (k2[i][1] / 2.0);
 	tht1_dot_k 	= pnd[i].tht1_dot 	+ (k2[i][2] / 2.0);
 	tht2_dot_k 	= pnd[i].tht2_dot 	+ (k2[i][3] / 2.0);

	k3[i][0] = h * f1(tht1_dot_k);
	k3[i][1] = h * f2(tht2_dot_k);
	k3[i][2] = h * f3(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
	k3[i][3] = h * f4(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
}

void K4(int i)
{
double 	tht1_k;
double 	tht2_k;
double 	tht1_dot_k;
double 	tht2_dot_k;

 	tht1_k 		= pnd[i].tht1 		+ k3[i][0];
 	tht2_k 		= pnd[i].tht2 		+ k3[i][1];
 	tht1_dot_k 	= pnd[i].tht1_dot 	+ k3[i][2];
 	tht2_dot_k 	= pnd[i].tht2_dot 	+ k3[i][3];

	k4[i][0] = h * f1(tht1_dot_k);
	k4[i][1] = h * f2(tht2_dot_k);
	k4[i][2] = h * f3(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
	k4[i][3] = h * f4(tht1_k, tht2_k, tht1_dot_k, tht2_dot_k, i);
}

double f1(double tht1_dot)
{
	return tht1_dot;
}

double f2(double tht2_dot)
{
	return tht2_dot;  
}

double f3(double tht1, double tht2, double tht1_dot, double tht2_dot, int i)
{
double 	num;				// numerator
double 	den;				// denominator
double 	m12; 				// (m1 / m2) ratio
double 	s_tht1;				// sin(tht1)
double 	s_tht2;				// sin(tht2)
double 	s_tht12;			// sin(tht1 - tht2)
double 	c_tht12;			// cos(tht1 - tht2)
double 	N1;					// first element of the numerator
double 	N2;					// second 	" " " "
double 	N3;					// third	" " " "	 
double 	l1;
double 	l2;

	l1 = pnd[i].l1;
	l2 = pnd[i].l2;

	m12 = (pnd[i].m1 / pnd[i].m2);

	s_tht1 = sin(tht1);
	s_tht2 = sin(tht2);
		// printf("sin tht1: %1.15lf\n", s_tht1);
		// printf("sin tht2: %1.15lf\n", s_tht2);

	s_tht12 = sin(tht1 - tht2);
	c_tht12 = cos(tht1 - tht2);
		// printf(ssin tht12: %1.15lf\n", s_tht12);
		// printf(scos tht12:%1.15lf\n", c_tht12);

	N1 = (G / l1) * ((s_tht2 * c_tht12) - (s_tht1 * (m12 + 1)));
	N2 = pow(tht1_dot, 2) * c_tht12 * s_tht12;
	N3 = pow(tht2_dot, 2) * (l2/l1) * s_tht12;

	num = N1 - N2 - N3;
	den = m12 + 1 - (pow (c_tht12, 2));

		// printf("robo: %1.15lf\n",((s_tht2 * c_tht12) - (s_tht1 * (m12 + 1)) ));
		// printf("robo: %1.15lf\n",(G / l1));
		// printf("N1: %1.15lf\n", N1);
		// printf("N2: %1.15lf\n", N2);
		// printf("N3: %1.15lf\n", N3);
		// printf("f3: %1.15lf\n", num/den);

	return (num / den);
}

double f4(double tht1, double tht2, double tht1_dot, double tht2_dot, int i)
{
double 	num;
double 	tht1_dotdot;
double 	s_tht2;				// sin(tht2)
double 	s_tht12;			// sin(tht1 - tht2)
double 	c_tht12;			// cos(tht1 - tht2)
double 	N1;					// 1st element of the numerator
double 	N2;					// 2nd element of the numerator
double 	N3;					// 3rd element of the numerator
double 	l1;
double 	l2;

	l1 = pnd[i].l1;
	l2 = pnd[i].l2;

	tht1_dotdot = f3(tht1, tht2, tht1_dot, tht2_dot, i);

	s_tht2 = sin(tht2);

	s_tht12 = sin(tht1 - tht2);
	c_tht12 = cos(tht1 - tht2);

	N1 = (l1 /l2) * pow(tht1_dot, 2) * s_tht12;
	N2 = (l1 /l2) * tht1_dotdot * c_tht12;
	N3 = (G / l2) * s_tht2; 

	num = N1 - N2 - N3;

	return num;
}

// -----------------------------------------------------------------------------
// THREAD MANAGEMENT FUNCTIONS
// -----------------------------------------------------------------------------
void task_init(void)
{
int 	i;

	for (i = 0; i <= npend; i++) {
		pthread_attr_init(&myatt[i]);
		pthread_attr_setinheritsched(&myatt[i], PTHREAD_EXPLICIT_SCHED);
		pthread_attr_setschedpolicy(&myatt[i], SCHED_FIFO);
	}
}

void task_create(int i, void *task_fun, int period, int dline, int pri)
{
struct 	sched_param mypar;

	tp[i].arg = i;
	tp[i].period = period;
	tp[i].deadline = dline;
	tp[i].priority = pri;

	pend[i] = (pthread_t) i;

	set_priority(i, &mypar);
	pthread_create(&pend[i], &myatt[i], task_fun, &tp[i]);
}

void set_priority(int i, struct sched_param *mypar)
{
	mypar->sched_priority = tp[i].priority;
	pthread_attr_setschedparam(&myatt[i], mypar);
}

int task_argument(void *arg)
{ 
struct 	task_par *tp;
	tp = (struct task_par *) arg;
	return tp->arg;
}

int task_period(int i)
{
	return tp[i].period;
}

void wait_for_task_end(int i)
{
	pthread_join(pend[i], NULL);
}

void set_period(int i) 
{
struct 	timespec t;

	clock_gettime(CLOCK_MONOTONIC, &t);
	time_copy(&(tp[i].at), t);
	time_copy(&(tp[i].dl), t);
	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].deadline);
}

// -----------------------------------------------------------------------------
// WAIT_FOR_PERIOD: It waits for the next period activation
// -----------------------------------------------------------------------------
void wait_for_period(int i)
{
	clock_nanosleep(CLOCK_MONOTONIC, TIMER_ABSTIME, &(tp[i].at), NULL);
	time_add_ms(&(tp[i].at), tp[i].period);
	time_add_ms(&(tp[i].dl), tp[i].period);
}

// -----------------------------------------------------------------------------
// DEADLINE_MISS: checks for deadline misses
// -----------------------------------------------------------------------------
int deadline_miss(int i) 
{
struct 	timespec now;

	clock_gettime(CLOCK_MONOTONIC, &now);
	if (time_cmp(now, tp[i].dl) > 0) {
		tp[i].dmiss++;
		return 1;
	}
	return 0;
}

// -----------------------------------------------------------------------------
// FUNCTIONS FOR TIME MANAGEMENT
// -----------------------------------------------------------------------------
void time_add_ms(struct timespec *t, int ms)
{
	t->tv_sec 	+= ms/1000;
	t->tv_nsec 	+= (ms%1000)*1000000;

	if (t->tv_nsec > 1000000000) {
		t->tv_nsec 	-= 1000000000;
		t->tv_sec 	+= 1;
	}
}

void time_copy(struct timespec *td, struct timespec ts) 
{
	td->tv_sec 	= ts.tv_sec;
	td->tv_nsec = ts.tv_nsec;
}

int time_cmp(struct timespec t1, struct timespec t2)
{
	if (t1.tv_sec > t2.tv_sec) 		return 1;
	if (t1.tv_sec < t2.tv_sec) 		return -1;
	if (t1.tv_nsec > t2.tv_nsec) 	return 1;
	if (t1.tv_nsec < t2.tv_nsec) 	return -1;
	return 0;
}
