#ifndef GRAPHLIB_H
#define GRAPHLIB_H

// -----------------------------------------------------------------------------
// TEXT RELATED CONSTANTS
// -----------------------------------------------------------------------------
#define LINE	80		// max length of a line of the .txt file
// -----------------------------------------------------------------------------
// TASK CONSTANTS
// -----------------------------------------------------------------------------
#define PER		15		// pendulum task period
#define DL		15		// pendulum task deadline
#define PPRI	80		// pendulum task priority
#define WPER	30		// deadline window period
#define WDL		30		// deadline window deadline
#define WPRI	50		// deadline window priority
// -----------------------------------------------------------------------------
// GRAPHICS CONSTANT
// -----------------------------------------------------------------------------
#define XWIN	800		// window x resolution
#define	YWIN	600		// window y resolution
#define	L_1P	600		// length of the image for max 1 pendulum 
#define	L_4P	300		// length of the image for max 4 pendulum
#define	L_9P	200		// length of the image for max 9 pendulum
// -----------------------------------------------------------------------------
// COLOUR
// -----------------------------------------------------------------------------
#define BKG		0		// background colour
#define WH 		15		// white colour
#define GREY 	8		// central pivot colour
#define BLUE 	1		// 1st mass colour
#define RED 	4		// 2nd mass colour
#define	GREEN	2		// pendulum colour
// -----------------------------------------------------------------------------
// STATUS WINDOW CONSTANT
// -----------------------------------------------------------------------------
#define LWIN	600		// left		status window coordinate
#define BWIN	599		// bottom 	status window coordinate
#define RWIN	799		// right	status window coordinate
#define UWIN	0		// upper 	status window coordinate
#define N_LINE	58		// number of line in the status window
#define TIT_P	40		// titol position
#define NUM_P	66		// position of the number of the page
#define	DIST_L	10		// vertical distance between lines
#define BEG_L	5		// orizontal position where the line begins	
// -----------------------------------------------------------------------------
// PENDULUM CONSTANT
// -----------------------------------------------------------------------------
#define MAX_DP	9					// max number of pendulums
#define R 		5					// radius of the mass dot
#define	PI 		3.141592653589793	// questo è necesario per avere una giusta approssimazione
#define TLEN	2					// trail point stored 
#define SH		5					// pixel radius of the trajectory shade

// -----------------------------------------------------------------------------
// STRUCTURE
// -----------------------------------------------------------------------------

// circular buffer to store the points of the trail
struct cbuf{
	int top;	// stores the index of the last element of the buffer
	int x[TLEN];
	int y[TLEN];
};

// point structure to be used to collocate the pendulum's bitmap 
struct point {
	int x;
	int y;
};

// structure to store the coordinates of the double pendulum in real dimension
struct point_real {
	float x1;
	float y1;
	float x2;
	float y2;
};

// structure to store the coordinates of the double pendulum in screen dimension
struct point_pixel {
	int x0;
	int y0;
	int x1;
	int y1;
	int x2;
	int y2;
};

// structure to contain the parameters of the pendulum
struct par {
	float 	l1;				// length 			of the 1st pendulum
	float 	m1;				// mass 			of the 1st pendulum
	double 	tht1;			// angle 			of the 1st pendulum
	float 	tht1_dot;		// angular velocity of the 1st pendulum
	float 	l2;				// length 			of the 2st pendulum
	float 	m2;				// mass 			of the 2st pendulum
	double 	tht2;			// angle 			of the 2st pendulum
	float 	tht2_dot;		// angular velocity of the 2st pendulum
};

// -----------------------------------------------------------------------------
// EXTERN VARIABLES
// -----------------------------------------------------------------------------
extern struct 	point_real 	pr[MAX_DP];
extern struct 	point_pixel pp[MAX_DP];
extern struct 	par 		pnd[MAX_DP];
extern int 					win;

// -----------------------------------------------------------------------------
// FILE MANAGEMENT FUNCTIONS
// -----------------------------------------------------------------------------
void 	read_data();
// -----------------------------------------------------------------------------
// GRAPHIC INITIALIZATION FUNCTIONS
// -----------------------------------------------------------------------------
void 	display_init();
// -----------------------------------------------------------------------------
// GRAPHIC MANAGEMENT FUNCTIONS
// -----------------------------------------------------------------------------
char 	get_scancode();
void 	execute_scan(char scan);
void 	draw_pend(int i);
void 	draw_dline();

#endif