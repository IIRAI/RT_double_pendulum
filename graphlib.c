/*
COMMENT HERE
*/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <allegro.h>
#include <math.h>
#include "graphlib.h"
#include "tasklib.h"
// -----------------------------------------------------------------------------
// FILE MANAGEMENT VARIABLES
// -----------------------------------------------------------------------------
FILE 	*fd;
int 	npend;							// number of double pendulum 
// int 	value[N_PAR * MAX_DP];			// contains the initial values 
//------------------------------------------------------------------------------
struct 		par 		pnd[MAX_DP];	
struct 		point 		ref[MAX_DP];	
struct 		point_real 	pr[MAX_DP];			
struct 		point_pixel pp[MAX_DP];			
struct 		cbuf 		trail[MAX_DP];	
int 		dline_flag;					// flag to draw the deadline window
int 		change_trail;				// flag to change the trajectory
float 		lag[MAX_DP];				// istantaneous Lagrangian value
// -----------------------------------------------------------------------------
// BITMAP VARIABLES
// -----------------------------------------------------------------------------
int 		l;							// length of the BITMAP to display the pendulums
BITMAP  	*bkg_pend[MAX_DP];			// background of the pendulum
BITMAP 		*light_traj[MAX_DP];		// contains the light 		mode trail 
BITMAP 		*lagr_traj[MAX_DP];			// contains the lagrangian 	mode trail
// BITMAP for the status window
BITMAP  	*keycmd;					// lists the keyboard inputs
BITMAP  	*init_data;					// lists the initial data
BITMAP  	*geom_data;					// lists the geometrical data
BITMAP  	*dline;						// lists the missed deadline
COLOR_MAP 	table;						// contains the shade values for light effects
// -----------------------------------------------------------------------------
// FILE MANAGEMENT FUNCTIONS
// -----------------------------------------------------------------------------
static void 	get_num_pend();
static void		get_theta_1();
static void		get_theta_2();
static void 	get_theta_dot_1();
static void 	get_theta_dot_2();
static void		get_length_1();
static void		get_length_2();
static void		get_mass_1();
static void		get_mass_2();
// -----------------------------------------------------------------------------
// GRAPHIC INITIALIZATION FUNCTIONS
// -----------------------------------------------------------------------------
static void 	pbox_init();
static void 	bkg_length();
static void		get_ref_1p();
static void		get_ref_4p();
static void		get_ref_9p();
static void 	create_bkg();
static void 	create_image(BITMAP **image, int w, int l);
static void 	pend_init();
static void 	status_win_init();
static void 	create_window_bitmap();
static void 	keycmd_page();
static void 	light_column();
static void 	lagr_column();
static void 	initdata_page();
static void 	geomdata_page();
static void 	deadline_page();
// -----------------------------------------------------------------------------
// GRAPHIC MANAGEMENT FUNCTIONS
// -----------------------------------------------------------------------------
static double 	deg2rad(float degree);
static float 	rad2deg(double radiant);
static void 	xy_real(int i);
static void 	xy_graph(int i);
static int 		real2graph(float value, int i);
// -----------------------------------------------------------------------------
// TRAJECTORY FUNCTIONS
// -----------------------------------------------------------------------------
static void 	store_trail(int i);
static void 	draw_trail(int i);
static void 	light_trail(int x2, int y2, int x1, int y1, int a);
static void 	lagr_trail(int x2, int y2, int x1, int y1, int a);
static float 	tan_tr(int i);
static float 	gauss(float var, float x);
static float 	lin_v (int i);
static float 	L(int i);

// -----------------------------------------------------------------------------
// FILE MANAGEMENT FUNCTIONS;
// -----------------------------------------------------------------------------
// READ_DATA: get the data from the .txt file
// -----------------------------------------------------------------------------
void read_data()
{
	get_num_pend();
	get_theta_1();
	get_theta_2();
	get_theta_dot_1();
	get_theta_dot_2();
	get_length_1();
	get_length_2();
	get_mass_1();
	get_mass_2();
}
// -----------------------------------------------------------------------------
// GET_NUM_PEND: searches and gets the number of the double pendulums 
// -----------------------------------------------------------------------------
void get_num_pend()
{
int 	start = 1;								// flag for the while loop
char 	*ref = "Number of Double Pendulums";	// line to be search
char 	row[LINE];								// stores the actual row	

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}			

	while (start) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) start = 0;

	}
	
	fgets(row, LINE, fd);
	npend = atoi(row);

	fclose(fd);

}

// -----------------------------------------------------------------------------
// GET_THETA_1: searches and gets the value of theta 1 for each double pendulum 
// -----------------------------------------------------------------------------

void get_theta_1()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Theta_1:";						// line to be search
char 	row[LINE];								// stores the actual row

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}					

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].tht1 = deg2rad(atoi(row));
			start++;
		}
	}

	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_THETA_2: searches and gets the value of theta 2 for each double pendulum 
// -----------------------------------------------------------------------------

void get_theta_2()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Theta_2:";						// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].tht2 = deg2rad(atoi(row));
			start++;
		}
	}
	
	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_THETA_DOT_1: searches and gets the value of the initial angular velocity
// of the first pendulum for each double pendulum 
// -----------------------------------------------------------------------------

void get_theta_dot_1()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Angular_velocity_1:";			// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].tht1_dot = (atoi(row));
			start++;
		}
	}
	
	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_THETA_DOT_2: searches and gets the value of the initial angular velocity
// of the second pendulum for each double pendulum 
// -----------------------------------------------------------------------------

void get_theta_dot_2()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Angular_velocity_2:";			// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].tht2_dot = (atoi(row));
			start++;
		}
	}
	
	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_LENGTH_1: searches and gets the length of the first pendulum for each
// double pendulum
// -----------------------------------------------------------------------------

void get_length_1()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Length_1:";						// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].l1 = atoi(row);
			start++;
		}
	}
	
	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_LENGTH_2: searches and gets the length of the second pendulum for each
// double pendulum
// -----------------------------------------------------------------------------

void get_length_2()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Length_2:";						// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].l2 = atoi(row);
			start++;
		}
	}

	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_MASS_1: searches and gets the mass of the first pendulum for each
// double pendulum
// -----------------------------------------------------------------------------

void get_mass_1()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Mass_1:";						// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].m1 = atoi(row);
			start++;
		}
	}

	fclose(fd);
}

// -----------------------------------------------------------------------------
// GET_MASS_2: searches and gets the mass of the second pendulum for each
// double pendulum
// -----------------------------------------------------------------------------

void get_mass_2()
{
int 	start = 0;								// flag for the while loop
char 	*ref = "Mass_2:";						// line to be search
char 	row[LINE];								// stores the actual row				

	fd = fopen("Data.txt", "r");	
	if (fd == NULL) {
		printf ("Error opening the file.\n");
		exit(1);
	}	

	while (start < npend) {

		fgets(row, LINE, fd);

		if (row == NULL) {
			printf("Error reading the file.\n");
		}

		if (strncmp(row, ref, (int)strlen(ref)) == 0) {
			fgets(row, LINE, fd);
			pnd[start].m2 = atoi(row);
			start++;
		}
	}
		fclose(fd);
}

// -----------------------------------------------------------------------------
// INIT_DISPLAY: initialize the display based on the number of pendulums
// -----------------------------------------------------------------------------

void display_init() 
{
	allegro_init();
	set_gfx_mode(GFX_AUTODETECT_WINDOWED, XWIN, YWIN, 0, 0);
	clear_to_color(screen, BKG);
	install_keyboard();

	// generate the table containing the shades of white
	create_light_table(&table, default_palette, 0, 0, 0, NULL);		

	pbox_init();				// draw the Pendulum's Box
	pend_init();				// initialize and draw the Pendulums
	status_win_init();					// draw the Data Window
}

// -----------------------------------------------------------------------------
// PBOX_INIT: initializes the boxes where the double pendulums are displayed
// -----------------------------------------------------------------------------

void pbox_init()
{
	bkg_length();
	create_bkg();
}

// -----------------------------------------------------------------------------
// BKG_LENGTH: calculates the side legth of the boxes, based on the number of
// double pendulums as input
// -----------------------------------------------------------------------------

void bkg_length()
{
	if (npend == 1) {
		l = L_1P;
		get_ref_1p();
	}		

	if (npend >= 2 && npend <= 4) {
		l = L_4P;
		get_ref_4p();
	} 	

	if (npend >= 5 && npend <= 9) {
		l = L_9P;
		get_ref_9p();
	}	
}

// -----------------------------------------------------------------------------
// GET_REF_1P; GET_REF_4P; GET_REF_9P: 
// defines the coordinates of the boxes, based on the number of double pendulums
// respectively for 1 pendulum, for 4 pendulums at most and 9 pendulums at most. 
// -----------------------------------------------------------------------------

void get_ref_1p()
{
	ref[0].x = 0;
	ref[0].y = 0;
}

void get_ref_4p()
{
	ref[0].x = 0;
	ref[0].y = 0;
	ref[1].x = l;
	ref[1].y = 0;
	ref[2].x = 0; 
	ref[2].y = l;
	ref[3].x = l;
	ref[3].y = l;	

}

void get_ref_9p()
{
	ref[0].x = 0;
	ref[0].y = 0;
	ref[1].x = l;
	ref[1].y = 0;
	ref[2].x = 0; 
	ref[2].y = l;
	ref[3].x = l;
	ref[3].y = l;
	ref[2].x = l * 2;
	ref[2].y = 0;
	ref[3].x = 0;
	ref[3].y = l;
	ref[4].x = l;
	ref[4].y = l;
	ref[5].x = l * 2;
	ref[5].y = l;
	ref[6].x = 0;
	ref[6].y = l * 2;
	ref[7].x = l;
	ref[7].y = l * 2;
	ref[8].x = l * 2;
	ref[8].y = l * 2;
}

// -----------------------------------------------------------------------------
// CREATE_BKG: creates all the BITMAP needed to display each double pendulum:
// -----------------------------------------------------------------------------

void create_bkg()
{
int 	i;

	for (i = 0; i < npend; i++) {
 		create_image(&bkg_pend[i], l, l);		// to display the pendulum
		create_image(&light_traj[i], l, l);		// to store the Light trail
		create_image(&lagr_traj[i], l, l);		// to store the Lagrangian trail
	}
}

void create_image(BITMAP **image, int w, int l)
{
		*image = create_bitmap(w, l);
		clear_bitmap(*image);
		rect(*image, 0, 0, w - 1, l - 1, BLUE);
}

// -----------------------------------------------------------------------------
// PEND_INIT: stores the first values inside the cycle buffer to draw the trail
// and draws the double pendulums in their initial position
// -----------------------------------------------------------------------------

void pend_init()
{
int 	i = 0;
int 	n;

	for (n = 0; n < npend * 6; n += 6) {
			// printf("grad1:   %2.15lf\n", pnd[i].tht1);
			// printf("grad1.1: %2.30lf\n", cos(pnd[i].tht1));
			// printf("grad1.2: %2.30lf\n", sin(pnd[i].tht1));
			// printf("grad2:   %2.15lf\n", pnd[i].tht2);
			// printf("grad2.1: %2.30lf\n", cos(pnd[i].tht1));
			// printf("grad2.3: %2.30lf\n", sin(pnd[i].tht1));
		xy_real(i);		
		xy_graph(i);
		store_trail(i);
		draw_pend(i);
		i ++;
	}
}

// -----------------------------------------------------------------------------
// STATUS_WIN_INIT: initializes the bitmaps to report the information about
// the program.
// -----------------------------------------------------------------------------			

void status_win_init()
{
	create_window_bitmap();
	keycmd_page();
	initdata_page();
	geomdata_page();
	deadline_page();
}

// -----------------------------------------------------------------------------
// CREATE_WINDOW_BITMAP: creates the bitmaps to be used for the status window
// -----------------------------------------------------------------------------

void create_window_bitmap()
{
	create_image(&keycmd, 		RWIN-LWIN+1, BWIN - UWIN + 1);
	create_image(&init_data, 	RWIN-LWIN+1, BWIN - UWIN + 1);
	create_image(&geom_data,	RWIN-LWIN+1, BWIN - UWIN + 1);
	create_image(&dline, 		RWIN-LWIN+1, BWIN - UWIN + 1);
}

// -----------------------------------------------------------------------------
// KEYCMD_PAGE: writes the keyboard command and the legend about the two 
// different trail used in the program
// -----------------------------------------------------------------------------

void keycmd_page()
{
	textout_ex(keycmd, font, "Command Window", 40, 10, WH, BKG);

	textout_ex(keycmd, font, "Press ENTER:", 5, 30, WH, BKG);
	textout_ex(keycmd, font, "run the simulation", 5, 40, WH, BKG);

	textout_ex(keycmd, font, "Press SPACE:", 5, 60, WH, BKG);
	textout_ex(keycmd, font, "change the trail", 5, 70, WH, BKG);

	textout_ex(keycmd, font, "Press 1:", 5, 90, WH, BKG);
	textout_ex(keycmd, font, "turn to Command Window", 5, 100, WH, BKG);

	textout_ex(keycmd, font, "Press 2:", 5, 120, WH, BKG);
	textout_ex(keycmd, font, "turn to Initial Data", 5, 130, WH, BKG);
	textout_ex(keycmd, font, "Window", 5, 140, WH, BKG);

	textout_ex(keycmd, font, "Press 3:", 5, 160, WH, BKG);
	textout_ex(keycmd, font, "turn to Geometrical", 5, 170, WH, BKG);
	textout_ex(keycmd, font, "Data Window", 5, 180, WH, BKG);

	textout_ex(keycmd, font, "Press 4:", 5, 200, WH, BKG);
	textout_ex(keycmd, font, "turn to Deadline Window", 5, 210, WH, BKG);

	textout_ex(keycmd, font, "Press ESC:", 5, 230, WH, BKG);
	textout_ex(keycmd, font, "exit the program", 5, 240, WH, BKG);

	textout_ex(keycmd, font, "Types of trail;", 5, 260, WH, BKG);
	textout_ex(keycmd, font, "left:  proportional", 5, 270, WH, BKG);
	textout_ex(keycmd, font, "to the linear velocity", 5, 280, WH, BKG);
	textout_ex(keycmd, font, "right: proportional", 5, 290, WH, BKG);
	textout_ex(keycmd, font, "to the Lagrangian value", 5, 300, WH, BKG);

	textout_ex(keycmd, font, "page 1/4", 66, BWIN - 15, WH, BKG);

	light_column();
	lagr_column();

	blit(keycmd, screen, 0, 0, LWIN, UWIN, keycmd->w, keycmd->h);
}

// -----------------------------------------------------------------------------
// LIGHT_COLUMN: draws a column in the command page to show the shades of the
// light trail adopted in the program
// -----------------------------------------------------------------------------

void light_column()
{
int 	i;
int 	a;
float 	shade[SH];

	for (i = 0; i < 255; i++) {

		for (a = 0; a < SH; a++) {

			shade[a] = gauss( (0.0005 * (i+1)), (0.05 * a)) / gauss( (0.0005 * (i+1)), 0);

			putpixel(keycmd, 66 + a, 320 + i, table.data[(int)(255 * shade[a])][WH]);
			putpixel(keycmd, 66 - a, 320 + i, table.data[(int)(255 * shade[a])][WH]);
		}
	}
}

// -----------------------------------------------------------------------------
// LAGR_COLUMN: draws a column in the command page to show the shade of colours
// of the lagrangian trail adopted in the program
// -----------------------------------------------------------------------------

void lagr_column()
{
int 	i;				// index of the for cycle
int 	a;				// index of the internal for cycle
float 	shade[SH];
int 	col;
int 	r, g, b;

	for (i = 0; i < SH; i++) {
		shade[i] = pow(1.5, -i);
	}

	for (i = 0; i < 255; i++) {

		hsv_to_rgb( i, 1, 1,   &r, &g, &b);		
		col = makecol(r, g, b);

		for (a = 0; a < SH; a++) {
			putpixel(keycmd, 132 + a, 320 + i, table.data[(int)(255 * shade[a])][col]);
			putpixel(keycmd, 132 - a, 320 + i, table.data[(int)(255 * shade[a])][col]);
		}
	}
}

// ------------------------------------------------------------------------------
// INITDATA_PAGE: writes the initial position and velcity of each double pendulum 
// ------------------------------------------------------------------------------

void initdata_page()
{	
int 	i;
char 	d[30];

	textout_ex(init_data, font, "Initial Data", 40, 10, WH, BKG);

	for (i=0; i<npend; i++) {

		sprintf(d, "Pendulum %d", i+1);
		textout_ex(init_data, font, d, 5, 30 + (60 * i), WH, BKG);

		sprintf(d, "   Theta1:  %g°", rad2deg(pnd[i].tht1));
		textout_ex(init_data, font, d, 5, 40 + (60 * i), WH, BKG);

		sprintf(d, "   Omega1:  %g [rad/s]", pnd[i].tht1_dot);
		textout_ex(init_data, font, d, 5, 50 + (60 * i), WH, BKG);

		sprintf(d, "   Theta2:  %g°", rad2deg(pnd[i].tht2));
		textout_ex(init_data, font, d, 5, 60 + (60 * i), WH, BKG);

		sprintf(d, "   Omega2:  %g [rad/s]", pnd[i].tht2_dot );
		textout_ex(init_data, font, d, 5, 70 + (60 * i), WH, BKG);
	}

	textout_ex(init_data, font, "page 2/4", 66, BWIN - 15, WH, BKG);
}

// -----------------------------------------------------------------------------
// GEOMDATA_PAGE: writes the geometrical data of each double pendulum 
// as length and mass
// -----------------------------------------------------------------------------

void geomdata_page()
{
int 	i;
char 	d[30];

	textout_ex(geom_data, font, "Geometrical Data", 40, 10, WH, BKG);

	for (i=0; i<npend; i++) {

		sprintf(d, "Pendulum %d", i+1);
		textout_ex(geom_data, font, d, 5, 30 + (60 * i), WH, BKG);

		sprintf(d, "   Length1: %g", pnd[i].l1 );
		textout_ex(geom_data, font, d, 5, 40 + (60 * i), WH, BKG);

		sprintf(d, "   Mass1:   %g", pnd[i].m1 );
		textout_ex(geom_data, font, d, 5, 50 + (60 * i), WH, BKG);

		sprintf(d, "   Length2: %g", pnd[i].l2 );
		textout_ex(geom_data, font, d, 5, 60 + (60 * i), WH, BKG);

		sprintf(d, "   Mass2:   %g", pnd[i].m2 );
		textout_ex(geom_data, font, d, 5, 70 + (60 * i), WH, BKG);

	}

	textout_ex(geom_data, font, "page 3/4", 66, BWIN - 15, WH, BKG);
}

// -----------------------------------------------------------------------------
// DEADLINE_PAGE: writes the deadline missed of every running thread 
// -----------------------------------------------------------------------------

void deadline_page()
{
char 	w[30];

	dline_flag = 0;

	textout_ex(dline, font, "Deadline Window", 40, 10, WH, BKG);

	task_create(npend, line_wintask, WPER, WDL, WPRI);

	textout_ex(dline, font, "page 4/4", 66, BWIN - 15, WH, BKG);

	if (dline_flag == 1) {
		blit(dline, screen, 0, 0, LWIN, UWIN, dline->w, dline->h);
	}
}

// -----------------------------------------------------------------------------
// GET_SCANCODE: returns the scancode of a pressed key
// -----------------------------------------------------------------------------

char get_scancode()
{
	if (keypressed()) 
		return readkey() >> 8;
	else return 0;
}

// -----------------------------------------------------------------------------
// EXECUTE_SCAN: execute different functions based on the key pressed
// -----------------------------------------------------------------------------

void execute_scan(char scan)
{
int 	n;
	switch(scan) {
		case KEY_ENTER:
			for (n = 0; n < npend; n++) {
				task_create(n, pend_task, PER, DL, PPRI);
			}
			break;
		case KEY_SPACE:
			if (change_trail == 0) change_trail = 1;
			else change_trail = 0;
			break;
		case KEY_1:
			dline_flag = 0;
			blit(keycmd, screen, 0, 0, LWIN, UWIN, keycmd->w, keycmd->h);
			break;
		case KEY_2:
			dline_flag = 0;
			blit(init_data, screen, 0, 0, LWIN, UWIN, init_data->w, init_data->h);
			break;
		case KEY_3:
			dline_flag = 0;
			blit(geom_data, screen, 0, 0, LWIN, UWIN, geom_data->w, geom_data->h);
			break;
		case KEY_4:
			dline_flag = 1;
			blit(dline, screen, 0, 0, LWIN, UWIN, dline->w, dline->h);
			break;
		default: 
			break;
	}
}

// -----------------------------------------------------------------------------
// DEG2RAD: converts values from degree  to radiant
// -----------------------------------------------------------------------------

double deg2rad(float degree)
{
double 	deg; 
double 	rad;

	deg = degree;
	rad = (deg / 180.0) * PI;

	return rad;
}

// -----------------------------------------------------------------------------
// RAD2DEG: converts values from radiant to degree
// -----------------------------------------------------------------------------

float rad2deg(double radiant)
{
float 	deg;
float 	rad;

	rad = radiant;
	deg = (rad / PI) * 180.0;

	return deg;
}

// -----------------------------------------------------------------------------
// XY_REAL: 
// evaluates the X and Y coordinates with respect to the real world reference
// -----------------------------------------------------------------------------

void xy_real(int i)
{
float 	l1;
float 	l2;
float 	s_tht1;
float 	c_tht1;
float 	s_tht2;
float 	c_tht2;

	l1 		= pnd[i].l1;
	l2 		= pnd[i].l2;
	s_tht1 	= sin(pnd[i].tht1);
	c_tht1 	= cos(pnd[i].tht1);
	s_tht2 	= sin(pnd[i].tht2);
	c_tht2 	= cos(pnd[i].tht2);

	pr[i].x1 = l1 * s_tht1;
	pr[i].y1 = l1 * c_tht1;
	pr[i].x2 = pr[i].x1 + l2 * s_tht2;
	pr[i].y2 = pr[i].y1 + l2 * c_tht2;
}

// -----------------------------------------------------------------------------
// XY_GRAPH:
// evaluates the X and Y coordinates with respect to the display reference
// -----------------------------------------------------------------------------

void xy_graph(int i)
{
	pp[i].x0 = (l / 2);
	pp[i].y0 = (l / 2);
	pp[i].x1 = pp[i].x0 + real2graph(pr[i].x1, i);
	pp[i].y1 = pp[i].y0 + real2graph(pr[i].y1, i);
	pp[i].x2 = pp[i].x0 + real2graph(pr[i].x2, i);
	pp[i].y2 = pp[i].y0 + real2graph(pr[i].y2, i);
}

// -----------------------------------------------------------------------------
// REAL2GRAPH: 
// converts the coordinates from the the real world reference 
// to the display reference
// -----------------------------------------------------------------------------

int real2graph(float value, int i)
{
int 	pixel;		// value of the graphic display [pixel]
float 	l12;		// sum of the length of the two pendulums
	
	l12 = pnd[i].l1 + pnd[i].l2;
	pixel = (int) ((value / l12) * (float)((l / 2) - R));
	return pixel;
}

/******************************************************************************
 The next two functions are executed inside a thread 
 and manage the graphical aspect
*******************************************************************************/

// -----------------------------------------------------------------------------
// DRAW_PEND: draw the double pendulum on the backround bitmap
// -----------------------------------------------------------------------------

void draw_pend(int i)
{
	xy_real(i);
	xy_graph(i);

	store_trail(i);
	draw_trail(i);

	line(bkg_pend[i], pp[i].x0, pp[i].y0, pp[i].x1, pp[i].y1, GREEN);
	line(bkg_pend[i], pp[i].x1, pp[i].y1, pp[i].x2, pp[i].y2, GREEN);
	circlefill(bkg_pend[i], pp[i].x0, pp[i].y0, R, GREY);
	circlefill(bkg_pend[i], pp[i].x1, pp[i].y1, R, BLUE);
	circlefill(bkg_pend[i], pp[i].x2, pp[i].y2, R, RED);

	blit(bkg_pend[i], screen, 0, 0, ref[i].x, ref[i].y, bkg_pend[i]->w, bkg_pend[i]->h);
}

// -----------------------------------------------------------------------------
// DRAW_DLINE: draw the deadline missed of the running threads
// -----------------------------------------------------------------------------

void draw_dline()
{
int 	i;				// "for" index
char 	w[30];

	for(i = 0; i < npend; i++) {
			sprintf(w, "Pendulum %d: %d", i + 1, tp[i].dmiss);
			textout_ex(dline, font, w, 5, 30 + (20 * i), WH, BKG);
	}

	sprintf(w, "Window task: %d", tp[npend].dmiss);
	textout_ex(dline, font, w, 5, 30 + (20 * npend), WH, BKG);
}

// -----------------------------------------------------------------------------
// STORE_TRAIL:
// stores the new position of the second pendulum in the cycle buffer
// -----------------------------------------------------------------------------

void store_trail(int i)
{
int 	k;

	k = trail[i].top;			// k takes the value of top
	trail[i].x[k] = pp[i].x2;
	trail[i].y[k] = pp[i].y2;
	k = (k + 1) % TLEN;			// k is updated to the next element
	trail[i].top = k;			// top is updated to the last element
}

// -----------------------------------------------------------------------------
// DRAW_TRAIL: 
// draw both trail (Light and Lagrangian) inside the respective bitmap
// -----------------------------------------------------------------------------

void draw_trail(int i)
{
int 	k;
int 	k_old;					// trail indexes
int 	x1, y1;					// graphics coordinates
int 	x2, y2;					// graphics coordinates

	k = trail[i].top;			// get the point indicated by j
	k_old = (k + 1) % TLEN;
	x1 = trail[i].x[k_old];
	y1 = trail[i].y[k_old];
	x2 = trail[i].x[k];
	y2 = trail[i].y[k];

	light_trail(x2, y2, x1, y1, i);	
	lagr_trail(x2, y2, x1, y1, i);

	if (change_trail == 0) {
		blit(light_traj[i], bkg_pend[i], 0, 0, 0, 0, light_traj[i]->w, light_traj[i]->h);
	}
	else {
		blit(lagr_traj[i], bkg_pend[i], 0, 0, 0, 0, lagr_traj[i]->w, lagr_traj[i]->h);
	}
}

// -----------------------------------------------------------------------------
// LIGHT_TRAIL: manages the representation of the Light trail
// -----------------------------------------------------------------------------

void light_trail(int x2, int y2, int x1, int y1, int i)
{
int 	a;
double 	shade[SH];		// array contenente il grado di sfumatura
float 	v;
float 	alpha;
float 	c_a;
float 	s_a;

	v = lin_v(i);

	alpha 	= tan_tr(i);
	c_a 	= cos(alpha);
	s_a 	= sin(alpha);

	for (a = 0; a < SH; a++) {

		shade[a] = gauss( fabs(1 / ((v) + 1)), (0.06 * a)) / gauss( fabs(1 / ((v) + 1)), 0);
	}

	for(a = 4; a >= 0; a--) {		// dal più scuro al più chiaro
		line(light_traj[i], x2 + (int)(a * s_a), y2 + (int)(a * c_a),  x1 + (int)(a * s_a), y1 + (int)(a * c_a), table.data[(int)(255 * shade[a])][WH]);
		line(light_traj[i], x2 - (int)(a * s_a), y2 - (int)(a * c_a),  x1 - (int)(a * s_a), y1 - (int)(a * c_a), table.data[(int)(255 * shade[a])][WH]);

		line(light_traj[i], x2 + (int)(a * s_a), y2 + (int)((a - 1) * c_a),  x1 + (int)(a * s_a), y1 + (int)((a - 1) * c_a), table.data[(int)(255 * shade[a])][WH]);
		line(light_traj[i], x2 - (int)(a * s_a), y2 - (int)((a - 1) * c_a),  x1 - (int)(a * s_a), y1 - (int)((a - 1) * c_a), table.data[(int)(255 * shade[a])][WH]);		
	}
	line(light_traj[i], x2, y2,  x1, y1, table.data[255][WH]);
}

// -----------------------------------------------------------------------------
// LAGR_TRAIL: manages the representation of the Lagrangian trail
// -----------------------------------------------------------------------------

void lagr_trail(int x2, int y2, int x1, int y1, int i)
{
int 	a;
double 	shade[SH];		// array contenente il grado di sfumatura
float 	alpha;
float 	c_a;
float 	s_a;
int 	col;
int 	r, g, b;

	hsv_to_rgb( (int)(L(i)) % 255, 1, 1,   &r, &g, &b);		
	col = makecol(r, g, b);

	alpha = tan_tr(i);
	c_a = cos(alpha);
	s_a = sin(alpha);

	for (a = 0; a < SH; a++) {
		shade[a] = pow(1.5, -a);			// proporzione esponenziale 
	}

	for(a = 4; a >= 0; a--) {		// dal più scuro al più chiaro
		line(lagr_traj[i], x2 + (int)(a * s_a), y2 + (int)(a * c_a),  x1 + (int)(a * s_a), y1 + (int)(a * c_a), table.data[(int)(255 * shade[a])][col]);
		line(lagr_traj[i], x2 - (int)(a * s_a), y2 - (int)(a * c_a),  x1 - (int)(a * s_a), y1 - (int)(a * c_a), table.data[(int)(255 * shade[a])][col]);

		line(lagr_traj[i], x2 + (int)(a * s_a), y2 + (int)((a - 1) * c_a),  x1 + (int)(a * s_a), y1 + (int)((a - 1) * c_a), table.data[(int)(255 * shade[a])][col]);
		line(lagr_traj[i], x2 - (int)(a * s_a), y2 - (int)((a - 1) * c_a),  x1 - (int)(a * s_a), y1 - (int)((a - 1) * c_a), table.data[(int)(255 * shade[a])][col]);		
	}
	line(lagr_traj[i], x2, y2,  x1, y1, table.data[200][col]);	
}

// -----------------------------------------------------------------------------
// TAN_TR: 
// computes the istant tangent angle of the trajectory of the second pendulum
// -----------------------------------------------------------------------------

float tan_tr(int i)
{
float 	alpha;					// istant tangent angle 
float 	tht1, tht2;
float 	tht1_dot, tht2_dot;
float 	l1, l2;
float 	c_th1, c_th2;
float 	s_th1, s_th2;
float 	x_dot, y_dot;

	tht1 		= pnd[i].tht1;
	tht2 		= pnd[i].tht2;
	tht1_dot 	= pnd[i].tht1_dot;
	tht2_dot 	= pnd[i].tht2_dot;
	l1 			= pnd[i].l1;
	l2 			= pnd[i].l2;

	c_th1 = cos(tht1);
	c_th2 = cos(tht2);
	s_th1 = sin(tht1);
	s_th2 = sin(tht2);

	x_dot = (l1 * tht1_dot * c_th1) + (l2 * tht2_dot * c_th2);
	y_dot = (l1 * tht1_dot * s_th1) + (l2 * tht2_dot * s_th2);

	alpha = atan2(y_dot, x_dot);

	return alpha;
}

// -----------------------------------------------------------------------------
// GAUSS: 
// computes the value of a Gaussian function.
// Variables:
// 	"var" stores the variance 
// 	"x" stores the value at which the Gaussian function is evaluated
// -----------------------------------------------------------------------------

float gauss(float var, float x)
{
float 	c, e;		// constant and exponent of the gaussian function

	c 	= 1 / (var * sqrt(2 * PI));
	e = -(x * x / (2 * var * var));

	return c * exp(e);
}

// -----------------------------------------------------------------------------
// LIN_V: evaluates the istantaneous linear velocity of the second mass attached
// at the end of the second pendulum
// -----------------------------------------------------------------------------

float lin_v (int i)		
{
float 	vel;
float 	tht1, tht2;
float 	tht1_dot, tht2_dot;
float 	l1, l2;
float 	c_th1th2;
float 	N1, N2, N3;

	tht1 		= pnd[i].tht1;
	tht2 		= pnd[i].tht2;
	tht1_dot 	= pnd[i].tht1_dot;
	tht2_dot 	= pnd[i].tht2_dot;
	l1 			= pnd[i].l1;
	l2 			= pnd[i].l2;

	c_th1th2 = cos(tht1 - tht2);

	N1 = pow( (l1 * tht1_dot), 2);
	N2 = pow( (l2 * tht2_dot), 2);
	N3 = 2 * l1 * l2 * tht1_dot * tht2_dot * c_th1th2;

	vel = sqrt( N1 + N2 + N3);
	
	return vel;
}

// -----------------------------------------------------------------------------
// L: evaluates the istantaneous value of the Lagrangian of the system
// -----------------------------------------------------------------------------

float L(int i)
{
float 	tht1, tht2;
float 	tht1_dot;
float 	v2;
float 	l1, l2;
float 	m1, m2;
float 	c_tht1, c_tht2;
float 	T1, T2;
float 	V1, V2;

	tht1 		= pnd[i].tht1;
	tht2 		= pnd[i].tht2;
	tht1_dot 	= pnd[i].tht1_dot;
	v2 			= lin_v(i);
	l1 			= pnd[i].l1;
	l2 			= pnd[i].l2;
	m1 			= pnd[i].m1;
	m2 			= pnd[i].m2;
	c_tht1 		= cos(tht1);
	c_tht2 		= cos(tht2);

	T1 = 0.5 * m1 * pow((l1 * tht1_dot), 2);
	T2 = 0.5 * m2 * pow(v2, 2);
	V1 = - m1 * G * l1 * c_tht1;
	V2 = - m2 * G * ((l1 * c_tht1) + (l2 * c_tht2));

	lag[i] = T1 + T2 - V1 - V2;
	
	return lag[i];
}