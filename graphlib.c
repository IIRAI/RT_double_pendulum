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
// VARIABLES
// -----------------------------------------------------------------------------
FILE 		*fd;
int 		npend;						// number of double pendulum 
struct 		par 		pnd[MAX_DP];	
struct 		point 		ref[MAX_DP];	
struct 		point_real 	pr[MAX_DP];			
struct 		point_pixel pp[MAX_DP];			
struct 		cbuf 		trail[MAX_DP];	
int 		dline_flag = 0;				// flag to draw the deadline window
int 		change_trail = 0;			// flag to change the trajectory
// -----------------------------------------------------------------------------
// BITMAP VARIABLES
// -----------------------------------------------------------------------------
int 		l;							// length of the BITMAP to display the pendulums
BITMAP  	*bkg_pend[MAX_DP];			// background of the pendulum
BITMAP 		*light_tr[MAX_DP];			// contains the Light trail 
BITMAP 		*lagr_tr[MAX_DP];			// contains the Lagrangian trail
// BITMAP for the status window
int 		pos[N_LINE];				// vertical position of the lines in the status window
BITMAP  	*keycmd;					// lists the keyboard inputs
BITMAP  	*init_data;					// lists the initial data
BITMAP  	*geom_data;					// lists the geometrical data
BITMAP  	*dline;						// lists the missed deadline
COLOR_MAP 	table;						// contains the shade values for light effects
// -----------------------------------------------------------------------------
// FILE MANAGEMENT FUNCTIONS
// -----------------------------------------------------------------------------
static void 	search_num_pend(char *row);
static void 	search_theta_1(char *row);
static void 	search_theta_2(char *row);
static void 	search_theta_1_dot(char *row);
static void 	search_theta_2_dot(char *row);
static void 	search_length_1(char *row);
static void 	search_length_2(char *row);
static void 	search_mass_1(char *row);
static void 	search_mass_2(char *row);
// -----------------------------------------------------------------------------
// GRAPHIC INITIALIZATION FUNCTIONS
// -----------------------------------------------------------------------------
static void 	pbox_init();
static void 	bkg_length();
static void		get_ref_1p();
static void		get_ref_4p();
static void		get_ref_9p();
static void 	bkg_init();
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
static double 	deg2rad(double degree);
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
static float 	norm_gauss(float var, float x);
static float 	lin_v (int i);
static float 	L(int i);

// -----------------------------------------------------------------------------
// FILE MANAGEMENT FUNCTIONS;
// -----------------------------------------------------------------------------
// READ_DATA: get the data from the .txt file
// -----------------------------------------------------------------------------
void read_data()
{
char 	row[LINE];
	
	fd = fopen("Data.txt", "r");	
		if (fd == NULL) {
			printf ("Error opening the file.\n");
			exit(1);
		}

	search_num_pend(row);

	do {
		fgets(row, LINE, fd);
		search_theta_1(row);
		search_theta_2(row);
		search_theta_1_dot(row);
		search_theta_2_dot(row);
		search_length_1(row);
		search_length_2(row);
		search_mass_1(row);
		search_mass_2(row);
	}	while (!feof(fd));

	fclose(fd);
}

void search_num_pend(char *row)
{
int 	start = 1;								// flag for the while loop
char 	*ref = "Number of Double Pendulums";	// line to be search

	while (start) {
		fgets(row, LINE, fd);
		if (strncmp(row, ref, (int)strlen(ref)) == 0) start = 0;
	}
	
	fgets(row, LINE, fd);
	npend = atoi(row);

	if (npend > MAX_DP) {
		printf("Error: too many double pendulum\n");
		exit(0);
	}
}

void search_theta_1(char *row)
{
static int 	index = 0;
char 		*ref = "Theta_1:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].tht1 = deg2rad(atoi(row));
		index++;
	}
	else return;
}

void search_theta_2(char *row)
{
static int 	index = 0;
char 		*ref = "Theta_2:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].tht2 = deg2rad(atoi(row));
		index++;
	}
	else return;
}

void search_theta_1_dot(char *row)
{
static int 	index = 0;
char 		*ref = "Angular_velocity_1:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].tht1_dot = atoi(row);
		index++;
	}
	else return;
}

void search_theta_2_dot(char *row)
{
static int 	index = 0;
char 		*ref = "Angular_velocity_1:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].tht2_dot = atoi(row);
		index++;
	}
	else return;
}

void search_length_1(char *row)
{
static int 	index = 0;
char 		*ref = "Length_1:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].l1  = atoi(row);
		index++;
	}
	else return;
}

void search_length_2(char *row)
{
static int 	index = 0;
char 		*ref = "Length_2:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].l2  = atoi(row);
		index++;
	}
	else return;
}

void search_mass_1(char *row)
{
static int 	index = 0;
char 		*ref = "Mass_1:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].m1  = atoi(row);
		index++;
	}
	else return;
}

void search_mass_2(char *row)
{
static int 	index = 0;
char 		*ref = "Mass_2:";		// line to be searched

	if (strncmp(row, ref, (int)strlen(ref)) == 0) {
		fgets(row, LINE, fd);
		pnd[index].m2  = atoi(row);
		index++;
	}
	else return;
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
	status_win_init();			// draw the Data Window
}

// -----------------------------------------------------------------------------
// PBOX_INIT: initializes the boxes where the double pendulums are displayed
// -----------------------------------------------------------------------------
void pbox_init()
{
	bkg_length();
	bkg_init();
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

	if (npend >= 5) {
		l = L_9P;
		get_ref_9p();
	}	
}

// -----------------------------------------------------------------------------
// GET_REF_1P; GET_REF_4P; GET_REF_9P: 
// defines the reference coordinates, stored in the ref[] array,  of the Bitmaps,
// where the pendulums will be displayed.
// There are three functions, each working for different set of double pendulums:
// respectively for 1 pendulum, for 4 pendulums at most and 9 pendulums at most. 
// -----------------------------------------------------------------------------
void get_ref_1p()
{
int 	i; 			// for cycle index

	ref[0].x = 0;
	ref[0].y = 0;

	// initializes the other values of the arrays to a default value
	for (i = 1; i < MAX_DP; i++) {
		ref[i].x = 0;
		ref[i].y = 0;
	}
}

void get_ref_4p()
{
int 	def;			// for cycle index, of the default values of the array
int 	row;			// for cycle index, represents the row number 
int 	clm;			// for cycle index, representes the column number
int 	ind = 0;		// index of the box's reference defined

	// creates a checked pattern
	for (row = 0; row < 2; row ++) {
		for (clm = 0; clm < 2; clm ++) {
			ref[ind].x = clm * l;
			ref[ind].y = row * l;
			ind ++; 
		}
	}

	// initializes the other values of the arrays to a default value
	for (def = 4; def < MAX_DP; def++) {
		ref[def].x = 0;
		ref[def].y = 0;
	}
}

void get_ref_9p()
{
int 	def;			// for cycle index, of the default values of the array
int 	row;			// for cycle index, represents the row number 
int 	clm;			// for cycle index, representes the column number
int 	ind = 0;		// index of the box's reference defined

	// creates a checked pattern
	for (row = 0; row < 3; row++) {
		for (clm = 0; clm < 3; clm++) {
			ref[ind].x = clm * l;
			ref[ind].y = row * l;
			ind ++; 
		}
	}
}

// -----------------------------------------------------------------------------
// BKG_INIT: creates all the BITMAP needed to display each double pendulum
// -----------------------------------------------------------------------------
void bkg_init()
{
int 	i;		// for cycle index

	for (i = 0; i < npend; i++) {
 		create_image(&bkg_pend[i], l, l);
		create_image(&light_tr[i], l, l);
		create_image(&lagr_tr[i], l, l);
	}
}

void create_image(BITMAP **image, int w, int l)
{
	*image = create_bitmap(w, l);
	clear_bitmap(*image);
	rect(*image, 0, 0, w - 1, l - 1, BLUE);
}

// -----------------------------------------------------------------------------
// PEND_INIT: stores the first values inside the circular buffer to draw the 
// trail and draws the double pendulums in their initial position
// -----------------------------------------------------------------------------
void pend_init()
{
int 	i;		// for cycle index

	for (i = 0; i < npend; i++) {
		xy_real(i);		
		xy_graph(i);
		store_trail(i);
		draw_pend(i);
	}
}

// -----------------------------------------------------------------------------
// STATUS_WIN_INIT: initializes the bitmaps to report the information about
// the program.
// -----------------------------------------------------------------------------			
void status_win_init()
{
int 	i;		// for cycle index

	// initialize the position value of the lines in the status window
	for (i = 0; i < N_LINE; i++) {
		pos[i] = DIST_L * (i + 1);
	}

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
	create_image(&keycmd, 	 RWIN - LWIN + 1, BWIN - UWIN + 1);
	create_image(&init_data, RWIN - LWIN + 1, BWIN - UWIN + 1);
	create_image(&geom_data, RWIN - LWIN + 1, BWIN - UWIN + 1);
	create_image(&dline, 	 RWIN - LWIN + 1, BWIN - UWIN + 1);
}

// -----------------------------------------------------------------------------
// KEYCMD_PAGE: writes the keyboard command and the legend about the two 
// different trail used in the program
// -----------------------------------------------------------------------------
void keycmd_page()
{
	textout_ex(keycmd, font, "Command Window", TIT_P, pos[0], WH, BKG);

	textout_ex(keycmd, font, "Press ENTER:", BEG_L, pos[2], WH, BKG);
	textout_ex(keycmd, font, "runs the simulation", BEG_L, pos[3], WH, BKG);

	textout_ex(keycmd, font, "Press SPACE:", BEG_L, pos[5], WH, BKG);
	textout_ex(keycmd, font, "changes the trail", BEG_L, pos[6], WH, BKG);

	textout_ex(keycmd, font, "Press 1:", BEG_L, pos[8], WH, BKG);
	textout_ex(keycmd, font, "turns to Command Window", BEG_L, pos[9], WH, BKG);

	textout_ex(keycmd, font, "Press 2:", BEG_L, pos[11], WH, BKG);
	textout_ex(keycmd, font, "turns to Initial Data", BEG_L, pos[12], WH, BKG);
	textout_ex(keycmd, font, "Window", BEG_L, pos[13], WH, BKG);

	textout_ex(keycmd, font, "Press 3:", BEG_L, pos[15], WH, BKG);
	textout_ex(keycmd, font, "turns to Geometrical", BEG_L, pos[16], WH, BKG);
	textout_ex(keycmd, font, "Data Window", BEG_L, pos[17], WH, BKG);

	textout_ex(keycmd, font, "Press 4:", BEG_L, pos[19], WH, BKG);
	textout_ex(keycmd, font, "turns to Deadline Window", BEG_L, pos[20], WH, BKG);

	textout_ex(keycmd, font, "Press ESC:", BEG_L, pos[22], WH, BKG);
	textout_ex(keycmd, font, "exits the program", BEG_L, pos[23], WH, BKG);

	textout_ex(keycmd, font, "Types of trail;", BEG_L, pos[25], WH, BKG);
	textout_ex(keycmd, font, "left:  proportional", BEG_L, pos[26], WH, BKG);
	textout_ex(keycmd, font, "to the linear velocity", BEG_L, pos[27], WH, BKG);
	textout_ex(keycmd, font, "right: proportional", BEG_L, pos[28], WH, BKG);
	textout_ex(keycmd, font, "to the Lagrangian value", BEG_L, pos[29], WH, BKG);

	textout_ex(keycmd, font, "page 1/4", NUM_P, pos[N_LINE - 1], WH, BKG);

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
int 	i;				// for cycle index, length in pixel of the column
int 	a;				// for cycle index, half width in pixel of the column
float 	shd[SH];		// intensity value of the light shade

	for (i = 0; i < 255; i++) {

		for (a = 0; a < SH; a++) {
			// shade value evaluated as a normalized gaussian function
			shd[a] = norm_gauss( VSCALE * (i+1), DIST_G * a);
			// shade value adapted to the 255 dimension of the  colour table
			shd[a] *= 255;		

			putpixel(keycmd, NUM_P + a, pos[31] + i, table.data[(int)shd[a]][WH]);
			putpixel(keycmd, NUM_P - a, pos[31] + i, table.data[(int)shd[a]][WH]);
		}
	}
}

// -----------------------------------------------------------------------------
// LAGR_COLUMN: draws a column in the command page to show the shade of colours
// of the lagrangian trail adopted in the program
// -----------------------------------------------------------------------------
void lagr_column()
{
int 	i;				// for cycle index, length in pixel of the column
int 	a;				// for cycle index, half width in pixel of the column
float 	shd[SH];		// shade value of the colour
int 	col;			// value associated to the Lagrangian colour
int 	r, g, b;		// red, green, blue component of the colour

	for (a = 0; a < SH; a++) {
		// value between [0, 1]
		shd[a] = pow(1.5, -a);
		// shade value scaled to the 255 dimension of the colour table
		shd[a] *= 255;
	}

	for (i = 0; i < 255; i++) {

		hsv_to_rgb( i, 1, 1,   &r, &g, &b);		
		col = makecol(r, g, b);

		for (a = 0; a < SH; a++) {
			putpixel(keycmd, (NUM_P * 2) + a, pos[31] + i, table.data[(int)(shd[a])][col]);
			putpixel(keycmd, (NUM_P * 2) - a, pos[31] + i, table.data[(int)(shd[a])][col]);
		}
	}
}

// ------------------------------------------------------------------------------
// INITDATA_PAGE: writes the initial position and velcity of each double pendulum 
// ------------------------------------------------------------------------------
void initdata_page()
{	
int 	i;			// for cycle index
char 	d[30];		// stores the line to be printed on the page

	textout_ex(init_data, font, "Initial Data", TIT_P, pos[0], WH, BKG);

	for (i = 0; i < npend; i++) {

		sprintf(d, "Pendulum %d", i+1);
		textout_ex(init_data, font, d, BEG_L, pos[2 + (6 * i)], WH, BKG);

		sprintf(d, "   Theta1:  %g°", rad2deg(pnd[i].tht1));
		textout_ex(init_data, font, d, BEG_L, pos[3 + (6 * i)], WH, BKG);

		sprintf(d, "   Omega1:  %g [rad/s]", pnd[i].tht1_dot);
		textout_ex(init_data, font, d, BEG_L, pos[4 + (6 * i)], WH, BKG);

		sprintf(d, "   Theta2:  %g°", rad2deg(pnd[i].tht2));
		textout_ex(init_data, font, d, BEG_L, pos[5 + (6 * i)], WH, BKG);

		sprintf(d, "   Omega2:  %g [rad/s]", pnd[i].tht2_dot );
		textout_ex(init_data, font, d, BEG_L, pos[6 + (6 * i)], WH, BKG);
	}

	textout_ex(init_data, font, "page 2/4", NUM_P, pos[N_LINE - 1], WH, BKG);
}

// -----------------------------------------------------------------------------
// GEOMDATA_PAGE: writes the geometrical data of each double pendulum 
// as length and mass
// -----------------------------------------------------------------------------
void geomdata_page()
{
int 	i;			// for cycle index
char 	d[30];		// stores the line to be printed on the page

	textout_ex(geom_data, font, "Geometrical Data", TIT_P, pos[0], WH, BKG);

	for (i = 0; i < npend; i++) {

		sprintf(d, "Pendulum %d", i+1);
		textout_ex(geom_data, font, d, BEG_L, pos[2 + (6 * i)], WH, BKG);

		sprintf(d, "   Length1: %g", pnd[i].l1 );
		textout_ex(geom_data, font, d, BEG_L, pos[3 + (6 * i)], WH, BKG);

		sprintf(d, "   Mass1:   %g", pnd[i].m1 );
		textout_ex(geom_data, font, d, BEG_L, pos[4 + (6 * i)], WH, BKG);

		sprintf(d, "   Length2: %g", pnd[i].l2 );
		textout_ex(geom_data, font, d, BEG_L, pos[5 + (6 * i)], WH, BKG);

		sprintf(d, "   Mass2:   %g", pnd[i].m2 );
		textout_ex(geom_data, font, d, BEG_L, pos[6 + (6 * i)], WH, BKG);

	}

	textout_ex(geom_data, font, "page 3/4", NUM_P, pos[N_LINE - 1], WH, BKG);
}

// -----------------------------------------------------------------------------
// DEADLINE_PAGE: writes the deadline missed of every running thread 
// -----------------------------------------------------------------------------
void deadline_page()
{
	dline_flag = 0;

	textout_ex(dline, font, "Deadline Window", TIT_P, pos[0], WH, BKG);

	// the lines that write the missed deadline must be updated by a thread
	task_create(npend, dline_wintask, WPER, WDL, WPRI);

	textout_ex(dline, font, "page 4/4", NUM_P, pos[N_LINE - 1], WH, BKG);
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
int 	n;			// for cycle index

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
double deg2rad(double deg)
{			
double 	rad;

	rad = (deg / 180.0) * PI;
	return rad;
}

// -----------------------------------------------------------------------------
// RAD2DEG: converts values from radiant to degree
// -----------------------------------------------------------------------------
float rad2deg(double rad)
{
float 	deg;

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
float 	s_tht1, c_tht1;		// sine and cosine of theta_1
float 	s_tht2, c_tht2;		// sine and cosine of theta_2

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
int 	pixel;		// value of the graphic display
float 	l12;		// sum of the length of the two pendulums
	
	l12 = pnd[i].l1 + pnd[i].l2;
	pixel = (int) ((value / l12) * (float)((l / 2) - R));
	return pixel;
}

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
int 	i;				// for cycle index
char 	d[30];			// stores the line to be printed on the page

	sprintf(d, "Window task: %d", tp[npend].dmiss);
	textout_ex(dline, font, d, BEG_L, pos[2], WH, BKG);

	for(i = 0; i < npend; i++) {
			sprintf(d, "Pendulum %d: %d", i + 1, tp[i].dmiss);
			textout_ex(dline, font, d, BEG_L, pos[4 + (2 * i)], WH, BKG);
	}

	if (dline_flag) {
		blit(dline, screen, 0, 0, LWIN, UWIN, dline->w, dline->h);
	}

}

// -----------------------------------------------------------------------------
// STORE_TRAIL:
// stores the new position of the second pendulum in the circular buffer
// -----------------------------------------------------------------------------
void store_trail(int i)
{
int 	k;		// value of the last element of the circular buffer

	k = trail[i].top;
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
int 	k;			// value of the last element of the circular buffer
int 	k_old;		// value of the second to last element of the circular buffer
int 	x[2], y[2];	// graphics coordinates of the last two points

	k = trail[i].top;
	k_old = (k + 1) % TLEN;
	x[0] = trail[i].x[k_old];
	y[0] = trail[i].y[k_old];
	x[1] = trail[i].x[k];
	y[1] = trail[i].y[k];

	light_trail(x[1], y[1], x[0], y[0], i);
	lagr_trail(x[1], y[1], x[0], y[0], i);

	if (change_trail) {
		blit(lagr_tr[i], bkg_pend[i], 0, 0, 0, 0, lagr_tr[i]->w, lagr_tr[i]->h);
	}
	else {
		blit(light_tr[i], bkg_pend[i], 0, 0, 0, 0, light_tr[i]->w, light_tr[i]->h);	
	}
}

// -----------------------------------------------------------------------------
// LIGHT_TRAIL: manages the representation of the Light trail
// -----------------------------------------------------------------------------
void light_trail(int x2, int y2, int x1, int y1, int i)
{
int 	n;				// for cycle index
float 	shd[SH];		// intensity value of the light shade
float 	v;				// instantaneous value of the linear velocity
float 	alpha;			// instantaneous tangent angle of the trajectory
int 	n_s_a;
int 	n_c_a;
int 	nn_c_a;

	v = lin_v(i);

	alpha = tan_tr(i);

	for (n = 0; n < SH; n++) {
		// shade value evaluated as a normalized gaussian function
		shd[n] = norm_gauss( fabs(1 / ((v) + 1)), (0.06 * n));
		// shade value scaled to the 255 dimension of the colour table
		shd[n] *= 255;
	}

	// draws the lines from the darkest to the lightest
	for(n = 4; n >= 0; n--) {
		n_s_a = (int) (n * sin(alpha));
		n_c_a = (int) (n * cos(alpha));
		nn_c_a = (int) ((n - 1) * cos(alpha));

		line(light_tr[i], x2 + n_s_a, y2 + n_c_a,  x1 + n_s_a, y1 + n_c_a, table.data[(int)(shd[n])][WH]);
		line(light_tr[i], x2 - n_s_a, y2 - n_c_a,  x1 - n_s_a, y1 - n_c_a, table.data[(int)(shd[n])][WH]);

		line(light_tr[i], x2 + n_s_a, y2 + nn_c_a,  x1 + n_s_a, y1 + nn_c_a, table.data[(int)(shd[n])][WH]);
		line(light_tr[i], x2 - n_s_a, y2 - nn_c_a,  x1 - n_s_a, y1 - nn_c_a, table.data[(int)(shd[n])][WH]);		
	}
	line(light_tr[i], x2, y2,  x1, y1, table.data[255][WH]);
}

// -----------------------------------------------------------------------------
// LAGR_TRAIL: manages the representation of the Lagrangian trail
// -----------------------------------------------------------------------------
void lagr_trail(int x2, int y2, int x1, int y1, int i)
{
int 	n;				// for cycle index
float 	shd[SH];		// intensity value of the light shade
float 	alpha;			// instantaneous tangent angle of the trajectory
int 	col;			// value associated to the Lagrangian colour
int 	r, g, b;		// red, green, blue component of the colour
int 	n_s_a;
int 	n_c_a;
int 	nn_c_a;

	alpha = tan_tr(i);
	
	// makes the color in function of the Lagrangian value normalized to 255
	hsv_to_rgb( (int)(L(i)) % 255, 1, 1,   &r, &g, &b);		
	col = makecol(r, g, b);

	// generates the shade vector as an exponential function
	for (n = 0; n < SH; n++) {
		shd[n] = pow(1.5, -n);
		shd[n] *= 255;	// shade value scaled to the 255 dimension of the colour tables
	}

	// draws the lines from the darkest to the lightest
	for(n = 4; n >= 0; n--) {	
		n_s_a = (int) (n * sin(alpha));
		n_c_a = (int) (n * cos(alpha));
		nn_c_a = (int) ((n - 1) * cos(alpha));

		line(lagr_tr[i], x2 + n_s_a, y2 + n_c_a,  x1 + n_s_a, y1 + n_c_a, table.data[(int)(shd[n])][col]);
		line(lagr_tr[i], x2 - n_s_a, y2 - n_c_a,  x1 - n_s_a, y1 - n_c_a, table.data[(int)(shd[n])][col]);

		line(lagr_tr[i], x2 + n_s_a, y2 + nn_c_a,  x1 + n_s_a, y1 + nn_c_a, table.data[(int)(shd[n])][col]);
		line(lagr_tr[i], x2 - n_s_a, y2 - nn_c_a,  x1 - n_s_a, y1 - nn_c_a, table.data[(int)(shd[n])][col]);		
	}
	line(lagr_tr[i], x2, y2,  x1, y1, table.data[255][col]);	
}

// -----------------------------------------------------------------------------
// TAN_TR: 
// computes the istant tangent angle of the trajectory of the second pendulum
// -----------------------------------------------------------------------------
float tan_tr(int i)
{
float 	alpha;					// istantaneous tangent angle of the trajectory 
float 	tht1, tht2;				
float 	tht1_dot, tht2_dot;		
float 	l1, l2;
float 	c_th1, c_th2;			// cosine of theta 1 and theta 2
float 	s_th1, s_th2;			// sine of theta 1 and theta 2
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
// computes the value of a Gaussian function, with mean equal to zero and the 
// maximum normalized to one, so that the return value is always between [0, 1].
// Variables:
// 	"var" stores the variance 
// 	"x" stores the value at which the Gaussian function is evaluated
// -----------------------------------------------------------------------------
float norm_gauss(float var, float x)
{
float 	e;		// exponent of the gaussian function

	e = -(x * x / (2 * var * var));

	return exp(e);
}

// -----------------------------------------------------------------------------
// LIN_V: evaluates the istantaneous linear velocity of the second mass attached
// at the end of the second pendulum
// -----------------------------------------------------------------------------
float lin_v (int i)		
{
float 	vel;					// istantaneous linear velocity
float 	tht1, tht2;
float 	tht1_dot, tht2_dot;
float 	l1, l2;
float 	c_th1th2;				// cos(theta_1 - theta2)
float 	N1, N2, N3;				// addends of the velocity expression

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
float 	vel;				// istantaneous linear velocity
float 	l1, l2;
float 	m1, m2;
float 	c_tht1, c_tht2;		// cosine of theta 1 and theta 2
float 	T1, T2;				// kinetic energy
float 	V1, V2;				// gravitational potential energy
float	lag;				// value of the Lagrangian

	tht1 		= pnd[i].tht1;
	tht2 		= pnd[i].tht2;
	tht1_dot 	= pnd[i].tht1_dot;
	vel			= lin_v(i);
	l1 			= pnd[i].l1;
	l2 			= pnd[i].l2;
	m1 			= pnd[i].m1;
	m2 			= pnd[i].m2;
	c_tht1 		= cos(tht1);
	c_tht2 		= cos(tht2);

	T1 = 0.5 * m1 * pow((l1 * tht1_dot), 2);
	T2 = 0.5 * m2 * pow(vel, 2);
	V1 = - m1 * G * l1 * c_tht1;
	V2 = - m2 * G * ((l1 * c_tht1) + (l2 * c_tht2));

	lag = T1 + T2 - V1 - V2;
	
	return lag;
}