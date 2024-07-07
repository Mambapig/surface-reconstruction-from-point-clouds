///*
//cc -o fdm_surface11 fdm_surface11.c -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
//This code corresponds to the Partial differential equation with 32 variables, which can
//also be regarded as version 2 of our original paper. In the previous paper, there are 16 
//variables, and we add another 16 variables in this paper to try to make it more powerful 
//regarding reconstructing more complicated 3D shapes.
//
//Zaiping Zhu 27/10/2021
//*/
//
//#include <GL/glut.h>
//#include <math.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <time.h>
//#include "Header.h"
//#define jb  2    //jb=10
//
//#define nu 15
//#define nv 15
//#define nuv 289
//#define nuvm 277
//
//#define nn 3
////#define icurve 10
//
//int ciggj(double**, int, double*);
//// int fitting(int, int, int, int, double, double, double**, double**, double**, double*);
//// int PDEcoor(int, int, int, int, double, double, double*, double**, double**, double**);
//int fitting_new(int, double, double,double,double, double*, double*, double*, double*);
//int PDEcoor_new(int, double, double,double, double, double*, double, double, double*);
//FILE* in1, * out, * ToMaya, * test, * real_points, * realPointsCoor;
//
//#ifndef CALLBACK
//#define CALLBACK
//#endif
//#define pi 3.1415926
//
//inline double f1(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(q4 * v) * cos(q2 * u) * cos(q4 * v);
//}
//inline double f2(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(q4 * v) * cos(q2 * u) * sin(q4 * v);
//}
//inline double f3(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(q4 * v) * sin(q2 * u) * cos(q4 * v);
//}
//inline double f4(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(q4 * v) * sin(q2 * u) * sin(q4 * v);
//}
//inline double f5(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(-q4 * v) * cos(q2 * u) * cos(q4 * v);
//}
//inline double f6(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(-q4 * v) * cos(q2 * u) * sin(q4 * v);
//}
//inline double f7(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(-q4 * v) * sin(q2 * u) * cos(q4 * v);
//}
//inline double f8(double q2, double q4, double u, double v)
//{
//	return exp(q2 * u) * exp(-q4 * v) * sin(q2 * u) * sin(q4 * v);
//}
//inline double f9(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(q4 * v) * cos(q2 * u) * cos(q4 * v);
//}
//inline double f10(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(q4 * v) * cos(q2 * u) * sin(q4 * v);
//}
//inline double f11(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(q4 * v) * sin(q2 * u) * cos(q4 * v);
//}
//inline double f12(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(q4 * v) * sin(q2 * u) * sin(q4 * v);
//}
//inline double f13(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q4 * v) * cos(q2 * u) * cos(q4 * v);
//}
//inline double f14(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q4 * v) * cos(q2 * u) * sin(q4 * v);
//}
//inline double f15(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q4 * v) * sin(q2 * u) * cos(q4 * v);
//}
//inline double f16(double q2, double q4, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q4 * v) * sin(q2 * u) * sin(q4 * v);
//}
//
//inline double f17(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(q6 * v) * cos(q3 * u) * cos(q6 * v);
//}
//inline double f18(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(q6 * v) * cos(q3 * u) * sin(q6 * v);
//}
//inline double f19(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(q6 * v) * sin(q3 * u) * cos(q6 * v);
//}
//inline double f20(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(q6 * v) * sin(q3 * u) * sin(q6 * v);
//}
//inline double f21(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(-q6 * v) * cos(q3 * u) * cos(q6 * v);
//}
//inline double f22(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(-q6 * v) * cos(q3 * u) * sin(q6 * v);
//}
//inline double f23(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(-q6 * v) * sin(q3 * u) * cos(q6 * v);
//}
//inline double f24(double q3, double q6, double u, double v)
//{
//	return exp(q3 * u) * exp(-q6 * v) * sin(q3 * u) * sin(q6 * v);
//}
//inline double f25(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(q6 * v) * cos(q3 * u) * cos(q6 * v);
//}
//inline double f26(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(q6 * v) * cos(q3 * u) * sin(q6 * v);
//}
//inline double f27(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(q6 * v) * sin(q3 * u) * cos(q6 * v);
//}
//inline double f28(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(q6 * v) * sin(q3 * u) * sin(q6 * v);
//}
//inline double f29(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q6 * v) * cos(q3 * u) * cos(q6 * v);
//}
//inline double f30(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q6 * v) * cos(q3 * u) * sin(q6 * v);
//}
//inline double f31(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q6 * v) * sin(q3 * u) * cos(q6 * v);
//}
//inline double f32(double q3, double q6, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q6 * v) * sin(q3 * u) * sin(q6 * v);
//}
//
//GLuint startList;
//
//void CALLBACK errorCallback(GLenum errorCode)
//{
//	const GLubyte* estring;
//
//	estring = gluErrorString(errorCode);
//	fprintf(stderr, "Quadric Error: %s\n", estring);
//	exit(0);
//}
//double dx32[32], dy32[32], dz32[32];
//void init(void)
//{
//	int i, is, j, js, k, n, n0, l, m, iu, iv, jj, iuu, icurve, ie, je, points;
//	double fk, du, dLj, L, ui, vj, mFaceWidth,
//		d32[32],
//		//dx16[16], dy16[16], dz16[16], 
//		q2, q4, aax[32][32], bbx[32], EoneM, EoneA,
//		EtwoM, EtwoA, EoneC, EtwoC, //x_coor[520], y_coor[520], z_coor[520];
//		q3, q6;
//	// X83[][0]=x component, X83[][1]=y component, X83[][2]=zx component; U82[][0]= u component, U82[][1]= v component.
//	//	char char1, char2; uij[100][100], vij[100][100],X100[100][100], 
//
//	//  clock_t start, finish, duration;
//	//GLUquadricObj* qobj;
//	//GLfloat mat_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
//	//GLfloat mat_specular[] = { 0.3, 0.3, 0.3, 1.0 };
//	//GLfloat mat_shininess[] = { 100.0 };
//	//GLfloat light_position[] = { 1000.0, 70000.0, -100000.0, 0.0 };
//	//GLfloat model_ambient[] = { 0.9, 0.9, 0.9, 1.0 };
//	//    glClearColor(0.382, 0.382, 0.382, 1.0);
//
//	//glClearColor(1.0, 1.0, 1.0, 1.0);
//
//	//glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
//	//glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//	//glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
//	//glLightfv(GL_LIGHT0, GL_POSITION, light_position);
//	//glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambient);
//
//	//glEnable(GL_LIGHTING);
//	//glEnable(GL_LIGHT0);
//	//glEnable(GL_DEPTH_TEST);
//
//
//	//startList = glGenLists(1);
//	//qobj = gluNewQuadric();
//	//gluQuadricDrawStyle(qobj, GLU_FILL);
//	//gluQuadricNormals(qobj, GLU_SMOOTH);
//	//glNewList(startList, GL_COMPILE);
//	double* x_coor = new double[3000];
//	double* y_coor = new double[3000];
//	double* z_coor = new double[3000];
//	double* x_coor_real = new double[3000];
//	double* y_coor_real = new double[3000];
//	double* z_coor_real = new double[3000];
//	//	n=2*jb+1;
//	n = 2 * jb - 1;
//	n0 = 2 * jb;
//	q2 = 0.1;
//	q4 = 0.15;
//
//	q3 = 0.2;
//	q6 = 0.3;
//	//in1 = fopen("In1.txt", "rt");
//	//out = fopen("out.txt", "wt");
//	out = fopen("CylinderWOLid1_v2/cylinderwolid1_output_v2.txt", "wt");
//	//ToMaya = fopen("curve.mel", "wt");
//	ToMaya = fopen("CylinderWOLid1_v2/cylinderwolid1_curve_v2.mel", "wt");
//	realPointsCoor = fopen("CylinderWOLid1_v2/cylinderwolid1realPointsCoorCal_v2.txt", "wt");
//	//test = fopen("test.txt", "rt");
//	// points after projection
//	test = fopen("CylinderWOLid1_v2/cylinderwolid1_LowPoly_v2.txt", "rt");
//
//	// Points before projection
//	real_points = fopen("CylinderWOLid1_v2/cylinderwolid1_HighPoly_v2.txt", "rt");
//	fscanf(test, "%d\n", &points);
//	printf("points= %d\n", points);
//	for (int i = 0; i <= points - 1; i++)
//	{
//		fscanf(test, "%le %le\n", &x_coor[i], &y_coor[i]);
//		printf("points[%d] is:  %le %le\n", i, x_coor[i], y_coor[i]);
//
//		fscanf(real_points, "%le %le %le\n", &x_coor_real[i], &y_coor_real[i], &z_coor_real[i]);
//		printf("real_points[%d] is:  %le %le %le\n", i, x_coor_real[i], y_coor_real[i], z_coor_real[i]);
//	}
//	double x_coor_max = x_coor[0];
//	double x_coor_min = x_coor[0];
//	double y_coor_max = y_coor[0];
//	double y_coor_min = y_coor[0];
//	double z_coor_max = z_coor[0];
//	double z_coor_min = z_coor[0];
//	double u_length = 0.0;
//	double v_length = 0.0;
//	double* Ui = new double[3000];
//	double* Vi = new double[3000];
//	double* temp = new double[3000];
//	double* x_cal = new double[3000];
//	double* y_cal = new double[3000];
//	double* z_cal = new double[3000];
//	/*for (int i = 0; i <= points - 1; i++)
//	{
//		Ui[i] = x_coor[i];
//		Vi[i] = y_coor[i];
//		fprintf(out, "i,j,ui= %d %le\n", i, Ui[i]);
//		fprintf(out, "i,j,vi= %d %le\n", i, Vi[i]);
//	}*/
//
//	for (int i = 0; i <= points - 1; i++)
//	{
//
//		if (x_coor[i] > x_coor_max)
//			x_coor_max = x_coor[i];
//		if (x_coor[i] < x_coor_min)
//			x_coor_min = x_coor[i];
//		if (y_coor[i] > y_coor_max)
//			y_coor_max = y_coor[i];
//		if (y_coor[i] < y_coor_min)
//			y_coor_min = y_coor[i];
//		if (z_coor[i] > z_coor_max)
//			z_coor_max = z_coor[i];
//		if (z_coor[i] < z_coor_min)
//			z_coor_min = z_coor[i];
//	}
//	u_length = x_coor_max - x_coor_min;
//	v_length = y_coor_max - y_coor_min;
//	for (int i = 0; i <= points - 1; i++)
//	{
//		Ui[i] = (x_coor[i] - x_coor_min) / u_length;
//		Vi[i] = y_coor[i];
//		//Vi[i] = (y_coor[i] - y_coor_min) / v_length;
//		fprintf(out, "i,j,ui= %d %le\n", i, Ui[i]);
//		fprintf(out, "i,j,vi= %d %le\n", i, Vi[i]);
//	}
//	//  x component
//	for (i = 0; i < points; i++)
//		temp[i] = x_coor_real[i];
//	fitting_new(points, q2, q4,q3, q6, Ui, Vi, temp, d32);
//	for (i = 0; i < 32; i++)
//		dx32[i] = d32[i];
//	//  y component
//	for (i = 0; i < points; i++)
//		temp[i] = y_coor_real[i];
//	fitting_new(points, q2, q4,q3, q6, Ui, Vi, temp, d32);
//	for (i = 0; i < 32; i++)
//		dy32[i] = d32[i];
//	//  z component
//	for (i = 0; i < points; i++)
//		temp[i] = z_coor_real[i];
//	fitting_new(points, q2, q4,q3, q6, Ui, Vi, temp, d32);
//	for (i = 0; i < 32; i++)
//		dz32[i] = d32[i];
//
//	// Calculate the coordinate values of PDE surfaces
//	//PDEcoor_new(points, q2, q4, dx16, Ui, Vi, temp);
//	//for (i = 0; i <= points-1; i++)
//	//	x_cal[i] = temp[i];
//	//PDEcoor_new(points, q2, q4, dy16, Ui, Vi, temp);
//	//for (i = 0; i <= points - 1; i++)
//	//	y_cal[i] = temp[i];
//	//PDEcoor_new(points, q2, q4, dz16, Ui, Vi, temp);
//	//for (i = 0; i <= points - 1; i++)
//	//	z_cal[i] = temp[i];
//
//	int I = 10;
//	int J = 50;
//	double delta_u =1.0 / I;
//	double delta_v = 1.0 / J;
//	for (i = 0; i <= I; i++)
//	{
//		double u_i = i * delta_u;
//		fprintf(ToMaya, "curve -d 3\n");
//		for (j = 0; j <= J; j++)
//		{
//			double v_j = j * delta_v;
//			PDEcoor_new(points, q2, q4,q3, q6, dx32, u_i, v_j, x_cal);
//			PDEcoor_new(points, q2, q4,q3, q6, dy32, u_i, v_j, y_cal);
//			PDEcoor_new(points, q2, q4,q3, q6, dz32, u_i, v_j, z_cal);
//			fprintf(ToMaya, "-p   %le %le %le\n", *x_cal, *y_cal, *z_cal);
//			//fprintf(realPointsCoor, "%le %le %le\n", *x_cal, *y_cal, *z_cal);
//		}
//		fprintf(ToMaya, ";\n");
//	}
//	for (i = 0; i < points; i++)
//	{
//		PDEcoor_new(points, q2, q4, q3, q6,dx32, Ui[i], Vi[i], x_cal);
//		PDEcoor_new(points, q2, q4, q3, q6,dy32, Ui[i], Vi[i], y_cal);
//		PDEcoor_new(points, q2, q4, q3, q6,dz32, Ui[i], Vi[i], z_cal);
//		fprintf(realPointsCoor, "%le %le %le\n", double(*x_cal), double(*y_cal), double(*z_cal));
//	}
//	mFaceWidth = 0.0;
//
//
//	delete[] Ui;
//	delete[] Vi;
//	delete[] temp;
//	delete[] x_cal;
//	delete[] y_cal;
//	delete[] z_cal;
//
//	delete[] x_coor;
//	delete[] y_coor;
//	delete[] z_coor;
//	delete[] x_coor_real;
//	delete[] y_coor_real;
//	delete[] z_coor_real;
//	fclose(out);
//	fclose(test);
//	fclose(ToMaya);
//	fclose(realPointsCoor);
//	//	printf("Pass 0\n");
//
//	//glEndList();
//
//}
//
//void display(void)
//{
//
//	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//	glPushMatrix();
//
//	glEnable(GL_COLOR_MATERIAL);
//	glEnable(GL_LIGHTING);
//	glShadeModel(GL_SMOOTH);
//
//	glColor3f(0.5, 0.0, 1.0);
//	glTranslatef(0.0, 0.0, 0.0);
//	glPushMatrix();
//	glScalef(0.15, 0.15, 0.15);
//	//   glRotatef(45.0,45.0, 45.0, 1.0);
//	glRotatef(90.0, 90.0, 0.0, 1.0);
//
//	glCallList(startList);
//	glPopMatrix();
//
//	glFlush();
//}
//
//void reshape(int w, int h)
//{
//	glViewport(0, 0, (GLsizei)w, (GLsizei)h);
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	if (w <= h)
//		glOrtho(-2.5, 2.5, -2.5 * (GLfloat)h / (GLfloat)w,
//			2.5 * (GLfloat)h / (GLfloat)w, -10.0, 10.0);
//	else
//		glOrtho(-2.5 * (GLfloat)w / (GLfloat)h,
//			2.5 * (GLfloat)w / (GLfloat)h, -2.5, 2.5, -10.0, 10.0);
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//}
//
//void keyboard(unsigned char key, int x, int y)
//{
//	switch (key) {
//	case 27:
//		exit(0);
//		break;
//	}
//}
//
//int main(int argc, char** argv)
//{
//	glutInit(&argc, argv);
//	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
//	glutInitWindowSize(500, 500);
//	glutInitWindowPosition(100, 100);
//	glutCreateWindow(argv[0]);
//	init();
//	glutDisplayFunc(display);
//	glutReshapeFunc(reshape);
//	glutKeyboardFunc(keyboard);
//	glutMainLoop();
//	return 0;
//}
//
//int ciggj(double** a, int n, double* b)
//{
//	int i, j, k, is, u, v;
//	int* js = new int[n];
//	double d, t;
//	for (k = 0; k <= n - 1; k++)
//	{
//		d = 0.0;
//		for (i = k; i <= n - 1; i++)
//			for (j = k; j <= n - 1; j++)
//			{
//				t = fabs(a[i][j]);
//				if (t > d) { d = t; js[k] = j; is = i; }
//			}
//		if (d + 1.0 == 1.0)
//		{
//			delete[] js; printf("fail\n"); return(0);
//		}
//		if (is != k)
//		{
//			for (j = k; j <= n - 1; j++)
//			{
//				u = k * n + j; v = is * n + j;
//				t = a[k][j]; a[k][j] = a[is][j]; a[is][j] = t;
//			}
//			t = b[k]; b[k] = b[is]; b[is] = t;
//		}
//		if (js[k] != k)
//			for (i = 0; i <= n - 1; i++)
//			{
//				u = i * n + k; v = i * n + js[k];
//				t = a[i][k]; a[i][k] = a[i][js[k]]; a[i][js[k]] = t;
//			}
//		t = a[k][k];
//		for (j = k + 1; j <= n - 1; j++)
//		{
//			u = k * n + j;
//			if (a[k][j] != 0.0)a[k][j] = a[k][j] / t;
//		}
//		b[k] = b[k] / t;
//		for (j = k + 1; j <= n - 1; j++)
//		{
//			u = k * n + j;
//			if (a[k][j] != 0.0)
//			{
//				for (i = 0; i <= n - 1; i++)
//				{
//					v = i * n + k;
//					if ((i != k) && (a[i][k] != 0.0))
//					{
//						is = i * n + j;
//						a[i][j] = a[i][j] - a[i][k] * a[k][j];
//					}
//				}
//			}
//		}
//		for (i = 0; i <= n - 1; i++)
//		{
//			u = i * n + k;
//			if ((i != k) && (a[i][k] != 0.0))
//				b[i] = b[i] - a[i][k] * b[k];
//		}
//	}
//	for (k = n - 1; k >= 0; k--)
//		if (k != js[k])
//		{
//			t = b[k]; b[k] = b[js[k]]; b[js[k]] = t;
//		}
//	delete[] js;
//	return(1);
//}
//int fitting_new(int points, double q2, double q4,double q3, double q6, double* Ui, double* Vi, double* temp, double* d32)
//{
//	int i, j, k;
//	double ui, vj, fk;
//	double* bb = new double[32]; // 
//	double* bbx = new double[32]; // 
//	double** aa = new double* [32]; // 
//	double** aax = new double* [32]; // 
//	for (i = 0; i < 32; i++)
//	{
//		aa[i] = new double[32];
//		aax[i] = new double[32];
//	}
//	//	printf("is,ie,js,je,q2,q4= %d %d %d %d %le %le\n", is, ie, js, je, q2, q4);
//	for (i = 0; i < 32; i++)
//	{
//		bb[i] = 0.0;
//		for (j = 0; j < 32; j++)
//		{
//			aa[i][j] = 0.0;
//		}
//	}
//	for (k = 1; k <= 32; k++)
//	{
//		for (i = 0; i <= points - 1; i++)
//		{
//			if (k == 1)fk = f1(q2, q4, Ui[i], Vi[i]);
//			if (k == 2)fk = f2(q2, q4, Ui[i], Vi[i]);
//			if (k == 3)fk = f3(q2, q4, Ui[i], Vi[i]);
//			if (k == 4)fk = f4(q2, q4, Ui[i], Vi[i]);
//			if (k == 5)fk = f5(q2, q4, Ui[i], Vi[i]);
//			if (k == 6)fk = f6(q2, q4, Ui[i], Vi[i]);
//			if (k == 7)fk = f7(q2, q4, Ui[i], Vi[i]);
//			if (k == 8)fk = f8(q2, q4, Ui[i], Vi[i]);
//			if (k == 9)fk = f9(q2, q4, Ui[i], Vi[i]);
//			if (k == 10)fk = f10(q2, q4, Ui[i], Vi[i]);
//			if (k == 11)fk = f11(q2, q4, Ui[i], Vi[i]);
//			if (k == 12)fk = f12(q2, q4, Ui[i], Vi[i]);
//			if (k == 13)fk = f13(q2, q4, Ui[i], Vi[i]);
//			if (k == 14)fk = f14(q2, q4, Ui[i], Vi[i]);
//			if (k == 15)fk = f15(q2, q4, Ui[i], Vi[i]);
//			if (k == 16)fk = f16(q2, q4, Ui[i], Vi[i]);
//			if (k == 17)fk = f17(q3, q6, Ui[i], Vi[i]);
//			if (k == 18)fk = f18(q3, q6, Ui[i], Vi[i]);
//			if (k == 19)fk = f19(q3, q6, Ui[i], Vi[i]);
//			if (k == 20)fk = f20(q3, q6, Ui[i], Vi[i]);
//			if (k == 21)fk = f21(q3, q6, Ui[i], Vi[i]);
//			if (k ==22)fk = f22(q3, q6, Ui[i], Vi[i]);
//			if (k == 23)fk = f23(q3, q6, Ui[i], Vi[i]);
//			if (k == 24)fk = f24(q3, q6, Ui[i], Vi[i]);
//			if (k == 25)fk = f25(q3, q6, Ui[i], Vi[i]);
//			if (k == 26)fk = f26(q3, q6, Ui[i], Vi[i]);
//			if (k == 27)fk = f27(q3, q6, Ui[i], Vi[i]);
//			if (k == 28)fk = f28(q3, q6, Ui[i], Vi[i]);
//			if (k == 29)fk = f29(q3, q6, Ui[i], Vi[i]);
//			if (k == 30)fk = f30(q3, q6, Ui[i], Vi[i]);
//			if (k == 31)fk = f31(q3, q6, Ui[i], Vi[i]);
//			if (k == 32)fk = f32(q3, q6, Ui[i], Vi[i]);
//			//				printf("k,ui,vj,fk= %d %le %le %le\n", k, ui, vj, fk);
//			aa[k - 1][0] = aa[k - 1][0] + f1(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][1] = aa[k - 1][1] + f2(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][2] = aa[k - 1][2] + f3(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][3] = aa[k - 1][3] + f4(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][4] = aa[k - 1][4] + f5(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][5] = aa[k - 1][5] + f6(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][6] = aa[k - 1][6] + f7(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][7] = aa[k - 1][7] + f8(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][8] = aa[k - 1][8] + f9(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][9] = aa[k - 1][9] + f10(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][10] = aa[k - 1][10] + f11(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][11] = aa[k - 1][11] + f12(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][12] = aa[k - 1][12] + f13(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][13] = aa[k - 1][13] + f14(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][14] = aa[k - 1][14] + f15(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][15] = aa[k - 1][15] + f16(q2, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][16] = aa[k - 1][16] + f17(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][17] = aa[k - 1][17] + f18(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][18] = aa[k - 1][18] + f19(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][19] = aa[k - 1][19] + f20(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][20] = aa[k - 1][20] + f21(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][21] = aa[k - 1][21] + f22(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][22] = aa[k - 1][22] + f23(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][23] = aa[k - 1][23] + f24(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][24] = aa[k - 1][24] + f25(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][25] = aa[k - 1][25] + f26(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][26] = aa[k - 1][26] + f27(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][27] = aa[k - 1][27] + f28(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][28] = aa[k - 1][28] + f29(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][29] = aa[k - 1][29] + f30(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][30] = aa[k - 1][30] + f31(q3, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][31] = aa[k - 1][31] + f32(q3, q6, Ui[i], Vi[i]) * fk;
//			bb[k - 1] = bb[k - 1] + temp[i] * fk;
//		}
//		//		printf("k,bb[k-1],aa[k-1][0-1]= %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
//		//			k, bb[k - 1], aa[k - 1][0], aa[k - 1][1], aa[k - 1][2], aa[k - 1][3], aa[k - 1][4], aa[k - 1][5], aa[k - 1][6],
//		//			aa[k - 1][7], aa[k - 1][8], aa[k - 1][9], aa[k - 1][10], aa[k - 1][11], aa[k - 1][12], aa[k - 1][13], aa[k - 1][14],
//		//			aa[k - 1][15]);
//	}
//	for (i = 0; i < 32; i++)
//	{
//		bbx[i] = bb[i];
//		//		fprintf(out, "i,bb[i]= %d %le\n", i, bb[i]);
//		for (j = 0; j < 32; j++)
//		{
//			aax[i][j] = aa[i][j];
//			//			fprintf(out, "i,j,aa[i][j]= %d %d %le\n", i, j, aa[i][j]);
//		}
//	}
//	//	n = 16;
//	ciggj(aa, 32, bb);
//	for (i = 0; i < 32; i++)
//	{
//		d32[i] = bb[i];
//	}
//	return(1);
//	delete[] bb;
//	delete[] bbx;
//	for (i = 0; i < 32; i++)
//	{
//		delete aa[i];
//		delete aax[i];
//	}
//	delete[] aa;
//	delete[] aax;
//}
//int PDEcoor_new(int points, double q2, double q4, double q3, double q6, double* d32, double Ui, double Vi, double* temp)
//{
//	int i, j;
//	//	printf("is,ie,js,je,q2,q4= %d %d %d %d %le %le\n", is, ie, js, je, q2, q4);
//	for (int i = 0; i < points; i++)
//	{
//		temp[i] = d32[0] * f1(q2, q4, Ui, Vi) + d32[1] * f2(q2, q4, Ui, Vi) + d32[2] * f3(q2, q4, Ui, Vi) + d32[3]
//			* f4(q2, q4, Ui, Vi) + d32[4] * f5(q2, q4, Ui, Vi) + d32[5] * f6(q2, q4, Ui, Vi) + d32[6] * f7(q2, q4, Ui, Vi)
//			+ d32[7] * f8(q2, q4, Ui, Vi) + d32[8] * f9(q2, q4, Ui, Vi) + d32[9] * f10(q2, q4, Ui, Vi) + d32[10]
//			* f11(q2, q4, Ui, Vi) + d32[11] * f12(q2, q4, Ui, Vi) + d32[12] * f13(q2, q4, Ui, Vi) + d32[13]
//			* f14(q2, q4, Ui, Vi) + d32[14] * f15(q2, q4, Ui, Vi) + d32[15] * f16(q2, q4, Ui, Vi) + d32[16] * f17(q3, q6, Ui,Vi)
//			+ d32[17] * f18(q3, q6, Ui, Vi) + d32[18] * f19(q3, q6, Ui, Vi) + d32[19] * f20(q3, q6, Ui, Vi)+
//			d32[20] * f21(q3, q6, Ui, Vi) + d32[21] * f22(q3, q6, Ui, Vi) + d32[22] * f23(q3, q6, Ui, Vi) + 
//			d32[23] * f24(q3, q6, Ui, Vi) + d32[24] * f25(q3, q6, Ui, Vi) + d32[25] * f26(q3, q6, Ui, Vi)+
//			d32[26] * f27(q3, q6, Ui, Vi) + d32[27] * f28(q3, q6, Ui, Vi) + d32[28] * f29(q3, q6, Ui, Vi)+
//			d32[29] * f30(q3, q6, Ui, Vi) + d32[30] * f31(q3, q6, Ui, Vi) + d32[31] * f32(q3, q6, Ui, Vi);
//	}
//	return(1);
//}
//
//
//
//
