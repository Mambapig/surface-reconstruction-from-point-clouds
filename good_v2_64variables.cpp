///*
//cc -o fdm_surface11 fdm_surface11.c -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
//This code corresponds to the Partial differential equation with 64 variables, which can
//also be regarded as version 3 of our original paper. In the previous paper, there are 16 
//variables, and we add another 48 variables in this paper to try to make it more powerful 
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
//inline double f33(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(q6 * v) * cos(q2 * u) * cos(q6 * v);
//}
//inline double f34(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(q6 * v) * cos(q2 * u) * sin(q6 * v);
//}
//inline double f35(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(q6 * v) * sin(q2 * u) * cos(q6 * v);
//}
//inline double f36(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(q6 * v) * sin(q2 * u) * sin(q6 * v);
//}
//inline double f37(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(-q6 * v) * cos(q2 * u) * cos(q6 * v);
//}
//inline double f38(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(-q6 * v) * cos(q2 * u) * sin(q6 * v);
//}
//inline double f39(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(-q6 * v) * sin(q2 * u) * cos(q6 * v);
//}
//inline double f40(double q2, double q6, double u, double v)
//{
//	return exp(q2 * u) * exp(-q6 * v) * sin(q2 * u) * sin(q6 * v);
//}
//inline double f41(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(q6 * v) * cos(q2 * u) * cos(q6 * v);
//}
//inline double f42(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(q6 * v) * cos(q2 * u) * sin(q6 * v);
//}
//inline double f43(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(q6 * v) * sin(q2 * u) * cos(q6 * v);
//}
//inline double f44(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(q6 * v) * sin(q2 * u) * sin(q6 * v);
//}
//inline double f45(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q6 * v) * cos(q2 * u) * cos(q6 * v);
//}
//inline double f46(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q6 * v) * cos(q2 * u) * sin(q6 * v);
//}
//inline double f47(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q6 * v) * sin(q2 * u) * cos(q6 * v);
//}
//inline double f48(double q2, double q6, double u, double v)
//{
//	return exp(-q2 * u) * exp(-q6 * v) * sin(q2 * u) * sin(q6 * v);
//}
//
//inline double f49(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(q4 * v) * cos(q3 * u) * cos(q4 * v);
//}
//inline double f50(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(q4 * v) * cos(q3 * u) * sin(q4 * v);
//}
//inline double f51(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(q4 * v) * sin(q3 * u) * cos(q4 * v);
//}
//inline double f52(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(q4 * v) * sin(q3 * u) * sin(q4 * v);
//}
//inline double f53(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(-q4 * v) * cos(q3 * u) * cos(q4 * v);
//}
//inline double f54(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(-q4 * v) * cos(q3 * u) * sin(q4 * v);
//}
//inline double f55(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(-q4 * v) * sin(q3 * u) * cos(q4 * v);
//}
//inline double f56(double q3, double q4, double u, double v)
//{
//	return exp(q3 * u) * exp(-q4 * v) * sin(q3 * u) * sin(q4 * v);
//}
//inline double f57(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(q4 * v) * cos(q3 * u) * cos(q4 * v);
//}
//inline double f58(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(q4 * v) * cos(q3 * u) * sin(q4 * v);
//}
//inline double f59(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(q4 * v) * sin(q3 * u) * cos(q4 * v);
//}
//inline double f60(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(q4 * v) * sin(q3 * u) * sin(q4 * v);
//}
//inline double f61(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q4 * v) * cos(q3 * u) * cos(q4 * v);
//}
//inline double f62(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q4 * v) * cos(q3 * u) * sin(q4 * v);
//}
//inline double f63(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q4 * v) * sin(q3 * u) * cos(q4 * v);
//}
//inline double f64(double q3, double q4, double u, double v)
//{
//	return exp(-q3 * u) * exp(-q4 * v) * sin(q3 * u) * sin(q4 * v);
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
//double dx64[64], dy64[64], dz64[64];
//void init(void)
//{
//	int i, is, j, js, k, n, n0, l, m, iu, iv, jj, iuu, icurve, ie, je, points;
//	double fk, du, dLj, L, ui, vj, mFaceWidth,
//		d64[64],
//		//dx16[16], dy16[16], dz16[16], 
//		q2, q4, aax[64][64], bbx[64], EoneM, EoneA,
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
//	out = fopen("CylinderWOLid1_v3/cylinderwolid1_output_v3.txt", "wt");
//	//ToMaya = fopen("curve.mel", "wt");
//	ToMaya = fopen("CylinderWOLid1_v3/cylinderwolid1_curve_v3.mel", "wt");
//	realPointsCoor = fopen("CylinderWOLid1_v3/cylinderwolid1realPointsCoorCal_v3.txt", "wt");
//	//test = fopen("test.txt", "rt");
//	// points after projection
//	test = fopen("CylinderWOLid1_v3/cylinderwolid1_LowPoly_v3.txt", "rt");
//
//	// Points before projection
//	real_points = fopen("CylinderWOLid1_v3/cylinderwolid1_HighPoly_v3.txt", "rt");
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
//	fitting_new(points, q2, q4,q3, q6, Ui, Vi, temp, d64);
//	for (i = 0; i < 64; i++)
//		dx64[i] = d64[i];
//	//  y component
//	for (i = 0; i < points; i++)
//		temp[i] = y_coor_real[i];
//	fitting_new(points, q2, q4,q3, q6, Ui, Vi, temp, d64);
//	for (i = 0; i < 64; i++)
//		dy64[i] = d64[i];
//	//  z component
//	for (i = 0; i < points; i++)
//		temp[i] = z_coor_real[i];
//	fitting_new(points, q2, q4,q3, q6, Ui, Vi, temp, d64);
//	for (i = 0; i < 64; i++)
//		dz64[i] = d64[i];
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
//			PDEcoor_new(points, q2, q4,q3, q6, dx64, u_i, v_j, x_cal);
//			PDEcoor_new(points, q2, q4,q3, q6, dy64, u_i, v_j, y_cal);
//			PDEcoor_new(points, q2, q4,q3, q6, dz64, u_i, v_j, z_cal);
//			fprintf(ToMaya, "-p   %le %le %le\n", *x_cal, *y_cal, *z_cal);
//			//fprintf(realPointsCoor, "%le %le %le\n", *x_cal, *y_cal, *z_cal);
//		}
//		fprintf(ToMaya, ";\n");
//	}
//	for (i = 0; i < points; i++)
//	{
//		PDEcoor_new(points, q2, q4, q3, q6,dx64, Ui[i], Vi[i], x_cal);
//		PDEcoor_new(points, q2, q4, q3, q6,dy64, Ui[i], Vi[i], y_cal);
//		PDEcoor_new(points, q2, q4, q3, q6,dz64, Ui[i], Vi[i], z_cal);
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
//int fitting_new(int points, double q2, double q4,double q3, double q6, double* Ui, double* Vi, double* temp, double* d64)
//{
//	int i, j, k;
//	double ui, vj, fk;
//	double* bb = new double[64]; // 
//	double* bbx = new double[64]; // 
//	double** aa = new double* [64]; // 
//	double** aax = new double* [64]; // 
//	for (i = 0; i < 64; i++)
//	{
//		aa[i] = new double[64];
//		aax[i] = new double[64];
//	}
//	//	printf("is,ie,js,je,q2,q4= %d %d %d %d %le %le\n", is, ie, js, je, q2, q4);
//	for (i = 0; i < 64; i++)
//	{
//		bb[i] = 0.0;
//		for (j = 0; j < 64; j++)
//		{
//			aa[i][j] = 0.0;
//		}
//	}
//	for (k = 1; k <= 64; k++)
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
//
//			if (k == 33)fk = f33(q2, q6, Ui[i], Vi[i]);
//			if (k == 34)fk = f34(q2, q6, Ui[i], Vi[i]);
//			if (k == 35)fk = f35(q2, q6, Ui[i], Vi[i]);
//			if (k == 46)fk = f36(q2, q6, Ui[i], Vi[i]);
//			if (k == 37)fk = f37(q2, q6, Ui[i], Vi[i]);
//			if (k == 38)fk = f38(q2, q6, Ui[i], Vi[i]);
//			if (k == 39)fk = f39(q2, q6, Ui[i], Vi[i]);
//			if (k == 40)fk = f40(q2, q6, Ui[i], Vi[i]);
//			if (k == 41)fk = f41(q2, q6, Ui[i], Vi[i]);
//			if (k == 42)fk = f42(q2, q6, Ui[i], Vi[i]);
//			if (k == 43)fk = f43(q2, q6, Ui[i], Vi[i]);
//			if (k == 44)fk = f44(q2, q6, Ui[i], Vi[i]);
//			if (k == 45)fk = f45(q2, q6, Ui[i], Vi[i]);
//			if (k == 46)fk = f46(q2, q6, Ui[i], Vi[i]);
//			if (k == 47)fk = f47(q2, q6, Ui[i], Vi[i]);
//			if (k == 48)fk = f48(q2, q6, Ui[i], Vi[i]);
//			if (k == 49)fk = f49(q3, q4, Ui[i], Vi[i]);
//			if (k == 50)fk = f50(q3, q4, Ui[i], Vi[i]);
//			if (k == 51)fk = f51(q3, q4, Ui[i], Vi[i]);
//			if (k == 52)fk = f52(q3, q4, Ui[i], Vi[i]);
//			if (k == 53)fk = f53(q3, q4, Ui[i], Vi[i]);
//			if (k == 54)fk = f54(q3, q4, Ui[i], Vi[i]);
//			if (k == 55)fk = f55(q3, q4, Ui[i], Vi[i]);
//			if (k == 56)fk = f56(q3, q4, Ui[i], Vi[i]);
//			if (k == 57)fk = f57(q3, q4, Ui[i], Vi[i]);
//			if (k == 58)fk = f58(q3, q4, Ui[i], Vi[i]);
//			if (k == 59)fk = f59(q3, q4, Ui[i], Vi[i]);
//			if (k == 60)fk = f60(q3, q4, Ui[i], Vi[i]);
//			if (k == 61)fk = f61(q3, q4, Ui[i], Vi[i]);
//			if (k == 62)fk = f62(q3, q4, Ui[i], Vi[i]);
//			if (k == 63)fk = f63(q3, q4, Ui[i], Vi[i]);
//			if (k == 64)fk = f64(q3, q4, Ui[i], Vi[i]);
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
//
//			aa[k - 1][32] = aa[k - 1][32] + f33(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][33] = aa[k - 1][33] + f34(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][34] = aa[k - 1][34] + f35(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][35] = aa[k - 1][35] + f36(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][36] = aa[k - 1][36] + f37(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][37] = aa[k - 1][37] + f38(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][38] = aa[k - 1][38] + f39(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][39] = aa[k - 1][39] + f40(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][40] = aa[k - 1][40] + f41(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][41] = aa[k - 1][41] + f42(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][42] = aa[k - 1][42] + f43(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][43] = aa[k - 1][43] + f44(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][44] = aa[k - 1][44] + f45(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][45] = aa[k - 1][45] + f46(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][46] = aa[k - 1][46] + f47(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][47] = aa[k - 1][47] + f48(q2, q6, Ui[i], Vi[i]) * fk;
//			aa[k - 1][48] = aa[k - 1][48] + f49(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][49] = aa[k - 1][49] + f50(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][50] = aa[k - 1][50] + f51(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][51] = aa[k - 1][51] + f52(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][52] = aa[k - 1][52] + f53(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][53] = aa[k - 1][53] + f54(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][54] = aa[k - 1][54] + f55(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][55] = aa[k - 1][55] + f56(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][56] = aa[k - 1][56] + f57(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][57] = aa[k - 1][57] + f58(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][58] = aa[k - 1][58] + f59(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][59] = aa[k - 1][59] + f60(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][60] = aa[k - 1][60] + f61(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][61] = aa[k - 1][61] + f62(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][62] = aa[k - 1][62] + f63(q3, q4, Ui[i], Vi[i]) * fk;
//			aa[k - 1][63] = aa[k - 1][63] + f64(q3, q4, Ui[i], Vi[i]) * fk;
//			bb[k - 1] = bb[k - 1] + temp[i] * fk;
//		}
//		//		printf("k,bb[k-1],aa[k-1][0-1]= %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
//		//			k, bb[k - 1], aa[k - 1][0], aa[k - 1][1], aa[k - 1][2], aa[k - 1][3], aa[k - 1][4], aa[k - 1][5], aa[k - 1][6],
//		//			aa[k - 1][7], aa[k - 1][8], aa[k - 1][9], aa[k - 1][10], aa[k - 1][11], aa[k - 1][12], aa[k - 1][13], aa[k - 1][14],
//		//			aa[k - 1][15]);
//	}
//	for (i = 0; i < 64; i++)
//	{
//		bbx[i] = bb[i];
//		//		fprintf(out, "i,bb[i]= %d %le\n", i, bb[i]);
//		for (j = 0; j < 64; j++)
//		{
//			aax[i][j] = aa[i][j];
//			//			fprintf(out, "i,j,aa[i][j]= %d %d %le\n", i, j, aa[i][j]);
//		}
//	}
//	//	n = 16;
//	ciggj(aa, 64, bb);
//	for (i = 0; i < 64; i++)
//	{
//		d64[i] = bb[i];
//	}
//	return(1);
//	delete[] bb;
//	delete[] bbx;
//	for (i = 0; i < 64; i++)
//	{
//		delete aa[i];
//		delete aax[i];
//	}
//	delete[] aa;
//	delete[] aax;
//}
//int PDEcoor_new(int points, double q2, double q4, double q3, double q6, double* d64, double Ui, double Vi, double* temp)
//{
//	int i, j;
//	//	printf("is,ie,js,je,q2,q4= %d %d %d %d %le %le\n", is, ie, js, je, q2, q4);
//	for (int i = 0; i < points; i++)
//	{
//		temp[i] = d64[0] * f1(q2, q4, Ui, Vi) + d64[1] * f2(q2, q4, Ui, Vi) + d64[2] * f3(q2, q4, Ui, Vi) + d64[3]
//			* f4(q2, q4, Ui, Vi) + d64[4] * f5(q2, q4, Ui, Vi) + d64[5] * f6(q2, q4, Ui, Vi) + d64[6] * f7(q2, q4, Ui, Vi)
//			+ d64[7] * f8(q2, q4, Ui, Vi) + d64[8] * f9(q2, q4, Ui, Vi) + d64[9] * f10(q2, q4, Ui, Vi) + d64[10]
//			* f11(q2, q4, Ui, Vi) + d64[11] * f12(q2, q4, Ui, Vi) + d64[12] * f13(q2, q4, Ui, Vi) + d64[13]
//			* f14(q2, q4, Ui, Vi) + d64[14] * f15(q2, q4, Ui, Vi) + d64[15] * f16(q2, q4, Ui, Vi) + d64[16] * f17(q3, q6, Ui,Vi)
//			+ d64[17] * f18(q3, q6, Ui, Vi) + d64[18] * f19(q3, q6, Ui, Vi) + d64[19] * f20(q3, q6, Ui, Vi)+
//			d64[20] * f21(q3, q6, Ui, Vi) + d64[21] * f22(q3, q6, Ui, Vi) + d64[22] * f23(q3, q6, Ui, Vi) + 
//			d64[23] * f24(q3, q6, Ui, Vi) + d64[24] * f25(q3, q6, Ui, Vi) + d64[25] * f26(q3, q6, Ui, Vi)+
//			d64[26] * f27(q3, q6, Ui, Vi) + d64[27] * f28(q3, q6, Ui, Vi) + d64[28] * f29(q3, q6, Ui, Vi)+
//			d64[29] * f30(q3, q6, Ui, Vi) + d64[30] * f31(q3, q6, Ui, Vi) + d64[31] * f32(q3, q6, Ui, Vi) + 
//			d64[32] * f33(q2, q6, Ui, Vi) + d64[33] * f34(q2, q6, Ui, Vi) + d64[34] * f35(q2, q6, Ui, Vi) + d64[35]
//			* f36(q2, q6, Ui, Vi) + d64[36] * f37(q2, q6, Ui, Vi) + d64[37] * f38(q2, q6, Ui, Vi) + d64[38] * f39(q2, q6, Ui, Vi)
//			+ d64[39] * f40(q2, q6, Ui, Vi) + d64[40] * f41(q2, q6, Ui, Vi) + d64[41] * f42(q2, q6, Ui, Vi) + d64[42]
//			* f43(q2, q6, Ui, Vi) + d64[43] * f44(q2, q6, Ui, Vi) + d64[44] * f45(q2, q6, Ui, Vi) + d64[45]
//			* f46(q2, q6, Ui, Vi) + d64[46] * f47(q2, q6, Ui, Vi) + d64[47] * f48(q2, q6, Ui, Vi) + d64[48] * f49(q3, q4, Ui, Vi)
//			+ d64[49] * f50(q3, q4, Ui, Vi) + d64[50] * f51(q3, q4, Ui, Vi) + d64[51] * f52(q3, q4, Ui, Vi) +
//			d64[52] * f53(q3, q4, Ui, Vi) + d64[53] * f54(q3, q4, Ui, Vi) + d64[54] * f55(q3, q4, Ui, Vi) +
//			d64[55] * f56(q3, q3, Ui, Vi) + d64[56] * f57(q3, q4, Ui, Vi) + d64[57] * f58(q3, q4, Ui, Vi) +
//			d64[58] * f59(q3, q4, Ui, Vi) + d64[59] * f60(q3, q4, Ui, Vi) + d64[60] * f61(q3, q4, Ui, Vi) +
//			d64[61] * f62(q3, q4, Ui, Vi) + d64[62] * f63(q3, q4, Ui, Vi) + d64[63] * f64(q3, q4, Ui, Vi);
//	}
//	return(1);
//}
//
//
//
//
