///*
//cc -o fdm_surface11 fdm_surface11.c -lglut -lGLU -lGL -lXmu -lXext -lX11 -lm
//*/
//
//#include <GL/glut.h>
//#include <math.h>
//#include <stdlib.h>
//#include <stdio.h>
//#include <time.h>
//
//#define jb  2    //jb=10
//
//#define nu 15
//#define nv 15
//#define nuv 289
//#define nuvm 277
//
//#define nn 3
//#define icurve 10
//
//int ciggj(double **,int,double *);
//int fitting(int, int, int, int, double, double, double**, double**, double**, double*);
//int PDEcoor(int, int, int, int, double, double, double*, double**, double**, double**);
//
//FILE *in1,*out,*ToMaya;
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
//GLuint startList;
//
//void CALLBACK errorCallback(GLenum errorCode)
//{
//   const GLubyte *estring;
//
//   estring = gluErrorString(errorCode);
//   fprintf(stderr, "Quadric Error: %s\n", estring);
//   exit(0);
//}
//
//void init(void)
//{
//	int i, is, j, js, k, n, n0, l, m, im[100], iu, iv, jj, iuu, icurve, ie, je;
//	double fk, du, dLj, L, ui, vj, mFaceWidth, xori[100][100], yori[100][100], zori[100][100], Lj[100],
//		d16[16], dx16[16], dy16[16], dz16[16], q2, q4,aax[16][16],bbx[16],xone[500], yone[500], zone[500], xtwo[200], ytwo[200], ztwo[200], EoneM, EoneA,
//		EtwoM, EtwoA,EoneC, EtwoC;
// X83[][0]=x component, X83[][1]=y component, X83[][2]=zx component; U82[][0]= u component, U82[][1]= v component.
//	char char1, char2; uij[100][100], vij[100][100],X100[100][100], 
//	char szLine[256];
//
//   clock_t start, finish, duration;
//    GLUquadricObj *qobj;
//    GLfloat mat_ambient[] = { 0.8, 0.8, 0.8, 1.0 };
//    GLfloat mat_specular[] = { 0.3, 0.3, 0.3, 1.0 };
//    GLfloat mat_shininess[] = { 100.0 };
//    GLfloat light_position[] = { 1000.0, 70000.0, -100000.0, 0.0 };
//    GLfloat model_ambient[] = { 0.9, 0.9, 0.9, 1.0 };
//    glClearColor(0.382, 0.382, 0.382, 1.0);
//	glClearColor(1.0, 1.0, 1.0, 1.0);
//
//    glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
//    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
//    glLightfv(GL_LIGHT0, GL_POSITION, light_position);
//    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, model_ambient);
//
//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);
//    glEnable(GL_DEPTH_TEST);
//
//
//    startList = glGenLists(1);
//    qobj = gluNewQuadric();
//    gluQuadricDrawStyle(qobj, GLU_FILL); 
//    gluQuadricNormals(qobj, GLU_SMOOTH);
//    glNewList(startList, GL_COMPILE);
//
//	double** uij = new double* [100]; // 
//	double** vij = new double* [100]; // 
//	double** X100 = new double* [100]; // 
//	double** xcal = new double* [100]; // 
//	double** ycal = new double* [100]; // 
//	double** zcal = new double* [100]; // 
//	for (i = 0; i < 100; i++)
//	{
//		uij[i] = new double[100];
//		vij[i] = new double[100];
//		X100[i] = new double[100];
//		xcal[i] = new double[100];
//		ycal[i] = new double[100];
//		zcal[i] = new double[100];
//	}
//
//	n=2*jb+1;
//	n = 2 * jb - 1;
//	n0=2*jb;
//
//    in1 = fopen("In1.txt","rt"); 
//    out = fopen("out.txt","wt");
//	ToMaya = fopen("curve.mel", "wt");
//
//	fscanf(in1, "%d\n", &icurve);
//	printf("icurve= %d\n", icurve);
//
//	for (i = 0; i < icurve; i++)
//	{
//		fscanf(in1, "%d\n", &im[i]);
//		printf("i,im[i]= %d %d\n", i,im[i]);
//		fprintf(ToMaya, "curve -d 1\n");
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//					fscanf(in1,"%le %s %le %s %le\n",&xori[j],char1,&yori[j],char2,&zori[j]);
//			fscanf(in1, "%le %le %le\n", &xori[i][j], &yori[i][j], &zori[i][j]);
//			fprintf(ToMaya, "-p   %le %le %le\n", xori[i][j], yori[i][j], zori[i][j]);
//			printf("i,j-x-y-z=  %d %d %le %le %le\n", i, j, xori[i][j], yori[i][j], zori[i][j]);
//		}
//		fprintf(ToMaya, ";\n");
//	}
//
//	for (i = 0; i < 100; i++)
//	{
//		for (j = 0; j < 100; j++)
//		{
//			uij[i][j] = 0.0;
//			vij[i][j] = 0.0;
//		}
//	}
//
//	du = 1.0 / 65;
//	mFaceWidth = 0.0;
//	for (i = 0; i < icurve; i++)
//	{
//		ui = i * du;
//		L = 0.0;
//		Lj[0] = 0.0; 
//		for (j = 1; j <= im[i] - 1; j++)
//		{
//			dLj = sqrt((xori[i][j]- xori[i][j-1])* (xori[i][j] - xori[i][j - 1]) + (yori[i][j] - yori[i][j - 1]) 
//				* (yori[i][j] - yori[i][j - 1]) + (zori[i][j] - zori[i][j - 1]) * (zori[i][j] - zori[i][j - 1]));
//			L = L + dLj;
//			Lj[j] = Lj[j - 1] + dLj;
//		}
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			uij[i][j] = ui;
//			vj = Lj[j] / L;
//			vij[i][j] = vj;
//		}
//		if (L > mFaceWidth)mFaceWidth = L;
//	}
//
//	for (i = 0; i < icurve; i++)
//	{
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			fprintf(out, "i,j,ui= %d %d %le\n", i,j,uij[i][j]);
//		}
//	}
//
//	for (i = 0; i < icurve; i++)
//	{
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			fprintf(out, "i,j,vj= %d %d %le\n", i, j, vij[i][j]);
//		}
//	}
//
// PDE reconstruction
//	q2 = pi / 2.0;
//	q4 = pi;
//	q2 = 0.1;
//	q4 = 0.2;
//	q2 = 0.05;
//	q4 = 0.1;
//	q2 = 0.075;
//	q4 = 0.15;
//	q2 = 0.07;
//	q4 = 0.07;
//	q2 = 0.15;
//	q4 = 0.15;
//	ie = 5;
//	je = 5;
//	q2 = 0.1;
//	q4 = 0.1;
// First test
// FourByFour
//	is = 12;
//	ie = 15;
//	js = 15;
//	je = 18;
//  FiveByFive
//	is = 12;
//	ie = 16;
//	js = 15;
//	je = 19;
//  SixBySix
//	is = 12;
//	ie = 17;
//	js = 14;
//	je = 19;
//
// second test
// FourByFour
//	is = 32;
//	ie = 35;
//	js = 18;
//	je = 21;
//
//  //FiveByFive
//	is = 32;
//	ie = 36;
//	js = 18;
//	je = 22;
//
// // SixBySix
//	is = 32;
//	ie = 37;
//	js = 17;
//	je = 22;
///*
//	  SixBySix
//	is = 12;
//	ie = 18;
//	js = 14;
//	je = 20;
//*/
//	for (i = is; i <= ie; i++)
//	{
//		fprintf(ToMaya, "curve -d 3\n");
//		for (j = js; j <= je; j++)
//		{
//			fprintf(ToMaya, "-p   %le %le %le\n", xori[i][j], yori[i][j], zori[i][j]);
//		}
//		fprintf(ToMaya, ";\n");
//	}
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			fprintf(out, "i,j,xo,yo,zo= %d %d %le %le %le\n", i, j, xori[i][j], 
//				yori[i][j], zori[i][j]);
//		}
//	}
//	exit(1);
//  //x component
//	for (i = 0; i < icurve; i++)
//	{
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			X100[i][j]=xori[i][j];
//		}
//	}
//	for (i = 0; i < icurve; i++)
//	{
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			fprintf(out, "i,j,X100= %d %d %le\n", i, j, X100[i][j]);
//		}
//	}
//	fitting(is, ie, js, je, q2, q4, uij, vij, X100, d16);
//	for (i = 0; i < 16; i++)
//	{
//		dx16[i] = d16[i];
//	}
//  //y component
//	for (i = 0; i < icurve; i++)
//	{
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			X100[i][j] = yori[i][j];
//		}
//	}
//	fitting(is, ie, js, je, q2, q4, uij, vij, X100, d16);
//	for (i = 0; i < 16; i++)
//	{
//		dy16[i] = d16[i];
//	}
// // z component
//	for (i = 0; i < icurve; i++)
//	{
//		for (j = 0; j <= im[i] - 1; j++)
//		{
//			X100[i][j] = zori[i][j];
//		}
//	}
//	fitting(is, ie, js, je, q2, q4, uij, vij, X100, d16);
//	for (i = 0; i < 16; i++)
//	{
//		dz16[i] = d16[i];
//	}
// //Calculate the coordinate values of PDE surfaces
//	PDEcoor(is, ie, js, je, q2, q4, dx16, uij, vij, X100);
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			xcal[i][j] = X100[i][j];
//		}
//	}
//	PDEcoor(is, ie, js, je, q2, q4, dy16, uij, vij, X100);
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			ycal[i][j] = X100[i][j];
//		}
//	}
//	PDEcoor(is, ie, js, je, q2, q4, dz16, uij, vij, X100);
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			zcal[i][j] = X100[i][j];
//		}
//	}
//	for (i = is; i <= ie; i++)
//	{
//		fprintf(ToMaya, "curve -d 3\n");
//		for (j = js; j <= je; j++)
//		{
//			fprintf(ToMaya, "-p   %le %le %le\n", xcal[i][j], ycal[i][j], zcal[i][j]);
//		}
//		fprintf(ToMaya, ";\n");
//	}
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			fprintf(out, "i,j,xo,xc,yo,yc,zo,zc= %d %d %le %le %le %le %le %le %le %le %le\n", i, j, xori[i][j], 
//				(xori[i][j] - xcal[i][j]), (xori[i][j] - xcal[i][j]) / mFaceWidth, yori[i][j], (yori[i][j] - ycal[i][j]),
//				(yori[i][j] - ycal[i][j]) / mFaceWidth, zori[i][j], (zori[i][j] - zcal[i][j]), (zori[i][j] - zcal[i][j])
//				/ mFaceWidth);
//			fprintf(out, "i,j,xo,xc,yo,yc,zo,zc= %d %d %le %le %le %le %le %le\n", i, j, xori[i][j], (xori[i][j]-xcal[i][j])/ mFaceWidth,
//				yori[i][j], (yori[i][j]-ycal[i][j])/ mFaceWidth, zori[i][j], (zori[i][j]-zcal[i][j])/ mFaceWidth);
//		}
//	}
//
//	double Err, ErrA, ErrM;
// //error calculation
//	ErrM = 0.0;
//	ErrA = 0.0;
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			Err = sqrtf((xori[i][j] - xcal[i][j]) * (xori[i][j] - xcal[i][j]) + (yori[i][j] - ycal[i][j])
//				* (yori[i][j] - ycal[i][j]) + (zori[i][j] - zcal[i][j]) * (zori[i][j] - zcal[i][j]));
//			if (Err > ErrM)ErrM = Err;
//			ErrA = ErrA + Err;
//		}
//	}
//	ErrA = ErrA / (ie + 1 - is) / (je + 1 - js); 
//	printf("ie + 1 - is,je + 1 - js,ErrM,ErrA= %d %d %le %le\n", ie + 1 - is, je + 1 - js, ErrM, ErrA);
//
//	for (i = 0; i < 100; i++)
//	{
//		delete uij[i];
//		delete vij[i];
//		delete X100[i];
//		delete xcal[i];
//		delete ycal[i];
//		delete zcal[i];
//	}
//	delete[] uij;
//	delete[] vij;
//	delete[] X100;
//	delete[] xcal;
//	delete[] ycal;
//	delete[] zcal;
//
//	fclose(in1);
//	fclose(out);
//	fclose(ToMaya);
//
//	printf("Pass 0\n");
//
//	glEndList();
//
//}
//
//void display(void)
//{
//
//   glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//   glPushMatrix();
//
//   glEnable(GL_COLOR_MATERIAL);
//   glEnable (GL_LIGHTING);
//   glShadeModel (GL_SMOOTH);
//
//   glColor3f(0.5, 0.0, 1.0);
//   glTranslatef(0.0, 0.0, 0.0);
//   glPushMatrix();
//   glScalef(0.15, 0.15, 0.15);
//   glRotatef(45.0,45.0, 45.0, 1.0);
//   glRotatef(90.0, 90.0, 0.0, 1.0);
//
//   glCallList(startList);
//   glPopMatrix();
//
//   glFlush();
//}
//
//void reshape (int w, int h)
//{
//   glViewport(0, 0, (GLsizei) w, (GLsizei) h);
//   glMatrixMode(GL_PROJECTION);
//   glLoadIdentity();
//   if (w <= h)
//      glOrtho(-2.5, 2.5, -2.5*(GLfloat)h/(GLfloat)w,
//         2.5*(GLfloat)h/(GLfloat)w, -10.0, 10.0);
//   else
//      glOrtho(-2.5*(GLfloat)w/(GLfloat)h,
//         2.5*(GLfloat)w/(GLfloat)h, -2.5, 2.5, -10.0, 10.0);
//   glMatrixMode(GL_MODELVIEW);
//   glLoadIdentity();
//}
//
//void keyboard(unsigned char key, int x, int y)
//{
//   switch (key) {
//      case 27:
//         exit(0);
//         break;
//   }
//}
//
//int main(int argc, char** argv)
//{
//   glutInit(&argc, argv);
//   glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB | GLUT_DEPTH);
//   glutInitWindowSize(500, 500); 
//   glutInitWindowPosition(100, 100);
//   glutCreateWindow(argv[0]);
//   init();
//   glutDisplayFunc(display); 
//   glutReshapeFunc(reshape);
//   glutKeyboardFunc(keyboard);
//   glutMainLoop();
//   return 0;
//}
//
//int ciggj(double **a,int n,double *b)      
//{
//	int i, j, k, is, u, v;
//	int *js = new int [n];
//	double d, t;
//	for (k=0;k<=n-1;k++)
//	{
//		d=0.0;
//		for (i=k;i<=n-1;i++)
//		for (j=k;j<=n-1;j++)
//			{t=fabs(a[i][j]);			
//			if(t>d) {d=t;js[k]=j;is=i;}
//			}
//			if(d+1.0==1.0)
//				{delete js;printf("fail\n");return(0);}
//			if (is!=k)
//				{for(j=k;j<=n-1;j++)
//					{u=k*n+j;v=is*n+j;
//					 t=a[k][j];a[k][j]=a[is][j];a[is][j]=t;
//					 }
//				  t=b[k];b[k]=b[is];b[is]=t;
//				 }
//			if (js[k]!=k)
//			for (i=0;i<=n-1;i++)
//				{u=i*n+k;v=i*n+js[k];
//				 t=a[i][k];a[i][k]=a[i][js[k]];a[i][js[k]]=t;
//				}
//			t=a[k][k];
//			for (j=k+1;j<=n-1;j++)
//				{u=k*n+j;
//				 if (a[k][j]!=0.0)a[k][j]=a[k][j]/t;
//				}
//			b[k]=b[k]/t;
//			for (j=k+1;j<=n-1;j++)
//				{u=k*n+j;
//				 if(a[k][j]!=0.0)
//				 	{for (i=0;i<=n-1;i++)
//				 		{v=i*n+k;
//				 		 if((i!=k) && (a[i][k]!=0.0))
//				 		 	{is=i*n+j;
//				 		 	 a[i][j]=a[i][j]-a[i][k]*a[k][j];
//				 		 	}
//				 		 }
//				 	}
//				 }
//				 for(i=0;i<=n-1;i++)
//				 	{u=i*n+k;
//				 	 if((i!=k) && (a[i][k]!=0.0))
//				 	 	b[i]=b[i]-a[i][k]*b[k];
//				 	 }
//			}
//			for (k=n-1;k>=0;k--)
//				if(k!=js[k])
//					{t=b[k];b[k]=b[js[k]];b[js[k]]=t;}
//				delete js;
//				return(1);
//			}
//
//int fitting(int is, int ie, int js, int je, double q2, double q4, double **uij, double **vij, double **X100, double *d16)
//{
//	int i, j, k;
//	double ui, vj, fk;
//	double* bb = new double[16]; // 
//	double* bbx = new double[16]; // 
//	double** aa = new double* [16]; // 
//	double** aax = new double* [16]; // 
//	for (i = 0; i < 16; i++)
//	{
//		aa[i] = new double[16];
//		aax[i] = new double[16];
//	}
//	printf("is,ie,js,je,q2,q4= %d %d %d %d %le %le\n", is, ie, js, je, q2, q4);
//	for (i = 0; i < 16; i++)
//	{
//		bb[i] = 0.0;
//		for (j = 0; j < 16; j++)
//		{
//			aa[i][j] = 0.0;
//		}
//	}
//	for (k = 1; k <= 16; k++)
//	{
//		for (i = is; i <= ie; i++)
//		{
//			for (j = js; j <= je; j++)
//			{
//				ui = (uij[i][j] - uij[is][j]) / (uij[ie][j] - uij[is][j]);
//				vj = (vij[i][j] - vij[i][js]) / (vij[i][je] - vij[i][js]);
//				if (k == 1)fk = f1(q2, q4, ui, vj);
//				if (k == 2)fk = f2(q2, q4, ui, vj);
//				if (k == 3)fk = f3(q2, q4, ui, vj);
//				if (k == 4)fk = f4(q2, q4, ui, vj);
//				if (k == 5)fk = f5(q2, q4, ui, vj);
//				if (k == 6)fk = f6(q2, q4, ui, vj);
//				if (k == 7)fk = f7(q2, q4, ui, vj);
//				if (k == 8)fk = f8(q2, q4, ui, vj);
//				if (k == 9)fk = f9(q2, q4, ui, vj);
//				if (k == 10)fk = f10(q2, q4, ui, vj);
//				if (k == 11)fk = f11(q2, q4, ui, vj);
//				if (k == 12)fk = f12(q2, q4, ui, vj);
//				if (k == 13)fk = f13(q2, q4, ui, vj);
//				if (k == 14)fk = f14(q2, q4, ui, vj);
//				if (k == 15)fk = f15(q2, q4, ui, vj);
//				if (k == 16)fk = f16(q2, q4, ui, vj);
//								printf("k,ui,vj,fk= %d %le %le %le\n", k, ui, vj, fk);
//				aa[k - 1][0] = aa[k - 1][0] + f1(q2, q4, ui, vj) * fk;
//				aa[k - 1][1] = aa[k - 1][1] + f2(q2, q4, ui, vj) * fk;
//				aa[k - 1][2] = aa[k - 1][2] + f3(q2, q4, ui, vj) * fk;
//				aa[k - 1][3] = aa[k - 1][3] + f4(q2, q4, ui, vj) * fk;
//				aa[k - 1][4] = aa[k - 1][4] + f5(q2, q4, ui, vj) * fk;
//				aa[k - 1][5] = aa[k - 1][5] + f6(q2, q4, ui, vj) * fk;
//				aa[k - 1][6] = aa[k - 1][6] + f7(q2, q4, ui, vj) * fk;
//				aa[k - 1][7] = aa[k - 1][7] + f8(q2, q4, ui, vj) * fk;
//				aa[k - 1][8] = aa[k - 1][8] + f9(q2, q4, ui, vj) * fk;
//				aa[k - 1][9] = aa[k - 1][9] + f10(q2, q4, ui, vj) * fk;
//				aa[k - 1][10] = aa[k - 1][10] + f11(q2, q4, ui, vj) * fk;
//				aa[k - 1][11] = aa[k - 1][11] + f12(q2, q4, ui, vj) * fk;
//				aa[k - 1][12] = aa[k - 1][12] + f13(q2, q4, ui, vj) * fk;
//				aa[k - 1][13] = aa[k - 1][13] + f14(q2, q4, ui, vj) * fk;
//				aa[k - 1][14] = aa[k - 1][14] + f15(q2, q4, ui, vj) * fk;
//				aa[k - 1][15] = aa[k - 1][15] + f16(q2, q4, ui, vj) * fk;
//				bb[k - 1] = bb[k - 1] + X100[i][j] * fk;
//			}
//		}
//		printf("k,bb[k-1],aa[k-1][0-1]= %d %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le %le\n",
//			k, bb[k - 1], aa[k - 1][0], aa[k - 1][1], aa[k - 1][2], aa[k - 1][3], aa[k - 1][4], aa[k - 1][5], aa[k - 1][6],
//			aa[k - 1][7], aa[k - 1][8], aa[k - 1][9], aa[k - 1][10], aa[k - 1][11], aa[k - 1][12], aa[k - 1][13], aa[k - 1][14],
//			aa[k - 1][15]);
//	}
//	for (i = 0; i < 16; i++)
//	{
//		bbx[i] = bb[i];
//		fprintf(out, "i,bb[i]= %d %le\n", i, bb[i]);
//		for (j = 0; j < 16; j++)
//		{
//			aax[i][j] = aa[i][j];
//			fprintf(out, "i,j,aa[i][j]= %d %d %le\n", i, j, aa[i][j]);
//		}
//	}
//		n = 16;
//	ciggj(aa, 16, bb);
//	for (i = 0; i < 16; i++)
//	{
//		d16[i] = bb[i];
//	}
//	return(1);
//	delete[] bb;
//	delete[] bbx;
//	for (i = 0; i < 16; i++)
//	{
//		delete aa[i];
//		delete aax[i];
//	}
//	delete[] aa;
//	delete[] aax;
//}
//
//int PDEcoor(int is, int ie, int js, int je, double q2, double q4, double* d16, double** uij, double** vij, double** X100)
//{
//	int i, j;
//	double ui, vj;
//	printf("is,ie,js,je,q2,q4= %d %d %d %d %le %le\n", is, ie, js, je, q2, q4);
//	for (i = is; i <= ie; i++)
//	{
//		for (j = js; j <= je; j++)
//		{
//			ui = (uij[i][j] - uij[is][j]) / (uij[ie][j] - uij[is][j]);
//			vj = (vij[i][j] - vij[i][js]) / (vij[i][je] - vij[i][js]);
//			X100[i][j] = d16[0] * f1(q2, q4, ui, vj) + d16[1] * f2(q2, q4, ui, vj) + d16[2] * f3(q2, q4, ui, vj) + d16[3]
//				* f4(q2, q4, ui, vj) + d16[4] * f5(q2, q4, ui, vj) + d16[5] * f6(q2, q4, ui, vj) + d16[6] * f7(q2, q4, ui, vj)
//				+ d16[7] * f8(q2, q4, ui, vj) + d16[8] * f9(q2, q4, ui, vj) + d16[9] * f10(q2, q4, ui, vj) + d16[10]
//				* f11(q2, q4, ui, vj) + d16[11] * f12(q2, q4, ui, vj) + d16[12] * f13(q2, q4, ui, vj) + d16[13]
//				* f14(q2, q4, ui, vj) + d16[14] * f15(q2, q4, ui, vj) + d16[15] * f16(q2, q4, ui, vj);
//		}
//	}
//	return(1);
//}
//
