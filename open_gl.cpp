//
///* Copyright (c) Mark J. Kilgard, 1995. */
//
///* This program is freely distributable without licensing fees
//   and is provided without guarantee or warrantee expressed or
//   implied. This program is -not- in the public domain. */
//
//   /* molehill uses the GLU NURBS routines to draw some nice surfaces. */
//
//#include <GL/glut.h>
//#include<math.h>
//#include <stdlib.h>
//#include<stdio.h>
//#include "Header.h"
//
//
//inline double f1(double q2, double q4, double u, double v);
//inline double f2(double q2, double q4, double u, double v);
//inline double f3(double q2, double q4, double u, double v);
//inline double f4(double q2, double q4, double u, double v);
//inline double f5(double q2, double q4, double u, double v);
//inline double f6(double q2, double q4, double u, double v);
//inline double f7(double q2, double q4, double u, double v);
//inline double f8(double q2, double q4, double u, double v);
//inline double f9(double q2, double q4, double u, double v);
//inline double f10(double q2, double q4, double u, double v);
//inline double f11(double q2, double q4, double u, double v);
//inline double f12(double q2, double q4, double u, double v);
//inline double f13(double q2, double q4, double u, double v);
//inline double f14(double q2, double q4, double u, double v);
//inline double f15(double q2, double q4, double u, double v);
//inline double f16(double q2, double q4, double u, double v);
//int PDEcoor_new(int points, double q2, double q4, double* d16, double Ui, double Vi, double* temp);
//void init(void);
//const float PI = 3.141592f;
//GLfloat mat_red_diffuse[] = { 0.7, 0.0, 0.1, 1.0 };
//GLfloat mat_green_diffuse[] = { 0.0, 0.7, 0.1, 1.0 };
//GLfloat mat_blue_diffuse[] = { 0.0, 0.1, 0.7, 1.0 };
//GLfloat mat_yellow_diffuse[] = { 0.7, 0.8, 0.1, 1.0 };
//GLfloat mat_specular[] = { 1.0, 1.0, 1.0, 1.0 };
//GLfloat mat_shininess[] = { 100.0 };
//GLfloat knots[8] = { 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0 };
//GLfloat pts1[4][4][3], pts2[4][4][3];
//GLfloat pts3[4][4][3], pts4[4][4][3];
//GLUnurbsObj* nurb;
//int u, v;
//
////void drawSphere(double r, int lats, int longs) {
////    int i, j;
////    for (i = 0; i <= lats; i++) {
////        double lat0 = PI * (-0.5 + (double)(i - 1) / lats);
////        double z0 = sin(lat0);
////        double zr0 = cos(lat0);
////
////        double lat1 = PI * (-0.5 + (double)i / lats);
////        double z1 = sin(lat1);
////        double zr1 = cos(lat1);
////
////        glBegin(GL_QUAD_STRIP);
////        for (j = 0; j <= longs; j++) {
////            double lng = 2 * PI * (double)(j - 1) / longs;
////            double x = cos(lng);
////            double y = sin(lng);
////
////            glNormal3f(x * zr0, y * zr0, z0);
////            glVertex3f(r * x * zr0, r * y * zr0, r * z0);
////            glNormal3f(x * zr1, y * zr1, z1);
////            glVertex3f(r * x * zr1, r * y * zr1, r * z1);
////        }
////        glEnd();
////    }
////}
//double x, y, z, alpha, beta; // Storage for coordinates and angles        
//double radius = 5.0f;
//int gradation = 20;
//double radius_gl = 0.0f;
//extern double dx16[16], dy16[16], dz16[16];
//void drawSphere() {
//    double* x_coor_gl= new double[1800];
//    double* y_coor_gl = new double[1800];
//    double* z_coor_gl = new double[1800];
//
//    int I = 20;
//    int J = 20;
//    double delta_u = 1.0 / I;
//    double delta_v = 1.0 / J;
//    for (int i = 0; i <= I; i++)
//    {
//        glBegin(GL_QUAD_STRIP);
//        double u_i = i * delta_u;
//        for (int j = 0; j <= J; j++)
//        {
//            double v_j = j * delta_v;
//            PDEcoor_new(1730, 0.1, 0.1, dx16, u_i, v_j, x_coor_gl);
//            PDEcoor_new(1730, 0.1, 0.1, dy16, u_i, v_j, y_coor_gl);
//            PDEcoor_new(1730, 0.1, 0.1, dz16, u_i, v_j, z_coor_gl);
//            radius_gl = sqrt(pow(*x_coor_gl, 2) + pow(*y_coor_gl, 2) + pow(*z_coor_gl, 2));
//            x = *x_coor_gl;
//            y = *y_coor_gl;
//            z = *z_coor_gl;
//            glNormal3f(-x / radius_gl, -y / radius_gl, -z / radius_gl);
//            glVertex3f(x*3.0, y*3.0, z*3.0);
//            PDEcoor_new(1730, 0.1, 0.1, dx16, u_i +delta_u, v_j, x_coor_gl);
//            PDEcoor_new(1730, 0.1, 0.1, dy16, u_i +delta_u, v_j , y_coor_gl);
//            PDEcoor_new(1730, 0.1, 0.1, dz16, u_i+delta_u, v_j , z_coor_gl);
//            radius_gl = sqrt(pow(*x_coor_gl, 2) + pow(*y_coor_gl, 2) + pow(*z_coor_gl, 2));
//            x = *x_coor_gl;
//            y = *y_coor_gl;
//            z = *z_coor_gl;
//            glNormal3f(-x / radius_gl, -y / radius_gl, -z / radius_gl);
//            glVertex3f(x*3.0,y*3.0, z*3.0);
//            //fprintf(realPointsCoor, "%le %le %le\n", *x_cal, *y_cal, *z_cal);
//        }
//        glEnd();
//    }
//    //for (alpha = 0.0; alpha < PI; alpha += PI / gradation)
//    //{
//    //    glBegin(GL_TRIANGLE_STRIP);
//    //    //glBegin(GL_QUAD_STRIP);
//    //    for (beta = 0.0; beta < 2.01 * PI; beta += PI / gradation)
//    //    {
//    //        x = radius * cos(beta) * sin(alpha);
//    //        y = radius * sin(beta) * sin(alpha);
//    //        z = radius * cos(alpha);
//    //        glNormal3f(x/radius, y/radius, z/radius);
//    //        glVertex3f(x, y, z);
//    //        x = radius * cos(beta) * sin(alpha + PI / gradation);
//    //        y = radius * sin(beta) * sin(alpha + PI / gradation);
//    //        z = radius * cos(alpha + PI / gradation);
//    //        glNormal3f(x / radius, y / radius, z / radius);
//    //        glVertex3f(x, y, z);
//    //    }
//    //    glEnd();
//    //}
//    delete[] x_coor_gl;
//    delete[] y_coor_gl;
//    delete[] z_coor_gl;
//  
//}
//static void
//display(void)
//{
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    glCallList(1);
//    glFlush();
//}
//
//int
//main(int argc, char** argv)
//{
//    glutInit(&argc, argv);
//    glutInitWindowSize(1080, 960);
//    glutCreateWindow("molehill");
//    glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
//    glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
//    glEnable(GL_LIGHTING);
//    glEnable(GL_LIGHT0);
//    //glEnable(GL_DEPTH_TEST);
//    glEnable(GL_AUTO_NORMAL);
//    glEnable(GL_NORMALIZE);
//    nurb = gluNewNurbsRenderer();
//    gluNurbsProperty(nurb, GLU_SAMPLING_TOLERANCE, 25.0);
//    gluNurbsProperty(nurb, GLU_DISPLAY_MODE, GLU_FILL);
//
//    /* Build control points for NURBS mole hills. */
//    for (u = 0; u < 4; u++) {
//        for (v = 0; v < 4; v++) {
//            /* Red. */
//            pts1[u][v][0] = 2.0 * ((GLfloat)u);
//            pts1[u][v][1] = 2.0 * ((GLfloat)v);
//            if ((u == 1 || u == 2) && (v == 1 || v == 2))
//                /* Stretch up middle. */
//                pts1[u][v][2] = 6.0;
//            else
//                pts1[u][v][2] = 0.0;
//
//            /* Green. */
//            pts2[u][v][0] = 2.0 * ((GLfloat)u - 3.0);
//            pts2[u][v][1] = 2.0 * ((GLfloat)v - 3.0);
//            if ((u == 1 || u == 2) && (v == 1 || v == 2))
//                if (u == 1 && v == 1)
//                    /* Pull hard on single middle square. */
//                    pts2[u][v][2] = 15.0;
//                else
//                    /* Push down on other middle squares. */
//                    pts2[u][v][2] = -2.0;
//            else
//                pts2[u][v][2] = 0.0;
//
//            /* Blue. */
//            pts3[u][v][0] = 2.0 * ((GLfloat)u - 3.0);
//            pts3[u][v][1] = 2.0 * ((GLfloat)v);
//            if ((u == 1 || u == 2) && (v == 1 || v == 2))
//                if (u == 1 && v == 2)
//                    /* Pull up on single middple square. */
//                    pts3[u][v][2] = 11.0;
//                else
//                    /* Pull up slightly on other middle squares. */
//                    pts3[u][v][2] = 2.0;
//            else
//                pts3[u][v][2] = 0.0;
//
//            /* Yellow. */
//            pts4[u][v][0] = 2.0 * ((GLfloat)u);
//            pts4[u][v][1] = 2.0 * ((GLfloat)v - 3.0);
//            if ((u == 1 || u == 2 || u == 3) && (v == 1 || v == 2))
//                if (v == 1)
//                    /* Push down front middle and right squares. */
//                    pts4[u][v][2] = -2.0;
//                else
//                    /* Pull up back middle and right squares. */
//                    pts4[u][v][2] = 5.0;
//            else
//                pts4[u][v][2] = 0.0;
//        }
//    }
//    /* Stretch up red's far right corner. */
//    pts1[3][3][2] = 6;
//    /* Pull down green's near left corner a little. */
//    pts2[0][0][2] = -2;
//    /* Turn up meeting of four corners. */
//    pts1[0][0][2] = 1;
//    pts2[3][3][2] = 1;
//    pts3[3][0][2] = 1;
//    pts4[0][3][2] = 1;
//
//    glMatrixMode(GL_PROJECTION);
//    gluPerspective(55.0, 1.0, 2.0, 24.0);
//    glMatrixMode(GL_MODELVIEW);
//    glTranslatef(0.0, 0.0, -15.0);
//    glRotatef(330.0, 1.0, 0.0, 0.0);
//
//    glNewList(1, GL_COMPILE);
//    /* Render red hill. */
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_red_diffuse);
//    gluBeginSurface(nurb);
//    //gluNurbsSurface(nurb, 8, knots, 8, knots,
//     //   4 * 3, 3, &pts1[0][0][0],
//       // 4, 4, GL_MAP2_VERTEX_3);
//    gluEndSurface(nurb);
//
//    /* Render green hill. */
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_green_diffuse);
//    gluBeginSurface(nurb);
//    //gluNurbsSurface(nurb, 8, knots, 8, knots,
//      //  4 * 3, 3, &pts2[0][0][0],
//       // 4, 4, GL_MAP2_VERTEX_3);
//    gluEndSurface(nurb);
//
//    /* Render blue hill. */
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_blue_diffuse);
//    gluBeginSurface(nurb);
//    //gluNurbsSurface(nurb, 8, knots, 8, knots,
//       // 4 * 3, 3, &pts3[0][0][0],
//        //4, 4, GL_MAP2_VERTEX_3);
//    gluEndSurface(nurb);
//    //drawSphere(4.0, 20, 20);
//    init();
//    drawSphere();
//    /* Render yellow hill. */
//    glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_yellow_diffuse);
//    gluBeginSurface(nurb);
//    //gluNurbsSurface(nurb, 8, knots, 8, knots,
//       // 4 * 3, 3, &pts4[0][0][0],
//        //4, 4, GL_MAP2_VERTEX_3);
//    gluEndSurface(nurb);
//    glEndList();
//
//    glutDisplayFunc(display);
//    glutMainLoop();
//    return 0;             /* ANSI C requires main to return int. */
//}