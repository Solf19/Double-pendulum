#include <iostream>
#include<cmath>
#include <iomanip>
#include <fstream>
using namespace std;
void dYdt(double, double*, double*);
void RK4(double , double *, double , int , void (*Func)
        (double, double*, double*));
void RK2(double , double *, double , int , void (*Func)
        (double, double*, double*));
#define STAGE 5
#define PARAM 4
double g_m = 1.; // m1 = m2 = m
double g_l = 1.; // l1 = l2 = l
double g_g = 9.8; 

int main(){
    double ti = 0., t = ti, tf;
    int neq = 4;
    int m;
    double Y[neq], Y1[neq];
   
#if STAGE == 1
/* comparison of two different time steps to see how the relative
    error of the total energy evolves */
    //Fixed initial conditions
    Y[0] = 0.9; Y[1] = 0.3; Y[2] = 0.2; Y[3] = 0.;
    //Time steps
    double dt1, dt2;
    dt1 = 0.05; dt2 = 0.01;
    ofstream file1("aa_stage1_dt1.txt");
    ofstream file2("aa_stage1_dt2.txt");
    tf = 200.;
    m = (tf - ti) / dt1;
    double en_i, en;
    en_i = g_m * g_l * g_l * (Y[2] * Y[2] + 0.5 * Y[3] * Y[3] + 
        Y[2] * Y[3] * cos(Y[0] - Y[1])) - g_m * g_g * g_l * 
        (cos(Y[1]) + 2. * cos(Y[0]));
    for(int i = 0; i < m; i++){
        RK4(t, Y, dt1, neq, dYdt);
        t += dt1;
        en = g_m * g_l * g_l * (Y[2] * Y[2] + 0.5 * Y[3] * Y[3] + 
             Y[2] * Y[3] * cos(Y[0] - Y[1])) - g_m * g_g * g_l * 
             (cos(Y[1]) + 2. * cos(Y[0]));
        file1 << t << "   " << fabs((en_i - en) / en_i) << endl;
    }

    Y[0] = 0.9; Y[1] = 0.3; Y[2] = 0.2; Y[3] = 0.;
    t = 0;
    m = (tf - ti) / dt2;
    for(int j = 0; j < m; j++){
        RK4(t, Y, dt2, neq, dYdt);
        t += dt2;
        en = g_m * g_l * g_l * (Y[2] * Y[2] + 0.5 * Y[3] * Y[3] + 
             Y[2] * Y[3] * cos(Y[0] - Y[1])) - g_m * g_g * g_l * 
             (cos(Y[1]) + 2. * cos(Y[0]));
        file2 << t << "   " << fabs((en_i - en) / en_i) << endl;
    }

    file1.close();
    file2.close();
#endif

#if STAGE == 2
/*This wants to show how relative energy at a fixed time t and 
  fixed intial conditions varies with different time steps*/
  //Fixed initial conditions
    Y[0] = 3.; Y[1] = -2.5; Y[2] = 1.; Y[3] = 0.;
    //Time steps
    double  dt_min = 0.0001, dt = dt_min, dt_max = 0.05;
    int m = 1000000, n = 100;
    ofstream file1("aa_stage2_en_err(dt).txt");
    double en_i, en;
    en_i = g_m * g_l * g_l * (Y[2] * Y[2] + 0.5 * Y[3] * Y[3] + 
           Y[2] * Y[3] * cos(Y[0] - Y[1])) - g_m * g_g * g_l * 
           (cos(Y[1]) + 2. * cos(Y[0]));
    while(dt < dt_max){
        for(int i = 0; i < m; i++){
            RK4(t, Y, dt, neq, dYdt);
            t += dt;
        }
        en = g_m * g_l * g_l * (Y[2] * Y[2] + 0.5 * Y[3] * Y[3] + 
             Y[2] * Y[3] * cos(Y[0] - Y[1])) - g_m * g_g * g_l * 
             (cos(Y[1]) + 2. * cos(Y[0]));
        file1 << dt << "   " << fabs((en_i - en) / en_i) << endl;
        t = 0.;
        Y[0] = 3.; Y[1] = -2.5; Y[2] = 1.; Y[3] = 0.;
        cout << dt <<endl;
        dt += ((dt_max - dt_min) / n);
    }
    file1.close();
#endif

#if STAGE == 3
    ofstream file1("moto_pendolo.txt");  
    //Initial conditions 
    Y[0] = 0.1; Y[1] = 0.1; Y[2] = 0.0; Y[3] = 0.;
    //Initial condition for chaotic motion
    Y1[0] = 3.0; Y1[1] = -2.5 ; Y1[2] = 1.; Y1[3] = 0.;
    double dt;      
    dt = 0.005;
    for(double t = 0.; t < 300.; t += dt){        
        RK4(t, Y, dt, neq, dYdt);
        RK4(t, Y1, dt, neq, dYdt);        
        file1 << t << "   " << Y[0] << "   " << Y[1] << "   " << Y[2] 
              << "   " << Y[3]  << "   " << Y1[0] << "   " << Y1[1] 
              << "   " << Y1[2] << "   " << Y1[3] << endl;
    }     
file1.close();
#endif

#if STAGE == 4
    double err = 1.e-8; 
    double x = 0.01, x_end, dx, parm;
    int N = 100, count = 0, n_points = 200;
    #if PARAM == 1
        //Lyapunov exponent for the angle of the first mass
        ofstream file1("esp_lyapunov_theta1.txt"); 
        x_end = M_PI;
        parm = 0;
        dx = x_end / n_points;
    #endif
    #if PARAM == 2
        //Lyapunov exponent for the angle of the first mass
        ofstream file1("esp_lyapunov_theta2.txt"); 
        x_end = M_PI;
        parm = 1;
        dx = x_end / n_points;
    #endif
    #if PARAM == 3
        //Lyapunov exponent for the angle of the first mass
        ofstream file1("esp_lyapunov_omega1.txt"); 
        x_end = 2 * M_PI;
        parm = 2;
        dx = x_end / n_points;
    #endif
    #if PARAM == 4
        //Lyapunov exponent for the angle of the first mass
        ofstream file1("esp_lyapunov_omega2.txt"); 
        x_end = 2 * M_PI;
        parm = 3;
        dx = x_end / n_points;
    #endif
    
    for(int i = 0; i < neq; i++){
        if(i == parm) Y[i] = x;        
        else Y[i] = 0.;
        
        if(i < 2) Y1[i] = Y[i] + err;        
        else Y1[i] = Y[i];
        
    }
    double dt = 0.005, t_end = 1000.;
    double esp_lyap1, esp_lyap2;
    double x1_0, x2_0, y1_0, y2_0, x1_1, x2_1, y1_1, y2_1;
    double d0_1 = 0., d0_2 = 0., d1_1 = 0., d1_2 = 0.;
    int n, n_step;
    t = 0.;
    n_step = t_end / dt;
    
    while(x < x_end){
        for(int i = 0; i < N; i++){
            //Cartesian coordinates for initial conditions
            x1_0 = g_l * sin(Y[0]);
            x2_0 = x1_0 + g_l * sin(Y[1]);
            y1_0 = -g_l * cos(Y[0]);
            y2_0 = y1_0 - g_l * cos(Y[1]);

            x1_1 = g_l * sin(Y1[0]);
            x2_1 = x1_1 + g_l * sin(Y1[1]);
            y1_1 = -g_l * cos(Y1[0]);
            y2_1 = y1_1 - g_l * cos(Y1[1]);
            
            //Euclidian distances for initial coonditions
            d0_1 += sqrt((x1_0 - x1_1) * (x1_0 - x1_1) + (y1_0 - y1_1) * (y1_0 - y1_1));
            d0_2 += sqrt((x2_0 - x2_1) * (x2_0 - x2_1) + (y2_0 - y2_1) * (y2_0 - y2_1));

            for(int k = 0; k < n_step; k++){        
                RK4(t, Y, dt, neq, dYdt);
                RK4(t, Y1, dt, neq, dYdt);
                t += dt;        
            }

            //Cartesian coordinates
            x1_0 = g_l * sin(Y[0]);
            x2_0 = x1_0 + g_l * sin(Y[1]);
            y1_0 = -g_l * cos(Y[0]);
            y2_0 = y1_0 - g_l * cos(Y[1]);

            x1_1 = g_l * sin(Y1[0]);
            x2_1 = x1_1 + g_l * sin(Y1[1]);
            y1_1 = -g_l * cos(Y1[0]);
            y2_1 = y1_1 - g_l * cos(Y1[1]);
            
            //Euclidian distances
            d1_1 += sqrt((x1_0 - x1_1) * (x1_0 - x1_1) + (y1_0 - y1_1) * (y1_0 - y1_1));
            d1_2 += sqrt((x2_0 - x2_1) * (x2_0 - x2_1) + (y2_0 - y2_1) * (y2_0 - y2_1));
            
            //Making the conditions close again
            for(int i = 0; i < neq; i++){     
                
                if(i < 2) Y1[i] = Y[i] + err;        
                else Y1[i] = Y[i];
            }
            t = 0.;
          
        } 
        count++;
        cout << count << endl;
        esp_lyap1 = log(fabs(d1_1 / d0_1)) / t_end;
        esp_lyap2 = log(fabs(d1_2 / d0_2)) / t_end;        
        file1 << x << "  " << esp_lyap1 << "   " << esp_lyap2 << endl;

        x += dx;
        //Initial conditions
        for(int i = 0; i < neq; i++){
            if(i == parm) Y[i] = x;        
            else Y[i] = 0.;
            
            if(i < 2) Y1[i] = Y[i] + err;        
            else Y1[i] = Y[i];            
        }
    }
    file1.close();
#endif

#if STAGE == 5
    ofstream file1("flip1.txt");
    ofstream file2("flip2.txt"); 
    double dt = 0.005, t_stop = 1000., tol = 1.e-7;
    int stop, dati = 720;
    int flip1 = 0, flip2 = 0;
    double theta1_old = 0., theta2_old = 0.;

    stop = t_stop / dt; 

    for(int i = 0; i < dati; i++){
        for(int j = 0; j < dati; j++){
            Y[0] = i * (2 * M_PI / dati) - M_PI;
            Y[1] = j * (2 * M_PI / dati) - M_PI;
            Y[2] = 0.; Y[3] = 0.;
            t = 0.;

            for(int m = 0; m < stop; m++){
                RK4(t, Y, dt, neq, dYdt);
                if(flip1 == 0 && sin(Y[0]) * sin(theta1_old) < 0. && t > dt && cos (Y [0]) < -1. + tol){
                    file1 << t << "   " << i * (2 * M_PI / dati) - M_PI << "   " << j * (2 * M_PI / dati) - M_PI << endl;
                    flip1 = 1;                    
                }
                if(flip2 == 0 && sin(Y[1]) * sin(theta2_old) < 0. && t > dt && cos (Y [1]) < -1. + tol){
                    file2 << t << "   " << i * (2 * M_PI / dati) - M_PI << "   " << j * (2 * M_PI / dati) - M_PI << endl;
                    flip2 = 1;                    
                }
                if(flip1 == 0) theta1_old = Y[0];
                if(flip2 == 0) theta2_old = Y[1];
                
                if(flip1 == 1 && flip2 == 1) break; 

                t += dt;

            }

            if(flip1 == 0){
                file1 << -1 << "   " << i * (2 * M_PI / dati) - M_PI << "   " << j * (2 * M_PI / dati) - M_PI << endl;
            }
            if(flip2 == 0){
                file2 << -1 << "   " << i * (2 * M_PI / dati) - M_PI << "   " << j * (2 * M_PI / dati) - M_PI << endl;
            }
            t = 0.;
            flip1 = 0;
            flip2 = 0;
            }
            cout << i << endl;
            file1 << endl;
            file2 << endl;
    }    
    file1.close();
    file2.close();
#endif


}

//=========================RK2================================

void RK2(double t, double *Y, double dt, int neq, 
         void (*RHS_Func)(double, double*, double*)){
    double Y1[neq], k1[neq], k2[neq];
    dYdt(t, Y, k1);
    for(int i = 0; i < neq ; i++){
        Y1[i] = Y[i] + 0.5 * dt * k1[i];
    }
    dYdt(t + 0.5 * dt, Y1, k2);
    for(int i = 0; i < neq ; i++){
        Y[i] += dt * k2[i];
    }
}

//========================RK4==========================

void RK4(double t, double *Y, double dt, int neq, 
         void (*RHS_Func)(double, double*, double*)){
    double Y1[neq], Y2[neq], Y3[neq];
    double k1[neq], k2[neq], k3[neq], k4[neq];
    dYdt(t, Y, k1);
    for(int i = 0; i < neq; i++){
        Y1[i] = Y[i] + 0.5 * dt * k1[i];
    }
    dYdt(t + 0.5 * dt, Y1, k2);
    for(int i = 0; i < neq ; i++){
        Y2[i] = Y[i] + 0.5 * dt * k2[i];
    }
    dYdt(t + 0.5 * dt, Y2, k3);
    for(int i = 0; i < neq ; i++){
        Y3[i] = Y[i] + dt * k3[i];
    }
    dYdt(t + dt, Y3, k4);
    for(int i = 0; i < neq ; i++){
        Y[i] +=  dt / 6. * (k1[i] + 2 * k2[i] + 2 *k3[i] + k4[i]);
    }

}

//========================dYdt================================

void dYdt(double t, double *Y, double *R){
    double theta1 = Y[0], theta2 = Y[1], v1 = Y[2], v2 = Y[3], den;
    den = g_l * (3. * g_m - g_m * cos((2. * theta1) - (2. * theta2)));
    R[0] = v1;
    R[1] = v2;
    R[2] = (- 3. * g_g * g_m * sin(theta1) - g_m * g_g * sin(theta1 - 2. * theta2) 
            - 2. * g_m * sin(theta1 - theta2) * (v2 * v2 * g_l + 
            v1 * v1 * g_l * cos(theta1 - theta2))) / den;
    R[3] = ((2. * sin(theta1 - theta2)) * (2. * g_m * g_l * v1 * v1 + 2.* g_g * g_m 
             * cos(theta1) + v2 * v2 * g_l * g_m * cos(theta1 - theta2))) / den;
    
}





