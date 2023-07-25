#include <stdio.h>
#define _USE_MATH_DEFINES 
#include <math.h>
#include <stdlib.h>
#define R (8)
#define N (128)
#define NMESH_x (N)　//number of xmesh
#define NMESH_y (N) //number of vmesh
#define Ngreen (1 + NMESH_x / 2)
#define alpha (4.0)
#define eps (0.0) // pow(10, -10)
#define MIN(i, j) (((i) < (j)) ? (i) : (j))
#define MAX(i, j) (((i) > (j)) ? (i) : (j))
#define SQR(n) (n * n)
#define LF_eps (0.0e0)
#define eps_pp (1.0e-30)
double sgn(double x) // 符号関数を定義
{
    if (x < 0.0)
    {
        return -1.0;
    }
    else if (x > 0.0)
    {
        return 1.0;
    }
    else
    {
        return 0.0;
    }
}

double minmod(double a1, double a2) // minmod関数を定義(符号が異なると0、同じであれば小さいものを返す)
{
    return 0.5 * (copysign(1.0, a1) + copysign(1.0, a2)) * fmin(fabs(a1), fabs(a2));
}

double median(double a3, double a4, double a5) // 中央値を求める関数を定義
{
    return a3 + minmod(a4 - a3, a5 - a3);
}

double MIN3(double b1, double b2, double b3) // 3つの数の最小値を求める関数を定義
{
    double bmin;
    bmin = fmin(b1, fmin(b2, b3));
    return bmin;
}

double MAX3(double b4, double b5, double b6) // 3つの数の最大値を求める関数を定義
{

    double bmax;
    bmax = fmax(b4, fmax(b5, b6));
    return bmax;
}

double MIN4(double c1, double c2, double c3, double c4) // 4つの数の最小値を求める関数を定義
{
    double cmin;
    cmin = fmin(c1, fmin(c2, fmin(c3, c4)));
    return cmin;
}

double minmod4(double d1, double d2, double d3, double d4) // 4つの数のminmod関数を定義
{
    return 0.125 * (copysign(1.0, d1) + copysign(1.0, d2)) * fabs((copysign(1.0, d1) + copysign(1.0, d3)) * (copysign(1.0, d1) + copysign(1.0, d4))) * MIN4(fabs(d1), fabs(d2), fabs(d3), fabs(d4));
}

double CSL7(double fjm3, double fjm2, double fjm1, double fj, double fjp1, double fjp2, double fjp3,
            double *cCSL, double delta_tCSL, double deltaCSL, int i, int j)
{
    double nuCSL;
    double C_0, C_1, C_2, C_3, C_4, C_5, C_6;
    nuCSL = fabs(cCSL[j]) * delta_tCSL / deltaCSL; ///

    C_0 = (-1.0 * fjm3 / 140.0) + (5.0 * fjm2 / 84.0) - (101.0 * fjm1 / 420.0) + (319.0 * fj / 420.0) + (107.0 * fjp1 / 210.0) - (19.0 * fjp2 / 210.0) + (fjp3 / 105.0);
    C_1 = (fjm2 / 180.0) - (5.0 * fjm1 / 72.0) + (49.0 * fj / 72.0) - (49.0 * fjp1 / 72.0) + (5.0 * fjp2 / 72.0) - (fjp3 / 180.0);
    C_2 = (7.0 * fjm3 / 720.0) - (19.0 * fjm2 / 240.0) + (7.0 * fjm1 / 24.0) - (23.0 * fj / 72.0) + (fjp1 / 48.0) + (7.0 * fjp2 / 80.0) - (fjp3 / 90.0);
    C_3 = (-1.0 * fjm2 / 144.0) + (11.0 * fjm1 / 144.0) - (7.0 * fj / 36.0) + (7.0 * fjp1 / 36.0) - (11.0 * fjp2 / 144.0) + (fjp3 / 144.0);
    C_4 = (-1.0 * fjm3 / 360.0) + (fjm2 / 48.0) - (13.0 * fjm1 / 240.0) + (23.0 * fj / 360.0) - (fjp1 / 30.0) + (fjp2 / 240.0) + (fjp3 / 720.0);
    C_5 = (fjm2 / 720.0) - (fjm1 / 144.0) + (fj / 72.0) - (fjp1 / 72.0) + (fjp2 / 144.0) - (fjp3 / 720.0);
    C_6 = (fjm3 / 5040.0) - (fjm2 / 840.0) + (fjm1 / 336.0) - (fj / 252.0) + (fjp1 / 336.0) - (fjp2 / 840.0) + (fjp3 / 5040.0);

    double nuCSL_2, nuCSL_3, nuCSL_4, nuCSL_5, nuCSL_6;
    nuCSL_2 = nuCSL * nuCSL;
    nuCSL_3 = nuCSL_2 * nuCSL;
    nuCSL_4 = nuCSL_3 * nuCSL;
    nuCSL_5 = nuCSL_4 * nuCSL;
    nuCSL_6 = nuCSL_5 * nuCSL;
    double f_CSL7 = (C_0 + nuCSL * C_1 + nuCSL_2 * C_2 + nuCSL_3 * C_3 + nuCSL_4 * C_4 + nuCSL_5 * C_5 + nuCSL_6 * C_6);

    return f_CSL7;
}

double MP7(double fjm3, double fjm2, double fjm1, double fj, double fjp1, double fjp2, double fjp3,
           double *cMP, double delta_tMP, double deltaMP, int i, int j)
{
    double VL;
    double VMP;
    double D;
    double DM;
    double DP;
    double DM4P;
    double DM4M;
    double VUL;
    double VAV;
    double VMD;
    double VLC;
    double Vmin;
    double Vmax;
    double f_half;

    VL = CSL7(fjm3, fjm2, fjm1, fj, fjp1, fjp2, fjp3, cMP, delta_tMP, deltaMP, i, j);
    VMP = fj + minmod(fjp1 - fj, alpha * (fj - fjm1));
    if ((VL - fj) * (VL - VMP) < eps)
    {
        f_half = VL;
    }
    else
    {
        DM = fjm2 - 2.0 * fjm1 + fj;
        D = fjm1 - 2.0 * fj + fjp1;
        DP = fj - 2.0 * fjp1 + fjp2;
        DM4M = minmod4(4.0 * DM - D, 4.0 * D - DM, DM, D);
        DM4P = minmod4(4.0 * D - DP, 4.0 * DP - D, D, DP);
        VUL = fj + alpha * (fj - fjm1);
        VAV = 0.5 * (fj + fjp1);
        VMD = VAV - 0.5 * DM4P;
        VLC = fj + 0.5 * (fj - fjm1) + 4.0 * DM4M / 3.0;
        Vmin = fmax(MIN3(fj, fjp1, VMD), MIN3(fj, VUL, VLC));
        Vmax = fmin(MAX3(fj, fjp1, VMD), MAX3(fj, VUL, VLC));
        f_half = median(VL, Vmin, Vmax);
    }

    return f_half;
}

double PP(double fim3, double fim2, double fim1, double fi, double fip1, double fip2, double fip3,
          double *cPP, double delta_tPP, double deltaPP, int i, int j)
{

    double f_up_R;
    double f_bar_R;
    double U_up_p;
    double U_up_m_p1;
    double Up;
    double Um_p1;
    double theta_p_R;
    double theta_m_R;
    double theta_R;
    double f_PP_R;
    double nu1;
    f_bar_R = MP7(fim3, fim2, fim1, fi, fip1, fip2, fip3, cPP, delta_tPP, deltaPP, i, j);
    nu1 = fabs(cPP[j]) * delta_tPP / deltaPP;
    Up = fi - 2.0 * nu1 * f_bar_R;
    Um_p1 = fip1 + 2.0 * nu1 * f_bar_R;

    f_up_R = fi;

    theta_m_R = 1.0;
    theta_p_R = 1.0;

    if (Up < LF_eps)
    {
        U_up_p = fi - 2.0 * nu1 * f_up_R;
        theta_p_R = (LF_eps - U_up_p) / (Up - U_up_p + eps_pp);
    }

    if (Um_p1 < 0.0)
    {
        U_up_m_p1 = fip1 + 2.0 * nu1 * f_up_R;
        theta_m_R = (LF_eps - U_up_m_p1) / (Um_p1 - U_up_m_p1 + eps_pp);
    }

    theta_R = MIN(theta_p_R, theta_m_R);

    f_PP_R = theta_R * f_bar_R + (1.0 - theta_R) * f_up_R;

    return f_PP_R;
}

void sweep_x(double *f_star1, double *fold,
             double *csweep_x, double delta_tsweepx, double deltasweep_x, int i, int j)
{
    int im4, im3, im2, im1, ip1, ip2, ip3, ip4;
    double fr, fl;
    double fr_star1, fl_star1, fr_star2, fl_star2;

    im4 = (i - 4 + N) % N;
    im3 = (i - 3 + N) % N;
    im2 = (i - 2 + N) % N;
    im1 = (i - 1 + N) % N;
    ip1 = (i + 1) % N;
    ip2 = (i + 2) % N;
    ip3 = (i + 3) % N;
    ip4 = (i + 4) % N;
    if (0.0 <= csweep_x[j])
    {

        fl = MP7(fold[j + NMESH_y * im4], fold[j + NMESH_y * im3], fold[j + NMESH_y * im2], fold[j + NMESH_y * im1],
                 fold[j + NMESH_y * i], fold[j + NMESH_y * ip1], fold[j + NMESH_y * ip2], csweep_x, delta_tsweepx, deltasweep_x, i, j);
        fr = MP7(fold[j + NMESH_y * im3], fold[j + NMESH_y * im2], fold[j + NMESH_y * im1], fold[j + NMESH_y * i],
                 fold[j + NMESH_y * ip1], fold[j + NMESH_y * ip2], fold[j + NMESH_y * ip3], csweep_x, delta_tsweepx, deltasweep_x, i, j);

        f_star1[j + NMESH_y * i] = fold[j + NMESH_y * i] - ((csweep_x[j]) * delta_tsweepx / deltasweep_x) * (fr - fl);

        if (f_star1[j + NMESH_y * i] < 0.0)
        {
            fl = PP(fold[j + NMESH_y * im4], fold[j + NMESH_y * im3], fold[j + NMESH_y * im2], fold[j + NMESH_y * im1],
                    fold[j + NMESH_y * i], fold[j + NMESH_y * ip1], fold[j + NMESH_y * ip2], csweep_x, delta_tsweepx, deltasweep_x, i, j);
            fr = PP(fold[j + NMESH_y * im3], fold[j + NMESH_y * im2], fold[j + NMESH_y * im1], fold[j + NMESH_y * i],
                    fold[j + NMESH_y * ip1], fold[j + NMESH_y * ip2], fold[j + NMESH_y * ip3], csweep_x, delta_tsweepx, deltasweep_x, i, j);
            f_star1[j + NMESH_y * i] = fold[j + NMESH_y * i] - ((csweep_x[j]) * delta_tsweepx / deltasweep_x) * (fr - fl);
        }
    }
    else
    {
        fr = MP7(fold[j + NMESH_y * ip4], fold[j + NMESH_y * ip3], fold[j + NMESH_y * ip2], fold[j + NMESH_y * ip1],
                 fold[j + NMESH_y * i], fold[j + NMESH_y * im1], fold[j + NMESH_y * im2], csweep_x, delta_tsweepx, deltasweep_x, i, j);
        fl = MP7(fold[j + NMESH_y * ip3], fold[j + NMESH_y * ip2], fold[j + NMESH_y * ip1], fold[j + NMESH_y * i],
                 fold[j + NMESH_y * im1], fold[j + NMESH_y * im2], fold[j + NMESH_y * im3], csweep_x, delta_tsweepx, deltasweep_x, i, j);
        f_star1[j + NMESH_y * i] = fold[j + NMESH_y * i] - ((csweep_x[j]) * delta_tsweepx / deltasweep_x) * (fr - fl);

        if (f_star1[j + NMESH_y * i] < 0.0)
        {
            fr = PP(fold[j + NMESH_y * ip4], fold[j + NMESH_y * ip3], fold[j + NMESH_y * ip2], fold[j + NMESH_y * ip1],
                    fold[j + NMESH_y * i], fold[j + NMESH_y * im1], fold[j + NMESH_y * im2], csweep_x, delta_tsweepx, deltasweep_x, i, j);
            fl = PP(fold[j + NMESH_y * ip3], fold[j + NMESH_y * ip2], fold[j + NMESH_y * ip1], fold[j + NMESH_y * i],
                    fold[j + NMESH_y * im1], fold[j + NMESH_y * im2], fold[j + NMESH_y * im3], csweep_x, delta_tsweepx, deltasweep_x, i, j);
            f_star1[j + NMESH_y * i] = fold[j + NMESH_y * i] - ((csweep_x[j]) * delta_tsweepx / deltasweep_x) * (fr - fl);
        }
    }
}

void sweep_y(double *f_star2, double *f_star1,
             double *csweep_y, double delta_tsweepy, double deltasweep_y, int i, int j)
{

    int jm4, jm3, jm2, jm1, jp1, jp2, jp3, jp4;
    double fr, fl;
    double fr_star1, fl_star1, fr_star2, fl_star2;

    jm4 = (j - 4);
    if (jm4 < 0)
        jm4 = 0;
    jm3 = (j - 3);
    if (jm3 < 0)
        jm3 = 0;
    jm2 = (j - 2);
    if (jm2 < 0)
        jm2 = 0;
    jm1 = (j - 1);
    if (jm1 < 0)
        jm1 = 0;
    jp1 = (j + 1);
    if (jp1 >= NMESH_y)
        jp1 = NMESH_y - 1;
    jp2 = (j + 2);
    if (jp2 >= NMESH_y)
        jp2 = NMESH_y - 1;
    jp3 = (j + 3);
    if (jp3 >= NMESH_y)
        jp3 = NMESH_y - 1;
    jp4 = (j + 4);
    if (jp4 >= NMESH_y)
        jp4 = NMESH_y - 1;
    if (0.0 <= csweep_y[i])
    {

        fl_star1 = MP7(f_star1[jm4 + NMESH_y * i], f_star1[jm3 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], f_star1[jm1 + NMESH_y * i],
                       f_star1[j + NMESH_y * i], f_star1[jp1 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);
        fr_star1 = MP7(f_star1[jm3 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], f_star1[jm1 + NMESH_y * i], f_star1[j + NMESH_y * i],
                       f_star1[jp1 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], f_star1[jp3 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);

        f_star2[j + NMESH_y * i] = f_star1[j + NMESH_y * i] - ((csweep_y[i]) * delta_tsweepy / deltasweep_y) * (fr_star1 - fl_star1);

        if (f_star2[j + NMESH_y * i] < 0.0)
        {
            fl_star1 = PP(f_star1[jm4 + NMESH_y * i], f_star1[jm3 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], f_star1[jm1 + NMESH_y * i],
                          f_star1[j + NMESH_y * i], f_star1[jp1 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);
            fr_star1 = PP(f_star1[jm3 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], f_star1[jm1 + NMESH_y * i], f_star1[j + NMESH_y * i],
                          f_star1[jp1 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], f_star1[jp3 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);
            f_star2[j + NMESH_y * i] = f_star1[j + NMESH_y * i] - ((csweep_y[i]) * delta_tsweepy / deltasweep_y) * (fr_star1 - fl_star1);
        }
    }
    else
    {
        fr_star1 = MP7(f_star1[jp4 + NMESH_y * i], f_star1[jp3 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], f_star1[jp1 + NMESH_y * i],
                       f_star1[j + NMESH_y * i], f_star1[jm1 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);
        fl_star1 = MP7(f_star1[jp3 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], f_star1[jp1 + NMESH_y * i], f_star1[j + NMESH_y * i],
                       f_star1[jm1 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], f_star1[jm3 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);

        f_star2[j + NMESH_y * i] = f_star1[j + NMESH_y * i] - ((csweep_y[i]) * delta_tsweepy / deltasweep_y) * (fr_star1 - fl_star1);

        if (f_star2[j + NMESH_y * i] < 0.0)
        {
            fr_star1 = PP(f_star1[jp4 + NMESH_y * i], f_star1[jp3 + NMESH_y * i], f_star1[jp2 + NMESH_y * i],
                          f_star1[jp1 + NMESH_y * i], f_star1[j + NMESH_y * i], f_star1[jm1 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);
            fl_star1 = PP(f_star1[jp3 + NMESH_y * i], f_star1[jp2 + NMESH_y * i], f_star1[jp1 + NMESH_y * i],
                          f_star1[j + NMESH_y * i], f_star1[jm1 + NMESH_y * i], f_star1[jm2 + NMESH_y * i], f_star1[jm3 + NMESH_y * i], csweep_y, delta_tsweepy, deltasweep_y, i, j);
            f_star2[j + NMESH_y * i] = f_star1[j + NMESH_y * i] - ((csweep_y[i]) * delta_tsweepy / deltasweep_y) * (fr_star1 - fl_star1);
        }
    }
}

void calc_distribution_function(double *DF_, double rho_bar_, double kwave_, double A_, double sigma_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        double x_ = i / (double)NMESH_x;
        for (int j = 0; j < NMESH_y; j++)
        {
            double v_ = 2.0 * j / (double)NMESH_y - 1.0;
            DF_[j + NMESH_y * i] = rho_bar_ * (1.0 + A_ * sin(kwave_ * x_)) * exp(-SQR(v_) / (2.0 * SQR(sigma_))) / sqrt(2.0 * M_PI * SQR(sigma_));
        }
    }
}

void calc_dens_field_1(double *rho_, double rho_bar_, double kwave_, double A_, double delta_x_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        double x_ = i / (double)NMESH_x;
        rho_[i] = rho_bar_ * (1.0 + A_ * cos(kwave_ * x_));
        // rho_[i] = rho_bar_ * (1.0 + A_ * sin(kwave_ * x_));
    }
    for (int i = 0; i < NMESH_x; i++)
    {
        rho_[i] *= delta_x_;
    }
}
void calc_dens_field_2(double *rho_, double *DF_, double delta_x_, double delta_y_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        rho_[i] = 0.0;
        for (int j = 0; j < NMESH_y; j++)
        {
            rho_[i] += DF_[j + NMESH_y * i] * delta_y_;
        }
    }
    for (int i = 0; i < NMESH_x; i++)
    {
        rho_[i] *= delta_x_;
    }
}

void calc_green_2nd(double *green_, double delta_x_, double G_)
{
    green_[0] = 0.0;
    for (int k = 1; k < Ngreen; k++)
    {
        double sinc;
        sinc = sin(M_PI * (double)k * delta_x_) / delta_x_;
        green_[k] = -G_ * M_PI / (SQR(sinc));
    }
}
void calc_green_4th(double *green_, double delta_x_, double G_)
{
    green_[0] = 0.0;
    for (int k = 1; k < Ngreen; k++)
    {
        green_[k] = 12.0 * M_PI * G_ * SQR(delta_x_) / (SQR(sin(2.0 * M_PI * (double)k * delta_x_)) - 16.0 * SQR(sin(M_PI * (double)k * delta_x_)));
    }
}
void calc_green_6th(double *green_, double delta_x_, double G_)
{
    green_[0] = 0.0;
    for (int k = 1; k < Ngreen; k++)
    {
        green_[k] = -180.0 * G_ * M_PI * SQR(delta_x_) / (2.0 * SQR(sin(3.0 * M_PI * (double)k * delta_x_)) - 27.0 * SQR(sin(2.0 * M_PI * (double)k * delta_x_)) + 270.0 * SQR(sin(M_PI * (double)k * delta_x_)));
    }
}

void calc_rho_DFT(double *rho_hat_real_, double *rho_hat_image_, double *rho_, double delta_x_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        rho_hat_real_[i] = 0.0;
        rho_hat_image_[i] = 0.0;
        for (int j = 0; j < NMESH_x; j++)
        {
            rho_hat_real_[i] += ((rho_[j]) * cos(2.0 * M_PI * ((double)i * (double)j) * delta_x_));
            rho_hat_image_[i] += -((rho_[j]) * sin(2.0 * M_PI * ((double)i * (double)j) * delta_x_));
        }
    }
}

void calc_phi_DFT(double *rho_hat_real_, double *rho_hat_image_, double *green_)
{
    for (int k = 0; k < NMESH_x / 2; k++)
    {
        rho_hat_real_[k] *= green_[k];
        rho_hat_image_[k] *= green_[k];
        rho_hat_real_[NMESH_x / 2 + k] *= green_[NMESH_x / 2 - k];
        rho_hat_image_[NMESH_x / 2 + k] *= green_[NMESH_x / 2 - k];
    }
}

void calc_phi_IDFT(double *phi_real_, double *phi_image_, double *rho_hat_real_, double *rho_hat_image_, double delta_x_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        phi_real_[i] = 0.0;
        phi_image_[i] = 0.0;
        for (int j = 0; j < NMESH_y; j++)
        {
            phi_real_[i] += rho_hat_real_[j] * cos(2.0 * M_PI * ((double)i * (double)j) * delta_x_) -
                            rho_hat_image_[j] * sin(2.0 * M_PI * ((double)i * (double)j) * delta_x_); //?
            phi_image_[i] += rho_hat_real_[j] * sin(2.0 * M_PI * ((double)i * (double)j) * delta_x_) +
                             rho_hat_image_[j] * cos(2.0 * M_PI * ((double)i * (double)j) * delta_x_);
        }
    }
}

void calc_pot(double *pot_, double *phi_real_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        pot_[i] = phi_real_[i];
    }
}

void calc_force_2nd(double *acc_2nd_, double *pot_, double delta_x_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        int im1, ip1;

        im1 = (i - 1 + N) % N;
        ip1 = (i + 1) % N;

        acc_2nd_[i] = -0.5 * (pot_[ip1] - pot_[im1]) / delta_x_;
    }
}
void calc_force_4th(double *acc_4th_, double *pot_, double delta_x_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        int im2, im1, ip1, ip2;

        im2 = (i - 2 + N) % N;
        im1 = (i - 1 + N) % N;
        ip1 = (i + 1) % N;
        ip2 = (i + 2) % N;

        acc_4th_[i] = ((pot_[ip2] - pot_[im2]) - 8.0 * (pot_[ip1] - pot_[im1])) / (12.0 * delta_x_);
    }
}
void calc_force_6th(double *acc_6th_, double *pot_, double delta_x_)
{
    for (int i = 0; i < NMESH_x; i++)
    {
        int im3, im2, im1, ip1, ip2, ip3;

        im3 = (i - 3 + N) % N;
        im2 = (i - 2 + N) % N;
        im1 = (i - 1 + N) % N;
        ip1 = (i + 1) % N;
        ip2 = (i + 2) % N;
        ip3 = (i + 3) % N;

        acc_6th_[i] = -((pot_[ip3] - pot_[im3]) - 9.0 * (pot_[ip2] - pot_[im2]) + 45.0 * (pot_[ip1] - pot_[im1])) / (60.0 * delta_x_);
    }
}

void calc_max(double max, double *acc_)
{

    static double num[NMESH_x];
    max = fabs(acc_[0]);

    for (int i = 0; i < NMESH_x; i++)
    {
        num[i] = fabs(acc_[i]);
        if (num[i] > max)
        {
            max = num[i];
        }
    }
}

void initialization1(double *v_ini, double kwave_, double rho_bar_, double A_, double sigma_)
{
    double x_, v_;
    double sigma2_ = sigma_ * sigma_;
    for (int i = 0; i < NMESH_x; i++)
    {
        x_ = i / (double)NMESH_x;
        for (int j = 0; j < NMESH_y; j++)
        {
            v_ = 2.0 * j / (double)NMESH_y - 1.0;

            v_ini[j + NMESH_y * i] = (rho_bar_ * (1.0 + A_ * cos(kwave_ * x_)) / sqrt(2.0 * M_PI * sigma2_)) * exp(-SQR(v_) / sigma2_);
        }
    }
}

void initialization2(double *v_ini, double L_, double T_, double rho_bar_, double A_)
{
    double x_, v_, rho, sigma_, sigma2, k_;
    sigma_ = L_ / (2.0 * sqrt(4.0 * M_PI) * T_);
    sigma2 = sigma_ * sigma_;
    k_ = 4.0 * M_PI / L_;
    for (int i = 0; i < NMESH_x; i++)
    {
        x_ = i / (double)NMESH_x;
        rho = rho_bar_ * (1.0 + A_ * cos(k_ * x_));
        for (int j = 0; j < NMESH_y; j++)
        {
            v_ = 2.0 * j / (double)NMESH_y - 1.0;

            v_ini[j + NMESH_y * i] = (rho / sqrt(2.0 * M_PI * sigma2)) * exp(-v_ * v_ / sigma2);
        }
    }
}

void calc_delta_t(double delta_t_, double *c_x_, double *c_y_, double delta_x_, double delta_y_, double nu_)
{
    double delta_t_x, delta_t_y, V_x_max, V_y_max;

    V_x_max = fabs(c_x_[0]);
    for (int i = 0; i < NMESH_y; i++)
    {
        if (fabs(c_x_[i]) > V_x_max)
        {
            V_x_max = fabs(c_x_[i]);
        }
    }

    V_y_max = fabs(c_y_[0]);
    for (int i = 0; i < NMESH_x; i++)
    {
        if (fabs(c_y_[i]) > V_y_max)
        {
            V_y_max = fabs(c_y_[i]);
        }
    }

    delta_t_x = nu_ * (delta_x_ / V_x_max);
    delta_t_y = nu_ * (delta_y_ / V_y_max);
    delta_t_ = fmin(delta_t_x, delta_t_y);
    // printf("%f %f\n", delta_t_x, delta_t_y);

    // printf("%f\n", delta_t_);
}

int main()
{
    FILE *fp_sin;
    char *fname_sin[256];
    FILE *fp_sin4;
    char *fname_sin4[256];
    FILE *fp_vlasov;
    char *fname_vlasov[256];

    static double v_vlasov[NMESH_x * NMESH_y];
    static double v_vlasov_star1[NMESH_x * NMESH_y];
    static double v_vlasov_star2[NMESH_x * NMESH_y];
    static double vnew_vlasov[NMESH_x * NMESH_y];

    double delta_x;
    double delta_y;
    double nu = (1.0 / ((1.0 + alpha)));
    double c_x[NMESH_y];
    double c_y[NMESH_y];
    double delta_t;
    double A;
    double L;
    double rho_bar;
    double T;
    double t;
    int count;
    double k_J, kwave, x, v, sigma;
    int wave_num;
    double G;
    delta_x = (1.0 / (double)(NMESH_x));
    delta_y = (2.0 / (double)(NMESH_y));
    G = 1.0;
    A = 0.01;
    L = 1.0;
    rho_bar = 1.0;
    T = 1.0 / sqrt(G * rho_bar);
    wave_num = 2;
    kwave = 2.0 * M_PI * (double)wave_num / L;
    k_J = kwave / 0.5;
    sigma = sqrt(4.0 * M_PI * G * rho_bar) / k_J;

    // calc_pot&force

    static double DF[NMESH_x * NMESH_y];
    // calc_distribution_function(DF, rho_bar, kwave, A, sigma);

    static double rho[NMESH_x];
    calc_dens_field_1(rho, rho_bar, kwave, A, delta_x);
    // calc_dens_field_2(rho, DF, delta_x, delta_y);

    static double green[Ngreen];
    // calc_green_2nd(green, delta_x, G);
    // calc_green_4th(green, delta_x, G);
    calc_green_6th(green, delta_x, G);

    static double rho_hat_real[NMESH_x], rho_hat_image[NMESH_x];
    calc_rho_DFT(rho_hat_real, rho_hat_image, rho, delta_x);

    calc_phi_DFT(rho_hat_real, rho_hat_image, green);

    static double phi_real[NMESH_x], phi_image[NMESH_x];
    calc_phi_IDFT(phi_real, phi_image, rho_hat_real, rho_hat_image, delta_x);

    static double pot[NMESH_x];
    calc_pot(pot, phi_real);

    static double acc[NMESH_x];
    // calc_force_2nd(acc, pot, delta_x);
    // calc_force_4th(acc, pot, delta_x);
    calc_force_6th(acc, pot, delta_x);

    static double acc_analy[NMESH_x];
    for (int i = 0; i < NMESH_x; i++)
    {
        x = (double)i / (double)NMESH_x;

        acc_analy[i] = -(4.0 * M_PI * G * rho_bar * A) * sin(kwave * x) / kwave;
    }

    // calc_e_acc
    double e_acc_6th;
    double acc_numerator_6th, acc_denominator;
    acc_numerator_6th = 0.0;
    acc_denominator = 0.0;
    for (int i = 0; i < NMESH_x; i++)
    {
        acc_numerator_6th += fabs(acc[i] - acc_analy[i]);
        acc_denominator += fabs(acc_analy[i]);
    }

    e_acc_6th = acc_numerator_6th / acc_denominator;

    for (int i = 0; i < NMESH_y; i++)
    {
        c_x[i] = 2.0 * (double)i / (double)NMESH_y - 1.0;
    }
    for (int i = 0; i < NMESH_x; i++)
    {
        c_y[i] = acc[i];
    }

    double delta_t_x, delta_t_y, V_x_max, V_y_max;
    for (int i = 0; i < NMESH_x; i++)
    {
        V_x_max = fabs(c_x[0]);
        if (fabs(c_x[i]) > V_x_max)
        {
            V_x_max = fabs(c_x[i]);
        }
    }

    for (int i = 0; i < NMESH_x; i++)
    {
        V_y_max = fabs(c_y[0]);
        if (fabs(c_y[i]) > V_y_max)
        {
            V_y_max = fabs(c_y[i]);
        }
    }

    delta_t_x = nu * (delta_x / V_x_max);
    delta_t_y = nu * (delta_y / V_y_max);
    delta_t = fmin(delta_t_x, delta_t_y);

    initialization1(v_vlasov, kwave, rho_bar, A, sigma);

    t = 0.0;
    count = 0.0;

    while (t < 4.0 * T)
    {
        if (count % 10 == 0)
        {

            double v;
            sprintf(fname_vlasov, "SL-MPP7_%d.dat", count);
            fp_vlasov = fopen(("%s", fname_vlasov), "w");
            for (int i = 0; i < NMESH_x; i++)
            {
                x = i / (double)NMESH_x;
                for (int j = 0; j < NMESH_y; j++)
                {
                    v = 2.0 * j / (double)NMESH_y - 1.0;
                    fprintf(fp_vlasov, "%f %f %f\n", 0.5 * delta_x + x, 0.5 * delta_y + v, v_vlasov[j + NMESH_y * i]);
                }
                fprintf(fp_vlasov, "\n");
            }
        }
        fclose(fp_vlasov);

        for (int i = 0; i < NMESH_x; i++)
        {

            for (int j = 0; j < NMESH_y; j++)
            {
                sweep_y(v_vlasov_star1, v_vlasov, c_y, 0.5 * delta_t, delta_y, i, j);
            }
        }

        for (int j = 0; j < NMESH_y; j++)
        {
            for (int i = 0; i < N; i++)
            {

                sweep_x(v_vlasov_star2, v_vlasov_star1, c_x, delta_t, delta_x, i, j);
            }
        }

        // calc_acc
        calc_dens_field_2(rho, v_vlasov_star2, delta_x, delta_y);
        calc_rho_DFT(rho_hat_real, rho_hat_image, rho, delta_x);
        calc_phi_DFT(rho_hat_real, rho_hat_image, green);
        calc_phi_IDFT(phi_real, phi_image, rho_hat_real, rho_hat_image, delta_x);
        calc_pot(pot, phi_real);
        calc_force_6th(c_y, pot, delta_x);
        for (int i = 0; i < NMESH_x; i++)
        {
            V_y_max = fabs(c_y[0]);
            if (fabs(c_y[i]) > V_y_max)
            {
                V_y_max = fabs(c_y[i]);
            }
        }
        delta_t_y = nu * (delta_y / V_y_max);
        delta_t = fmin(delta_t_x, delta_t_y);

        for (int i = 0; i < NMESH_x; i++)
        {
            for (int j = 0; j < NMESH_y; j++)
            {
                sweep_y(vnew_vlasov, v_vlasov_star2, c_y, 0.5 * delta_t, delta_y, i, j);
            }
        }
        for (int j = 0; j < NMESH_y; j++)
        {
            for (int i = 0; i < NMESH_x; i++)
            {

                v_vlasov[j + NMESH_y * i] = vnew_vlasov[j + NMESH_y * i];
            }
        }
        t += delta_t;
        count++;
    }
    printf("%d %f %f %f", count, t, delta_t, T);
}
