#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "thermodynamics.h"


// ch=1~4 are (1)  annih_bb.spec, (2) annih_ee.spec, (3) decay_bb.spec, (4) decay_ee.spec.
// T_or_Chi: integer, 0 for dT/dZ, else for dChie/dz
// mx: DM mass (GeV)
// z : redshift.

// Store the loaded data outside

static double KNS_mx[50], KNS_z[800], KNS_spec[50][800][2]; // mx,z,(0: dT/dz and 1: dChi/dz)
static int KNS_loaded;





// Please see the main as the example of usage.
/*
void main(){

   double zin,ansT,ansChi;
   int i;

   for ( i=0; i<800; i=i+1 ) {
      zin = KNS_z[i];
      ansT=Calc_dChidz_dTdz(2, 0, 0.1, zin);
      ansChi=Calc_dChidz_dTdz(2, 1, 0.1, zin);
      printf("%E %E %E\n",zin, ansT,ansChi);
   };


}
*/


// ----------------------------------------------------------------


void LoadKNS(int ch ){
    FILE *infile ;
    double  a, b, c;
    int i,j,k;
    char buf[10];
    double vmin;  // min DM mass
    double vmax;  // max DM mass
    double r;     // ratio


    if (ch==1) {
       infile = fopen("/work/CLASS_CMB/CLASS_CMB_mod/class/spec/annih_bb.spec", "r");
       vmin=5.0;
       };

    if (ch==2) {
       infile = fopen("/work/CLASS_CMB/CLASS_CMB_mod/class/spec/annih_ee.spec", "r");
       vmin=0.1;
       };


    if (ch==3) {
       infile = fopen("/work/CLASS_CMB/CLASS_CMB_mod/class/spec/decay_bb.spec", "r");
       vmin=11.0;
       };

    if (ch==4) {
       infile = fopen("/work/CLASS_CMB/CLASS_CMB_mod/class/spec/decay_ee.spec", "r");
       vmin=0.1;
       };


    vmax=1e3;
    r=pow(vmax/vmin, 1.0/49.0);
    KNS_mx[0]=vmin;

    fscanf(infile,"%*[^\n]");
    for ( i=0; i<50; i=i+1 ) {
       KNS_mx[i]=KNS_mx[0]*pow(r,i);
       for (j=0; j<800; j=j+1) {
          fscanf(infile, "%lf %lf %lf", &KNS_z[j], &KNS_spec[i][j][0], &KNS_spec[i][j][1]);
          //printf("%E %E %E %E\n", KNS_mx[i], KNS_z[j], KNS_spec[i][j][0], KNS_spec[i][j][1]);
       };
    };

    fclose(infile);
    KNS_loaded=1;
    //printf("Jui-Lin Kuo: here ... 1?\n");
    } ;

// ----------------------------------------------------------------



double Calc_dChidz_dTdz(int ch, int T_or_Chi, double mx, double z){
    int i0,i1,j0,j1,j;
    double x0,x1,y0,y1,z0,z1;
    double dChidz[2], dTdz[2]; // for two different masses
    double ans,r;

    //printf("Jui-Lin Kuo: ch = %d, mx = %f\n",ch, mx);
    ans=0.0;
    if (KNS_loaded==0) {LoadKNS(ch);};
    //printf("Jui-Lin Kuo: KNS_mxmin = %f\n",KNS_mx[0]);
    // NOTE: we have no data for the mass below KNS_mx[0].
    if (mx<KNS_mx[0]){
      printf("We have no data for the mass below %E GeV\n",KNS_mx[0]);
      return 0.0;
      };

    r=pow(KNS_mx[49]/KNS_mx[0], 1.0/49.0);

    // --- find index of mx ---
    if (mx>=KNS_mx[49]){
      i0=48;
      i1=49;
      }
    else {
      i0=log10(mx/KNS_mx[0])/log10(r);
      i1=i0+1;
      };

    //printf("%d %d\n",i0, i1);


    // --- find index of z ---
    if (z<KNS_z[0]){
      j0=0;
      j1=1;
      }
    else if (z>=KNS_z[799]){
      j0=798;
      j1=799;
      }
    else {
      for (j=0; j<799; j=j+1) {
         if (KNS_z[j]<z && KNS_z[j+1]>=z) {
           j0=j;
           j1=j+1;
           };
      };
    };



    x0=KNS_mx[i0];
    x1=KNS_mx[i1];

    //printf("%E %E\n",x0, x1);

    // --- for mass i0 ---
    y0=KNS_z[j0];
    y1=KNS_z[j1];

    if (T_or_Chi==0) {
      z0=KNS_spec[i0][j0][0];
      z1=KNS_spec[i0][j1][0];
      dTdz[0]=two_pts_interpolation(1, y0, y1, z0, z1, z);

//printf("y0,y1,z= %E %E %E\n",y0,y1,z);
//printf("nraws= %d\n",i0*800+j1+1);
//printf("z0,z1,result %E %E %E\n", z0, z1,dTdz[0]);
//exit(1);

      }
    else{
      z0=KNS_spec[i0][j0][1] ;
      z1=KNS_spec[i0][j1][1];
      dChidz[0]=two_pts_interpolation(1, y0, y1, z0, z1, z);
      };

    // --- for mass i1 and interpolation---
    if (T_or_Chi==0) {
      z0=KNS_spec[i1][j0][0] ;
      z1=KNS_spec[i1][j1][0];
      dTdz[1]=two_pts_interpolation(1, y0, y1, z0, z1, z);
      ans=two_pts_interpolation(1, x0, x1, dTdz[0], dTdz[1],mx);
      }
    else{
      z0=KNS_spec[i1][j0][1] ;
      z1=KNS_spec[i1][j1][1];
      dChidz[1]=two_pts_interpolation(1, y0, y1, z0, z1, z);
      ans=two_pts_interpolation(1, x0, x1, dChidz[0], dChidz[1],mx);
      };


    return (ans);
    } ;









// ----------------------------------------------------------------

double two_pts_interpolation(int use_log, double x0_, double x1_, double y0_, double y1_, double xin_ )
{
   double x0,x1,y0, y1, xin;
   double m,ans;



   if (use_log==1){
     x0=log10(x0_);
     x1=log10(x1_);
     if (y0_<=0.0) {y0_=1e-90;};
     if (y1_<=0.0) {y1_=1e-90;};
     y0=log10(y0_);
     y1=log10(y1_);
     xin=log10(xin_);
//printf("interpolation-1 = %E %E %E %E %E\n",x0_,x1_,y0_,y1_,xin_);
//printf("interpolation-2 = %E %E %E %E %E\n",x0,x1,y0,y1,xin);
     }
   else{
     x0=x0_;
     x1=x1_;
     y0=y0_;
     y1=y1_;
     xin=xin_;
     };

    if(x0==x1){
        ans=y0;
    }
    else{
        m=(y1-y0)/(x1-x0);
        ans=y0+m*(xin-x0);
    }
   if (use_log==1){ans=pow(10.0,ans); };

   return (ans);
   };
