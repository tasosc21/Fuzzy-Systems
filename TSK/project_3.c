/************************************************************************
*    -FNN using FUNCOM-LSE.                                             *
*    -The system to be identified is a SISO one.                        *
*    -The data file is created within the program.                      *
*    -Each Gaussian membership function has TWO (2) fitting parameters. *
*    -Function to maximize: current - average of previous epoch.        *
*    -step size adaptation mechanism                                    *
*************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
/*------------------------------------------------------*/
#define rules 5           /* number of fuzzy rules */
#define g 0.5             /* overlapping factor */
#define iterat 200        /* number of training epochs */
#define dpparh 0.01       /* initial step size */
#define precision 0.0001 /* RMSE */
/*------------------------------------------------------*/
#define MY_PI 3.1415
#define steps_back 3
#define samples 361
#define low1 0.999999
#define up_con1 1.005
#define up_prem1 1.005
#define low_prem2 0.8
#define low_con2 0.8
#define up2 1.005
#define et 0.1
#define ksi 0.9
#define false 0
#define true 1
/*----------------------------------------------------------------------*/
float upper_in,lower_in,upper_out,lower_out;
float mini_in,maxi_in,mini_out,maxi_out;

float degful[samples][rules],fuzout[samples][rules],membshp[samples][rules];

float wprm[rules],dwprm[rules],wprevm[steps_back][rules];
float wprdev[rules],dwprdev[rules],wprevdev[steps_back][rules];
float w[2][rules];

float lagrdof[samples][rules],lagrfout[samples][rules],lagr_out[samples];
float lfm[rules],lfdev[rules],lem[rules],ledev[rules];
float iffp,iefp,ieep,rmse,dpp;

float inpat[samples],outpat[samples],output[samples];
FILE *dfiniter;
/*----------------------------------------------------------------------*/
/* Function declaration */

float sqr(float x);
void create_data(void);
void normal_val(void);
void denormal_val(int end);
void memb_value(int iter);
void deg_ful(void);
void find_con_w(int M);
void fuz_out(void);
void out_put(void);
void lagr_mult(void);
void aux_le_lf(void);
void aux_i_ef(void);
void perturb(float dep);
void cor_weights(void);
void prev_epoch(void);
void initialization(void);
void features(void);
/*-----------------------------------------------------------------------*/
int main()
{
  int counter=0;
  int end=false; /* if (end==true) the datumtel.dat will be written */
  int i,j,k,iter,finish,M;
  float counter1,dep,counter2,fi,mse,e0,e1,e2,e3;
  FILE *derror;
  FILE *dfuzout;
  FILE *f1;
  dfuzout=fopen("dfuzout.dat","w");
  dfiniter=fopen("dfiniter.dat","w");
  derror=fopen("derror.dat","w");
  create_data();
  normal_val();
  initialization();
  counter2=0.0;
  iter=0;
  e0=20;
  e1=10;
  e2=5;
  finish=false;
  dpp=dpparh;
  features();
  M=2*rules;
  while ((iter<iterat) && (!finish))
  {
    memb_value(iter);
    deg_ful();
    find_con_w(M);
    fuz_out();
    out_put();
    denormal_val(end);
    printf( "iter=%d\n",iter );
    counter2=counter2+1.0;
    if (iter==0) fprintf( dfiniter,"initial error:%f\n",rmse );
    iter=iter+1;
    e3=rmse;
    if (rmse<precision)
    {
      finish=true;
      counter++;
      printf( "precision was attained in the %d iteration\n",iter );
      fprintf( dfiniter,"precision was attained in the %d iteration\n",iter );
      fprintf( dfiniter,"iter:%d\n",iter );
    }
    else
    {
      lagr_mult();
      aux_le_lf();
      aux_i_ef();
      dep=-ksi*dpp*sqrt((double)ieep);
      if ( !ieep && !iefp )
      {
	    printf( "error!!!!!!!!!!!!   iee=0\n" );
	    scanf( "%c" );
      }
      perturb(dep);
      prev_epoch();
      cor_weights();
      if (e3>e2)
      {
	    if ((e2>e1) || (e1>e0)) dpp=low_prem2*dpp;
	    else
	    {
	      if ((e1>e2) && (e2>e3))
	      {
	        if (((e1-e2)/e1)>((e2-e3)/e2)) dpp=low1*dpp;
	        else dpp=up_prem1*dpp;
	      }
	    }
      }
    }
    e0=e1;
    e1=e2;
    e2=e3;
    printf( "new dpp=%f\n",dpp );
    fprintf( derror,"%d     %f     %f\n",iter,rmse,dpp );
  }
  if (iter>=iterat) printf( "There was no convergence pattern after %d iterations.\n",iterat );
  for (j=0;j<samples;j++)
  {
    for (i=0;i<rules;i++)
    {
      if (membshp[j][i]>0.01) fprintf( dfuzout,"%f     ",fuzout[j][i] );
      else fprintf( dfuzout,"%f     ",0.0 );
    }
    fprintf( dfuzout,"%f\n",inpat[j] );
  }
  printf( "Total number of iterations:%f\n",counter2 );
  counter1=0;

  end=true;             /*  to present the final mapping  */
  denormal_val(end);    /*  in the actual space           */


  /* Copy the final membership values to file "inal_memb.dat" */
  f1=fopen("final_memb.dat","w");
  for (k=0;k<samples;k++)
    for (j=0;j<rules;j++)
    {
      if (-0.5*sqr(wprm[j]-inpat[k])/sqr(wprdev[j])<-10.0) membshp[k][j]=0;
      else membshp[k][j]=exp(-0.5*sqr(wprm[j]-inpat[k])/sqr(wprdev[j]));
      if (j==(rules-1)) fprintf(f1,"%f\n",membshp[k][j]);
	  else fprintf(f1,"%f ",membshp[k][j]);
    }
  fclose(f1);


  for (j=0;j<rules;j++)
    fprintf( dfiniter,"wprm[%d]=%f     wprdev[%d]=%f\n",j,wprm[j],j,wprdev[j] );
  fprintf( dfiniter,"RMSE=%f\n",rmse );
  for (j=0;j<samples;j++) counter1=counter1+fabs(outpat[j]-output[j]);
  fprintf( dfiniter," Mean Error=%f\n",counter1/(samples+0.0) );
  printf( "mean error:%f\n",counter1/(samples+0.0) );
  fprintf( dfiniter,"Total number of iterations:%f\n",counter2 );
  fclose(dfuzout);
  fclose(dfiniter);
  fclose(derror);

  return(0);
}
/*-----------------------------------------------------------------------*/
float sqr(float x)
{
  x=x*x;
  return(x);
}
/*-----------------------------------------------------------------------*/
/* creation of input and output samples */

void create_data(void)
{
  int j,i;
  float step1,step2,x0,x1;
  for (i=0;i<samples;i++)
  {
    inpat[i]=i;
    outpat[i]=0.5*sin(2*MY_PI*i/360)+0.5*sqr(sin(2*MY_PI*(i+30)/360));
  }
}
/*----------------------------------------------------------------------*/
/*  normalization of the input-output samples [-1,1]                    */

void normal_val(void)
{
  int i,j;
  lower_in=-1.0;
  upper_in=1.0;
  lower_out=-1.0;
  upper_out=1.0;
  mini_in=inpat[0];
  maxi_in=inpat[0];
  mini_out=outpat[0];
  maxi_out=outpat[0];
  for (i=1;i<samples;i++)
  {
    if (inpat[i]>maxi_in) maxi_in=inpat[i];
    else if (inpat[i]<mini_in) mini_in=inpat[i];
    if (outpat[i]>maxi_out) maxi_out=outpat[i];
    else if (outpat[i]<mini_out) mini_out=outpat[i];
  }
  for (i=0;i<samples;i++)
  {
    inpat[i]=((upper_in-lower_in)*inpat[i]+maxi_in*lower_in-mini_in*upper_in)/(maxi_in-mini_in);
    outpat[i]=((upper_out-lower_out)*outpat[i]+maxi_out*lower_out-mini_out*upper_out)/(maxi_out-mini_out);
  }
}
/*-----------------------------------------------------------------------*/
/*  return to the actual space, samples belonging to [lower,upper]       */

void denormal_val(int end)
{
  int i,j;
  float in[samples],out[samples],outp[samples];
  float sum;
  FILE *datumtel;
  for (i=0;i<samples;i++)
  {
    in[i]=inpat[i];
    out[i]=((maxi_out-mini_out)*outpat[i]+(upper_out*mini_out-lower_out*maxi_out))/(upper_out-lower_out);
    outp[i]=((maxi_out-mini_out)*output[i]+(upper_out*mini_out-lower_out*maxi_out))/(upper_out-lower_out);
  }
  rmse=0.0;
  for (i=0;i<samples;i++) rmse=rmse+(1/(samples+0.0))*sqr(outp[i]-out[i]);
  rmse=sqrt((double)rmse);
  printf( "RMSE=%f\n",rmse );
  if (end==true)
  {
    datumtel=fopen("datumtel.dat","w");
    for (i=0;i<samples;i++)
    {
      in[i]=((maxi_in-mini_in)*in[i]+(upper_in*mini_in-lower_in*maxi_in))/(upper_in-lower_in);
      fprintf( datumtel,"%f     %f     %f\n",in[i],out[i],outp[i] );
    }
    fclose( datumtel );
  }
}
/*----------------------------------------------------------------------*/
void memb_value(int iter)
{
  int i,j,k;
  FILE *f1;
  if (iter==0) f1=fopen("init_memb.dat","w");
  for (k=0;k<samples;k++)
    for (j=0;j<rules;j++)
    {
      if (-0.5*sqr(wprm[j]-inpat[k])/sqr(wprdev[j])<-10.0) membshp[k][j]=0;
      else membshp[k][j]=exp(-0.5*sqr(wprm[j]-inpat[k])/sqr(wprdev[j]));
      if (iter==0)
      {
        if (j==(rules-1)) fprintf(f1,"%f\n",membshp[k][j]);
		else fprintf(f1,"%f ",membshp[k][j]);
      }
    }
    if (iter==0) fclose(f1);
}
/*-----------------------------------------------------------------------*/
void deg_ful(void)
{
  int i,j,k;
  float x,x1;
  for (k=0;k<samples;k++)
    for (j=0;j<rules;j++) degful[k][j]=membshp[k][j];
}
/*-----------------------------------------------------------------------*/
void find_con_w(int M)
{
  int i,j,k;
  float sum,sum_rules,sum1;
  float wprev[2][rules];
  float a[2*rules],Sprev[2*rules][2*rules];
  float Scur[2*rules][2*rules],Siai1[2*rules];
  float Siai1sqr[2*rules][2*rules];
  float Si1ai1[2*rules],ai1wprev,gamma;
  gamma=1000.0;

  /*  initial Sprev,wprev  */
  for (j=0;j<2;j++)
    for (k=0;k<rules;k++) wprev[j][k]=0.0;
  for (j=0;j<M;j++)
    for (k=0;k<M;k++)
    {
      if (k==j) Sprev[k][j]=gamma;
      else Sprev[k][j]=0.0;
    }

  for (i=0;i<samples;i++)
  {
    sum_rules=0.0;
    for (j=0;j<rules;j++) sum_rules=sum_rules+degful[i][j];
    for (j=0;j<rules;j++)
    {
      a[j*2]=degful[i][j]/sum_rules;
      a[j*2+1]=degful[i][j]*inpat[i]/sum_rules;
    }


    /*  matrix Si*a[i+1]  */
    for (j=0;j<M;j++)
    {
      Siai1[j]=0.0;
      for (k=0;k<M;k++) Siai1[j]=Siai1[j]+Sprev[j][k]*a[k];
    }

    /*  matrix Si*a[i+1]*transposed(Si*a[i+1])  */
    /* IMPORTANT COMMENT:
       Matrix S is symmetrical, thus:
       transposed(a[i+1])*Si=transposed(a[i+1])*transposed(Si)
			    =transposed(Si*a[i+1])=transposed(Siai1) */
    for (j=0;j<M;j++)
      for (k=0;k<M;k++) Siai1sqr[j][k]=Siai1[j]*Siai1[k];

    /*  matrix transposed(a[i+1])*Si*a[i+1]  */
    sum1=0.0;
    for (j=0;j<M;j++) sum1=sum1+a[j]*Siai1[j];

    /*  matrix S[i+1]  */
    for (j=0;j<M;j++)
      for (k=0;k<M;k++) Scur[j][k]=Sprev[j][k]-Siai1sqr[j][k]/(1+sum1);

    /*  matrix S[i+1]*a[i+1]  */
    for (j=0;j<M;j++)
    {
      Si1ai1[j]=0;
      for (k=0;k<M;k++) Si1ai1[j]=Si1ai1[j]+Scur[j][k]*a[k];
    }

    /*  transposed(a[i+1])*wprev  */
    ai1wprev=0;
    for (j=0;j<rules;j++)
      for (k=0;k<2;k++) ai1wprev=ai1wprev+a[j*2+k]*wprev[k][j];

    /*  w  */
    for (j=0;j<rules;j++)
      for (k=0;k<2;k++)
	{
	  w[k][j]=wprev[k][j]+Si1ai1[j*2+k]*(outpat[i]-ai1wprev);
	}

    /*  update  */
    for (j=0;j<M;j++)
      for (k=0;k<M;k++) Sprev[j][k]=Scur[j][k];
    for (j=0;j<rules;j++)
      for (k=0;k<2;k++) wprev[k][j]=w[k][j];
  }
}
/*-----------------------------------------------------------------------*/
void fuz_out(void)
{
  int i,j,k,l;
  float x;
  for (k=0;k<samples;k++)
    for (j=0;j<rules;j++) fuzout[k][j]=w[0][j]+w[1][j]*inpat[k];
}
/*-----------------------------------------------------------------------*/
void out_put(void)
{
  int i,j,k;
  float x,y;
  for (k=0;k<samples;k++)
  {
    x=0.0;
    y=0.0;
    for (j=0;j<rules;j++)
    {
      x=x+degful[k][j];
      y=y+degful[k][j]*fuzout[k][j];
    }
    output[k]=y/x;
  }
}
/*-----------------------------------------------------------------------*/
void lagr_mult(void)
{
  int j,i,k;
  float x,y;
  for (i=0;i<samples;i++)
  {
    lagr_out[i]=-2*(outpat[i]-output[i])/samples;
    x=0.0;
    for (j=0;j<rules;j++) x=x+degful[i][j];
    for (j=0;j<rules;j++)
    {
      y=(outpat[i]-output[i])*(fuzout[i][j]-output[i]);
      lagrdof[i][j]=-2*y/(x*samples);
      lagrfout[i][j]=-2*(outpat[i]-output[i])*degful[i][j]/(x*samples);
    }
  }
}
/*-----------------------------------------------------------------------*/
void aux_le_lf(void)  /* These are practically the gradients */
{
  int j,i,l;
  float x,y,xm,xdev;
  for (j=0;j<rules;j++)
  {
    xm=0; xdev=0;
    for (l=0;l<steps_back;l++)
    {
      xm=xm+wprevm[l][j];
      xdev=xdev+wprevdev[l][j];
    }
    lfm[j]=2*(wprm[j]-xm/(steps_back+0.0));
    lfdev[j]=2*(wprdev[j]-xdev/(steps_back+0.0));
  }
  for (j=0;j<rules;j++)
  {
    x=0;
    y=0;
    for (i=0;i<samples;i++)
    {
      x=x+lagrdof[i][j]*degful[i][j]*(inpat[i]-wprm[j])/sqr(wprdev[j]);
      y=y+lagrdof[i][j]*degful[i][j]*sqr(inpat[i]-wprm[j])/(sqr(wprdev[j])*wprdev[j]);
    }
    lem[j]=x;
    ledev[j]=y;
  }
}
/*-----------------------------------------------------------------------*/
void aux_i_ef(void)
{
  int k,i,j;
  float x1=0; float x3=0; float x5=0;
  for (i=0;i<rules;i++)
  {
    x1=x1+sqr(lfm[i])+sqr(lfdev[i]);
    x3=x3+sqr(lem[i])+sqr(ledev[i]);
    x5=x5+lfm[i]*lem[i]+lfdev[i]*ledev[i];
  }
  iffp=x1;
  ieep=x3;
  iefp=x5;
}
/*-----------------------------------------------------------------------*/
void perturb(float dep)
{
  int j,i,k;
  float xp;
  xp=sqrt((double)((ieep*sqr(dpp)-sqr(dep))/(iffp*ieep-sqr(iefp))));
  for (j=0;j<rules;j++)
  {
    dwprm[j]=xp*(lfm[j]-lem[j]*iefp/ieep)+lem[j]*dep/ieep;
    dwprdev[j]=xp*(lfdev[j]-ledev[j]*iefp/ieep)+ledev[j]*dep/ieep;
  }
}
/*-----------------------------------------------------------------------*/
void cor_weights(void)
{
  int i,k,j;
  for (j=0;j<rules;j++)
  {
    wprm[j]=wprm[j]+dwprm[j];
    if (wprm[j]<-1.0) wprm[j]=-1.0;
    if (wprm[j]>1.0) wprm[j]=1.0;
    wprdev[j]=wprdev[j]+dwprdev[j];
    if (wprdev[j]<0.01) wprdev[j]=0.01;
  }
}
/*-----------------------------------------------------------------------*/
void prev_epoch(void)
{
  int j,i,k,l;
  for (j=0;j<rules;j++)
  {
    for (l=(steps_back-1);l>0;l--)
    {
      wprevm[l][j]=wprevm[l-1][j];
      wprevdev[l][j]=wprevdev[l-1][j];
    }
    wprevm[0][j]=wprm[j];
    wprevdev[0][j]=wprdev[j];
  }
}
/*-----------------------------------------------------------------------*/
void initialization(void)
{
  int j,i,l;
  for (j=0;j<rules;j++) /* uniform partition of the iput space */
  {
    wprm[j]=-1+2*j/(rules-1.0);
    wprevm[0][j]=wprm[j]+0.001;
    for (l=1;l<steps_back;l++)  wprevm[l][j]=0;
  }

  for (j=0;j<rules;j++) /* uniform partition of the iput space */
  {
//    wprdev[j]=0.5/(rules-1.0);
    wprdev[j]=0.5*(wprm[1]-wprm[0])/(sqrt(-2*log(g)));
    printf("\nmean[%d]=%f    sigma[%d]=%f\n",j,wprm[j],j,wprdev[j]);
    wprevdev[0][j]=wprdev[j]+0.001;
    for (l=1;l<steps_back;l++)  wprevdev[l][j]=0;
  }
}
/*-----------------------------------------------------------------------*/
void features(void)
{
  int i,j;
  fprintf( dfiniter,"System and Learning Characteristics:\n" );
  fprintf( dfiniter,"-----------------------------------\n" );
  fprintf( dfiniter,"number of rules:%d\n",rules );
  fprintf( dfiniter,"samples:%d     maximum number of iterations:%d\n",samples,iterat );
  fprintf( dfiniter,"low1:%f     up2:%f     low_con2:%f     low_prem2:%f\n",low1,up2,low_con2,low_prem2 );
  fprintf( dfiniter,"dpparh:%f",dpparh );
  fprintf( dfiniter,"ksi:%f     g:%f\n",ksi,g );
  for (i=0;i<2;i++) fprintf( dfiniter,"\n" );
  fprintf( dfiniter,"INITIAL MEANS AND DEVIATIONS\n" );
  for (j=0;j<rules;j++)
    fprintf( dfiniter,"wprm[%d]=%f     wprdev[%d]=%f\n",j,wprm[j],j,wprdev[j] );
}

