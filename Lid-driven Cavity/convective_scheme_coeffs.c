#include<math.h>

void convective_scheme_coeffs(int i, int j, double De, double Dw, double Dn, double Ds, double Fe, double Fw, double Fn, double Fs, double **Sp, int convec_scheme, double **aP, double **aE, double **aW, double **aN, double **aS)
{

if (convec_scheme==1) { // upwind scheme    
    aE[i][j] = De + fmax(0.0, -Fe);         aN[i][j] = Dn + fmax(0.0, -Fn);
    aW[i][j] = Dw + fmax(Fw,  0.0);         aS[i][j] = Ds + fmax(Fs,  0.0);
 }
      
else if(convec_scheme==2) { //C-D scheme
    aE[i][j] = De - 0.5*Fe;              aN[i][j] = Dn - 0.5*Fn;
    aW[i][j] = Dw + 0.5*Fw;              aS[i][j] = Ds + 0.5*Fs; 
  } 
/*
 else if (convec_scheme==3) {  //PowerLaw scheme
    aE[i][j] = De * fmax(0, pow((1.0-0.1*fabs(Fe/De)),5.0)) + fmax(0.0, -Fe);
    aW[i][j] = Dw * fmax(0, pow((1.0-0.1*fabs(Fw/Dw)),5.0)) + fmax(0.0,  Fw);
    aN[i][j] = Dn * fmax(0, pow((1.0-0.1*fabs(Fn/Dn)),5.0)) + fmax(0.0, -Fn);
    aS[i][j] = Ds * fmax(0, pow((1.0-0.1*fabs(Fs/Ds)),5.0)) + fmax(0.0, -Fs);
  }
  */
  // calculate aP

aP[i][j] = aE[i][j] + aW[i][j] + aS[i][j] + aN[i][j] + (Fe-Fw)+(Fn-Fs) - Sp[i][j];
}

