/* ageadj.c   Apply age-month of calving and days open adjustments */
/*            agebase = 0 adjusts to mature, 36 to 36 months, etc  */
/*                        based on factors of Schutz, 1994         */
/*            Oct 1994 George Wiggans - Original program           */
/*            Sep 2004 Paul VanRaden  - 36 month option added      */
/*-----------------------------------------------------------------*/
#define DOMEAN 120

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void ageadj(char * breed, short * agefrsh, short * fryear, short * frmonth
   ,short * parity, short * state, short * doprev, short * agebase
   ,double * agefacm, double * agefacf, double * agefacp) {

#include "age.h"
/*
#include "dosol1.c"
#include "age.equ"
#include "region.def"
#include "reg.mean"
*/
  char breeds[] = "ABGHJMW", *p, DEE='D'
   ,regbrd[] = { 0,0,0,1,2,3,2};
  double  pdoeffm, pdoefff, pdoeffp
     ,facm, facf, facp, do150, do150p, age2, doprev2
     ,mregtm, fregtm, pregtm
     ,avgfacm, avgfacf, avgfacp;
  int brd, brddo, reg, mo, yrgp, yrgpp
     ,age, lac, lacdo, pdo, regdo, i, j, pass
     ,month, months;

static int agelim[6][2] = {16,36, 28,56, 40,74, 52,86, 64,98, 76,165}
    ,regdot[] = {0,0,1,4,4,1,3,2,3,5,5,6};

  if ((p = strchr(breeds,*breed)) == NULL)
    brd = 3;    /* default to Holstein */
  else
    brd =  (int) (p - breeds);
  {
    int brdreg;
    brdreg =  regbrd[brd];
    reg = region[*state][brdreg] - 1;  /* Storage starts at 0 */
  }
  if (parity == 0) {
    fprintf(stderr,"Parity = 0 in function ageadj\n");
    exit(25);
  }
  lac = *parity - 1;
  if (lac > 5) lac = 5;
 /* If age outside of range for parity, use age limit */
  age = *agefrsh;
  if (age < agelim[lac][0])
    age =  agelim[lac][0];
  else if (age > agelim[lac][1])
    age = agelim[lac][1];
  mo = *frmonth -1;
  pdo = *doprev;
  if ((*parity > 1 && pdo < 20) || pdo > 305)
    pdo = DOMEAN;

/*   ME factor on first pass, 36-month factor on 2nd pass */
for (pass = 1; pass <= 2; pass++) {
  months = 1;
  if (pass == 2) {
    age = *agebase;
    if (age == 0) age = 36;
    lac = (*agebase - 18) / 13;
    pdo = 90;
    if(*agebase > 0) months = 12;
    avgfacm = 0;
    avgfacf = 0;
    avgfacp = 0;
    }
for (month = 1; month <= months; month++) {
  if (pass == 2) mo = month -1;
/* Calculate previous Days Open Values */
  if (lac > 0) {
    lacdo = (lac > 2) ? 2: lac; /* Values 0, 1, 2 for lactations 1, 2, 3 */
    lacdo--;
    regdo = regdot[reg];
    doprev2 = pdo * pdo;
                         /* Holstein */
    if (brd == 3) {  /* Holstein */
      pdoeffm = hmdop[regdo][lacdo][0] + hmdop[regdo][lacdo][1] * pdo
       + hmdop[regdo][lacdo][2] * doprev2;
      pdoefff = hfdop[regdo][lacdo][0] + hfdop[regdo][lacdo][1] * pdo
       + hfdop[regdo][lacdo][2] * doprev2;
      pdoeffp = hpdop[regdo][lacdo][0] + hpdop[regdo][lacdo][1] * pdo
       + hpdop[regdo][lacdo][2] * doprev2;
      if (pdo > 150) {
        do150p = (pdo - 150)*(pdo - 150);
        pdoeffm += hmdop[regdo][lacdo][3] * do150p;
        pdoefff += hfdop[regdo][lacdo][3] * do150p;
        pdoeffp += hpdop[regdo][lacdo][3] * do150p;
      }
    } else {             /*  Other Breeds */
      brddo = (brd < 3) ? brd : brd - 1;
      pdoeffm = mdop[brddo][lacdo][0]+ mdop[brddo][lacdo][1] * pdo
        + mdop[brddo][lacdo][2] * doprev2;
      pdoefff = fdop[brddo][lacdo][0] + fdop[brddo][lacdo][1] * pdo
        + fdop[brddo][lacdo][2] * doprev2;
      pdoeffp = pdop[brddo][lacdo][0] + pdop[brddo][lacdo][1] * pdo
        + pdop[brddo][lacdo][2] * doprev2;
      if (pdo > 150) {
        do150p = (pdo - 150)*(pdo - 150);
        pdoeffm += mdop[brddo][lacdo][3] * do150p;
        pdoefff += fdop[brddo][lacdo][3] * do150p;
        pdoeffp += pdop[brddo][lacdo][3] * do150p;
      }
    }
  } else {
    pdoeffm = 0.;
    pdoefff = 0.;
    pdoeffp = 0.;
  }  /* End of Previous Days Open */

  if (*fryear >= 1987) {
    yrgp = 4;
    yrgpp = 1;
  } else {
    yrgp = (*fryear >= 1969) ?  (*fryear - 1963)/6 : 0;
    yrgpp = 0;
  }
  lacdo = (lac > 2) ? 2: lac;
  age2 = age*age;
  switch (brd) {
  case 3: /* Holstein */
/*   regtm approximates region means for earlier year groups */
    mregtm = hmregt[reg] * mratio[yrgp];
    fregtm = hfregt[reg] * fratio[yrgp];
    pregtm = hpregt[reg] * pratio[yrgp];
    facm = mregtm/(mregtm + pdoeffm 
      + hmequ0[yrgp][reg][lac][mo]
      + hmequ[yrgp][reg][lac][0] * age
      + hmequ[yrgp][reg][lac][1] * age2);
    facf = fregtm/(fregtm + pdoefff
      + hfequ0[yrgp][reg][lac][mo]
      + hfequ[yrgp][reg][lac][0] * age
      + hfequ[yrgp][reg][lac][1] * age2);
    facp = pregtm/(pregtm + pdoeffp
      + hpequ0[yrgpp][reg][lac][mo]
      + hpequ[yrgpp][reg][lac][0] * age
      + hpequ[yrgpp][reg][lac][1] * age2);
    break;

  case 4: /* Jersey   */
    mregtm = jmregt[reg] * mratio[yrgp];
    fregtm = jfregt[reg] * fratio[yrgp];
    pregtm = jpregt[reg] * pratio[yrgp];
    facm = mregtm/(mregtm + pdoeffm 
      + jmequ0[yrgp][reg][lac][mo]
      + jmequ[yrgp][reg][lac][0] * age
      + jmequ[yrgp][reg][lac][1] * age2);
    facf = fregtm/(fregtm + pdoefff
      + jfequ0[yrgp][reg][lac][mo]
      + jfequ[yrgp][reg][lac][0] * age
      + jfequ[yrgp][reg][lac][1] * age2);
    facp = pregtm/(pregtm + pdoeffp
      + jpequ0[yrgpp][reg][lac][mo]
      + jpequ[yrgpp][reg][lac][0] * age
      + jpequ[yrgpp][reg][lac][1] * age2);
    break;

  case 2: /* Guernsey */
    mregtm = gmregt[reg] * mratio[yrgp];
    fregtm = gfregt[reg] * fratio[yrgp];
    pregtm = gpregt[reg] * pratio[yrgp];
    facm = mregtm/(mregtm + pdoeffm 
      + gmequ0[yrgp][reg][lac][mo]
      + gmequ[yrgp][reg][lac][0] * age
      + gmequ[yrgp][reg][lac][1] * age2);
    facf = fregtm/(fregtm + pdoefff
      + gfequ0[yrgp][reg][lac][mo]
      + gfequ[yrgp][reg][lac][0] * age
      + gfequ[yrgp][reg][lac][1] * age2);
    facp = pregtm/(pregtm + pdoeffp
      + gpequ0[yrgpp][reg][lac][mo]
      + gpequ[yrgpp][reg][lac][0] * age
      + gpequ[yrgpp][reg][lac][1] * age2);
    break;

  case 1: /* Brown Swiss */
    mregtm = bmregt[reg] * mratio[yrgp];
    fregtm = bfregt[reg] * fratio[yrgp];
    pregtm = bpregt[reg] * pratio[yrgp];
    facm = mregtm/(mregtm + pdoeffm 
      + bmequ0[yrgp][reg][lac][mo]
      + bmequ[yrgp][reg][lac][0] * age
      + bmequ[yrgp][reg][lac][1] * age2);
    facf = fregtm/(fregtm + pdoefff
      + bfequ0[yrgp][reg][lac][mo]
      + bfequ[yrgp][reg][lac][0] * age
      + bfequ[yrgp][reg][lac][1] * age2);
    facp = pregtm/(pregtm + pdoeffp
      + bpequ0[yrgpp][reg][lac][mo]
      + bpequ[yrgpp][reg][lac][0] * age
      + bpequ[yrgpp][reg][lac][1] * age2);
    break;

  case 0: /* Ayrshire */
    mregtm = amregt[reg] * mratio[yrgp];
    fregtm = afregt[reg] * fratio[yrgp];
    pregtm = apregt[reg] * pratio[yrgp];
    facm = mregtm/(mregtm + pdoeffm 
      + amequ0[yrgp][reg][lac][mo]
      + amequ[yrgp][reg][lac][0] * age
      + amequ[yrgp][reg][lac][1] * age2);
    facf = fregtm/(fregtm + pdoefff
      + afequ0[yrgp][reg][lac][mo]
      + afequ[yrgp][reg][lac][0] * age
      + afequ[yrgp][reg][lac][1] * age2);
    facp = pregtm/(pregtm + pdoeffp
      + apequ0[yrgpp][reg][lac][mo]
      + apequ[yrgpp][reg][lac][0] * age
      + apequ[yrgpp][reg][lac][1] * age2);
    break;

  case 5: /* Milking Shorthorn */
    mregtm = mmregt * mratio[yrgp];
    fregtm = mfregt * fratio[yrgp];
    pregtm = mpregt * pratio[yrgp];
    facm = mregtm/(mregtm + pdoeffm 
      + mmequ0[yrgp][lac][mo]
      + mmequ[yrgp][lac][0] * age
      + mmequ[yrgp][lac][1] * age2);
    facf = fregtm/(fregtm + pdoefff
      + mfequ0[yrgp][lac][mo]
      + mfequ[yrgp][lac][0] * age
      + mfequ[yrgp][lac][1] * age2);
    facp = pregtm/(pregtm + pdoeffp
      + mpequ0[yrgpp][lac][mo]
      + mpequ[yrgpp][lac][0] * age
      + mpequ[yrgpp][lac][1] * age2);
  }
  if (pass == 1) {
    *agefacm = facm;
    *agefacf = facf;
    *agefacp = facp;
    }
/*       Obtain average across 12 months of year */
  if (pass == 2) {
    avgfacm = avgfacm + facm / months;
    avgfacf = avgfacf + facf / months;
    avgfacp = avgfacp + facp / months;
    }
/*       Divide by 36 month average factor */
  if (month == 12 && *agebase > 0) {
    *agefacm = *agefacm / avgfacm;
    *agefacf = *agefacf / avgfacf;
    *agefacp = *agefacp / avgfacp;
    }
  }
  }
}
/* printf("facm %f",facm); */
/*      Standardize lactation yield
  adjyld[0] = facm*yield[0] + .5;
  adjyld[1] = facf*yield[1] + .5;
  adjyld[2] = facp*yield[2] + .5;
*/
