/* ageadjs.c file contains four C language subroutines:                */
/* ageadj1.c for age adjustment of milk, fat, and protein              */
/* ageadj4.c to call the above routine from Fortran                    */
/* adjscs.c  for age adjustment of somatic cell score                  */
/* adjscs4.c to call the above routine from Fortran                    */
/*                                                                     */
/*---------------------------------------------------------------------*/
/* ageadj1.c   Apply age, month of calving and prev days open adjustments */
/*-----------------------------------------------------------------------*/
#define DOMEAN 120

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

static int agelim[6][2] = {16,36, 28,56, 40,74, 52,86, 64,98, 76,165}
    ,regdot[] = {0,0,1,4,4,1,3,2,3,5,5,6};

#include "ageadj.h"

void ageadj1(char * breed, char * parity, char * fmoptr, char * state
   ,short * ageptr, short * dim, short * doprev
   ,short * docurrent, unsigned short * yield,unsigned short * adjyld) {

  char breeds[] = "ABGHJM", *p, regbrd[] = {0,0,0,1,2,3,2}
    ,parity1;
  float facm, facf, facp, do150p, age2, doprev2
    ,pdoeffm,pdoefff,pdoeffp,agem,agef,agep;
  int brd, brddo, reg, lac, lacdo, pdo, regdo, age, fmo;

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
  parity1 = *parity;
  if (parity1 > 6)
    parity1 = 6;
  lac = parity1 - 1;

  age = *ageptr;
 /* If age outside of range for parity, use age limit */
  if (age < agelim[lac][0])
    age =  agelim[lac][0];
  else if (age > agelim[lac][1])
    age = agelim[lac][1];

  pdo = *doprev;
  if ((parity1 > 1 && pdo < 20) || pdo > 305)
    pdo = DOMEAN;

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
      if (yield[2] > 0)
        pdoeffp = hpdop[regdo][lacdo][0] + hpdop[regdo][lacdo][1] * pdo
         + hpdop[regdo][lacdo][2] * doprev2;
      else
        pdoeffp = 0.;
      if (pdo > 150) {
        do150p = (pdo - 150)*(pdo - 150);
        pdoeffm += hmdop[regdo][lacdo][3] * do150p;
        pdoefff += hfdop[regdo][lacdo][3] * do150p;
        if (yield[2] > 0)
          pdoeffp += hpdop[regdo][lacdo][3] * do150p;
      }
    } else {             /*  Other Breeds */
      brddo = (brd < 3) ? brd : brd - 1;
      pdoeffm = mdop[brddo][lacdo][0]+ mdop[brddo][lacdo][1] * pdo
        + mdop[brddo][lacdo][2] * doprev2;
      pdoefff = fdop[brddo][lacdo][0] + fdop[brddo][lacdo][1] * pdo
        + fdop[brddo][lacdo][2] * doprev2;
      if (yield[2] > 0)
        pdoeffp = pdop[brddo][lacdo][0] + pdop[brddo][lacdo][1] * pdo
          + pdop[brddo][lacdo][2] * doprev2;
      else
        pdoeffp = 0.;
      if (pdo > 150) {
        do150p = (pdo - 150)*(pdo - 150);
        pdoeffm += mdop[brddo][lacdo][3] * do150p;
        pdoefff += fdop[brddo][lacdo][3] * do150p;
        if (yield[2] > 0)
          pdoeffp += pdop[brddo][lacdo][3] * do150p;
      }
    }
  } else {
    pdoeffm = 0.;
    pdoefff = 0.;
    pdoeffp = 0.;
  }  /* End of Previous Days Open */


/* Apply Prev Days Open, Month, and Age Values to create Factor */
  fmo = *fmoptr - 1;
  age2 = age*age;

  switch (brd) {
  case 3: /* Holstein */
    agem = hmequ[reg][lac][0] * age
      + hmequ[reg][lac][1] * age2 + hmregtm[reg];
    facm = hmregtm[reg]/(pdoeffm + hmequ0[reg][lac][fmo] + agem);
    agef = hfequ[reg][lac][0] * age
      + hfequ[reg][lac][1] * age2 + hfregtm[reg];
    facf = hfregtm[reg]/(pdoefff + hfequ0[reg][lac][fmo] + agef);
    if (yield[2] > 0) {
      agep = hpequ[reg][lac][0] * age
        + hpequ[reg][lac][1] * age2 + hpregtm[reg];
      facp = hpregtm[reg]/(pdoeffp + hpequ0[reg][lac][fmo] + agep);
    }
    break;

  case 4: /* Jersey   */
    agem = jmequ[reg][lac][0] * age
      + jmequ[reg][lac][1] * age2 + jmregtm[reg];
    facm = jmregtm[reg]/(pdoeffm + jmequ0[reg][lac][fmo] + agem);
    agef = jfequ[reg][lac][0] * age
      + jfequ[reg][lac][1] * age2 + jfregtm[reg];
    facf = jfregtm[reg]/(pdoefff + jfequ0[reg][lac][fmo] + agef);
    if (yield[2] > 0) {
      agep = jpequ[reg][lac][0] * age
         + jpequ[reg][lac][1] * age2 + jpregtm[reg];
      facp = jpregtm[reg]/(pdoeffp + jpequ0[reg][lac][fmo] + agep);
    }
    break;

  case 2: /* Guernsey */
    agem = gmequ[reg][lac][0] * age
      + gmequ[reg][lac][1] * age2 + gmregtm[reg];
    facm = gmregtm[reg]/(pdoeffm + gmequ0[reg][lac][fmo] + agem);
    agef = gfequ[reg][lac][0] * age
      + gfequ[reg][lac][1] * age2 + gfregtm[reg];
    facf = gfregtm[reg]/(pdoefff + gfequ0[reg][lac][fmo] + agef);
    if (yield[2] > 0) {
      agep = gpequ[reg][lac][0] * age
        + gpequ[reg][lac][1] * age2 + gpregtm[reg];
      facp = gpregtm[reg]/(pdoeffp + gpequ0[reg][lac][fmo] + agep);
    }
    break;

  case 1: /* Brown Swiss */
    agem = bmequ[reg][lac][0] * age
      + bmequ[reg][lac][1] * age2 + bmregtm[reg];
    facm = bmregtm[reg]/(pdoeffm + bmequ0[reg][lac][fmo] + agem);
    agef = bfequ[reg][lac][0] * age
        + bfequ[reg][lac][1] * age2 + bfregtm[reg];
    facf = bfregtm[reg]/(pdoefff + bfequ0[reg][lac][fmo] + agef);
    if (yield[2] > 0) {
      agep = bpequ[reg][lac][0] * age
        + bpequ[reg][lac][1] * age2 + bpregtm[reg];
      facp = bpregtm[reg]/(pdoeffp + bpequ0[reg][lac][fmo] + agep);
    }
    break;

  case 0: /* Ayrsyire */
    agem = amequ[reg][lac][0] * age
      + amequ[reg][lac][1] * age2 + amregtm[reg];
    facm = amregtm[reg]/(pdoeffm + amequ0[reg][lac][fmo] + agem);
    agef = afequ[reg][lac][0] * age
      + afequ[reg][lac][1] * age2 + afregtm[reg];
    facf = afregtm[reg]/(pdoefff + afequ0[reg][lac][fmo] + agef);
    if (yield[2] > 0) {
      agep = apequ[reg][lac][0] * age
        + apequ[reg][lac][1] * age2 + apregtm[reg];
      facp = apregtm[reg]/(pdoeffp + apequ0[reg][lac][fmo] + agep);
    }
    break;

  case 5: /* Milking Shorthorn */
    agem = mmequ[lac][0] * age
      + mmequ[lac][1] * age2 + mmregtm;
    facm = mmregtm/(pdoeffm + mmequ0[lac][fmo] + agem);
    agef = mfequ[lac][0] * age
      + mfequ[lac][1] * age2 + mfregtm;
    facf = mfregtm/(pdoefff + mfequ0[lac][fmo] + agef);
    if (yield[2] > 0) {
      agep = mpequ[lac][0] * age
         + mpequ[lac][1] * age2 + mpregtm;
      facp = mpregtm/(pdoeffp + mpequ0[lac][fmo] + agep);
    }
  }
  adjyld[0] = facm*yield[0] + .5;
  if (yield[1] > 0)
    adjyld[1] = facf*yield[1] + .5;
  else
    adjyld[1] = 0;
  if (yield[2] > 0)
    adjyld[2] = facp*yield[2] + .5;
  else
    adjyld[2] = 0;
}
/****************************************************************/
/* Function: ageadj4.c                                           */
/* Purpose:  Converts 2-byte integers to one-byte characters    */
/* Preconditions:  Must be called by the fmt4dcr.f program,     */
/*     which passes 10 variables to ageadj4.c.                   */
/* Postconditions: All ten variables are passed on to the       */ 
/*     ageadj1.c function for adjustment. Only adjyld1 is passed */
/*     back to the fmt4dcr.f function                           */
/*                                        Written by Lori Smith */
/****************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <true.h>

void ageadj1(         char*  breed1,          char*   newpar,
                      char*  newmon,          char*   newstst,
                      short* age1,            short*  fill,
                      short* prev,            short*  cur,
             unsigned short* yld,    unsigned short*  adjyld1);

int ageadj4(          char*   breed1,        short* par,
                      short*  mon,           short* st,
                      short*  age1,          short* fill,
                      short*  prev,          short* cur,
             unsigned short*  yld,  unsigned short* adjyld1)
   {
   char newpar, newmon, newst ;
   newpar = *par;
   newmon = *mon;
   newst  = *st;
/*  
   printf("Varpak pointer test = %d %d %d \n", *par, *mon, *st);
   printf("Varpak int pointer -> char test = %d %d %d \n", newpar, newmon, newst);
*/  
   ageadj1(breed1, &newpar, &newmon, &newst,   age1,
           fill,   prev,   cur,    yld,     adjyld1);
   return (0);
   }
/* adjscs.c      Apply adjustments for DIM, age, month of calving */
/*               Inputs, outputs, and adjustments are 100 * SCS   */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "adjscs.h"

int adjscs(char breed, char parity, char fresh_mo, char state ,short age
   ,short dim, short scs) {

  static int first = TRUE, hidim = 0, lodim = 0;
  int i, j, lact ,facbrd, reg, adjscs;
  double rint(double x);
  static double dimfac[306][2][2],agefac[121][2],monfac[13][4][2]; 
  FILE *fopen(), *fadj;

  if (first) {
    if ((fadj = fopen ("adjust.scs","r")) == NULL) {
      fprintf(stderr,"adjscs.c failed to open adjust.scs\n");
      exit(62);
    }

    /*
    printf("ageadjs.c breed    = %c\n", breed);
    printf("ageadjs.c parity   = %c\n", parity);
    printf("ageadjs.c fresh_mo = %c\n", fresh_mo);
    printf("ageadjs.c state    = %c\n", state);
    printf("ageadjs.c age      = %d\n", age);
    printf("ageadjs.c dim      = %d\n", dim);
    printf("ageadjs.c scs      = %d\n", scs);
    */

    /* LOAD ADJUSTMENT FACTORS FOR DIM, AGE, AND MONTH.
      DIMENSIONS:
                 dimfac(DIM,lactation,breed)
                 agefact(age,breed)
                 monfac(month,region,breed)
    */
       /* Days in Milk Adjustments */
    for (j=0;j<2;j++) {
      for (i=15;i<=305;i++) {
        fscanf(fadj,"%lf %lf",&dimfac[i][j][0],&dimfac[i][j][1]);
          /*  printf("%d %d %8.1f %8.1f\n"
           ,j,i,dimfac[i][j][0], dimfac[i][j][1]);  */
      }
    }
       /* Age at calving Adjustments */
    for (j=18;j<=120;j++) {
      fscanf(fadj,"%lf %lf",&agefac[j][0],&agefac[j][1]);
      /* printf("%d %8.1f  %8.1f\n",j,agefac[j][0], agefac[j][1]); */
    }

       /* Month of calving Adjustments */
    for (j=0;j<4;j++) {
      for (i=1;i<=12;i++) {
        fscanf(fadj,"%lf %lf",&monfac[i][j][0],&monfac[i][j][1]);
      }
    }
    first = FALSE;
  }
  if (parity == 1) {
    lact = 0;
  } else {
    lact = 1;
  }
  if (breed == 'J' || breed == 'G') facbrd = 1;
  else facbrd = 0;

        /*  calculate age at calving and get freshmo */
  if (age < 18) age = 18;
  else if (age > 120) age = 120;
 
  if (state <= 23) reg = 0;
  else if (state <= 48) reg = 1;
  else if (state <= 74) reg = 2;
      else reg = 3;
  if (dim > 305) {
    hidim++;
    if (hidim <= 5) printf("adjscs.c DIM=%d\n",dim);
    dim = 305;
  }  else if (dim < 15) {
    lodim++;
    if (lodim <= 5) printf("adjscs.c DIM=%d\n",dim);
    dim = 15;
  }
/*                                 Additive SCS adjustment factor */
  adjscs = rint(scs - dimfac[dim][lact][facbrd]
                    - monfac[fresh_mo][reg][facbrd]
                    + (agefac[age][facbrd] - 1) * 300);

/*                     Previous multiplicative factor used 1994-2009 
  adjscs =  rint((scs - dimfac[dim][lact][facbrd]
         - monfac[fresh_mo][reg][facbrd]) * agefac[age][facbrd]); */
  if (adjscs < 1) {
    adjscs = 1;
  }
  if (adjscs > 999) {
    adjscs = 999;
  }
  return adjscs;
}
/****************************************************************/
/* Function: adjscs4.c                                           */
/* Purpose:  Converts 2-byte integers to one-byte characters    */
/* Preconditions:  Must be called by the fmt4dcr.f program,     */
/*     which passes 7 values to adjscs4.c.                       */
/* Postconditions: All seven variables are passed on to the     */ 
/*     adjscs.c function for adjustment. Returns an integer to  */
/*     the fmt4dcr.f program.                                   */
/*                                        Written by Lori Smith */
/****************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int adjscs  (char  brd,     char  newpar,  char   newmon,
            char  newst,   short age1,    short  dim,
            short somcell);

int adjscs4 (char  *breed1, short *parity, short  *month,
            short *state,  short *age2,   short  *filler,
            short *scs )         
   {
   char brd     = *breed1;
   char newpar  = *parity;
   char newmon  = *month;
   char newst   = *state;
   short age1    = *age2;
   short days    = *filler;
   short somcell = *scs;
   int adjust;
/*   
   printf("adjscs4 pointer test = %d %d %d \n", *parity, *month, *state);
   printf("adjscs4 int pointer -> char value test = %d %d %d \n", newpar, newmon, newst);
*/  
   adjust = adjscs(brd, newpar, newmon, newst, age1, days, somcell);
   return adjust;
   }

