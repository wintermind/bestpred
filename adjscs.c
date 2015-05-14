/* adjscs.c   Apply adjustments for DIM, age, month of calving */
/*            Inputs, outputs, and adjustments are 100 * SCS   */

/* Gary C. Fok  03/22/2001  Added <true.h>. */
/* Gary C. Fok  02/11/2002  Replaced nearest() with rint().  Added <math.h>. */
/* Paul VanRaden 7/30/2009  Additive adjustment for SCS instead of multiplicative */
#include <stdlib.h>
#include <stdio.h>
#include <true.h>
#include <math.h>
/* ADJUSTMENT FACTORS FOR DIM, AGE, AND MONTH.
  DIMENSIONS:
             dimfac(DIM,lactation,breed)
             agefact(age,breed)
             monfac(month,region,breed)
*/
#include <adjscs.h>

int32_t adjscs(char breed, char parity, char fresh_mo, char state ,short age
   ,short dim, short scs) {

  static int32_t first = TRUE, hidim = 0, lodim = 0;
  int32_t i, j, lact ,facbrd, reg, adjscs;
  double rint(double x);

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
 
  if (state == 66) reg = 3;
  else if (state <= 23) reg = 0;
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
/* JBC 08/17/2009 Rolled back because the adjustments are not actually applied
                  applied to the TD data in this function. Actually applying
                  then additively rather than multiplicatively will require
                  downstream changes in the routines that actually use the
                  age adjustments, e.g., bestpred.f90.
/*                                 Additive SCS adjustment factor */
/*
  adjscs = rint(scs - dimfac[dim][lact][facbrd]
                    - monfac[fresh_mo][reg][facbrd]
                    + (agefac[age][facbrd] - 1) * 300);
*/

/*                     Previous multiplicative factor used 1994-2009 */
  adjscs = rint((scs - dimfac[dim][lact][facbrd]
                     - monfac[fresh_mo][reg][facbrd])
                     * agefac[age][facbrd]);

/*
  printf("adjscs.c: scs    = %d\n", scs);
  printf("adjscs.c: dimfac = %d\n", dimfac[dim][lact][facbrd]);
  printf("adjscs.c: agefac = %d\n", monfac[fresh_mo][reg][facbrd]);
*/  

  if (adjscs < 1) {
    adjscs = 1;
  }
  if (adjscs > 999) {
    adjscs = 999;
  }
  return adjscs;
}
