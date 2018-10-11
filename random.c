// http://www.math.keio.ac.jp/~matumoto/ver980409.html

// This is the ``Mersenne Twister'' random number generator MT19937,
// which generates pseudorandom integers uniformly distributed in
// 0..(2^32 - 1) starting from any odd seed in 0..(2^32 - 1).  This
// version is a recode by Shawn Cokus (Cokus@math.washington.edu) on
// March 8, 1998 of a version by Takuji Nishimura (who had suggestions
// from Topher Cooper and Marc Rieffel in July-August 1997).
//
// Effectiveness of the recoding (on Goedel2.math.washington.edu, a DEC
// Alpha running OSF/1) using GCC -O3 as a compiler: before recoding:
// 51.6 sec. to generate 300 million random numbers; after recoding: 24.0
// sec. for the same (i.e., 46.5% of original time), so speed is now
// about 12.5 million random number generations per second on this
// machine.
//
// According to the URL <http://www.math.keio.ac.jp/~matumoto/emt.html>
// (and paraphrasing a bit in places), the Mersenne Twister is ``designed
// with consideration of the flaws of various existing generators,'' has
// a period of 2^19937 - 1, gives a sequence that is 623-dimensionally
// equidistributed, and ``has passed many stringent tests, including the
// die-hard test of G. Marsaglia and the load test of P. Hellekalek and
// S.  Wegenkittl.''  It is efficient in memory usage (typically using
// 2506 to 5012 bytes of static data, depending on data type sizes, and
// the code is quite short as well).  It generates random numbers in
// batches of 624 at a time, so the caching and pipelining of modern
// systems is exploited.  It is also divide- and mod-free.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as
// published by the Free Software Foundation (either version 2 of the
// License or, at your option, any later version).  This library is
// distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY, without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Library General Public
// License for more details.  You should have received a copy of the GNU
// Library General Public License along with this library; if not, write
// to the Free Software Foundation, Inc., 59 Temple Place, Suite 330,
// Boston, MA 02111-1307, USA.
//
// The code as Shawn received it included the following notice:
//
//   Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.  When you
//   use this, send an e-mail to <matumoto@math.keio.ac.jp> with an
//   appropriate reference to your work.
//
// It would be nice to CC: <Cokus@math.washington.edu> when you write.
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//
// uint32 must be an unsigned integer type capable of holding at least
// 32 bits; exactly 32 should be fastest, but 64 is better on an Alpha
// with GCC at -O3 optimization so try your options and see what's best
// for you
//

typedef unsigned long uint32;

#define N              (624)                // length of state vector
#define M              (397)                // a period parameter
#define K              (0x9908B0DFU)        // a magic constant
#define hiBit(u)       ((u) & 0x80000000U)  // mask all but highest bit
#define loBit(u)       ((u) & 0x00000001U)  // mask all but lowest bit
#define loBits(u)      ((u) & 0x7FFFFFFFU)  // mask the highest bit
#define mixBits(u, v)  (hiBit(u)|loBits(v)) // move hi bit u to hi bit v

static uint32   state[N+1]; // state vector + 1 extra to respect ANSI C
static uint32   *next;      // next random value is computed from here
static int      left = -1;  // can *next++ so many times before reloading


void
seedMT // VISIBLE 
(
   uint32 seed
)
//
// We initialize state[0..(N-1)] via the generator
//
//   x_new = (69069 * x_old) mod 2^32
//
// from Line 15 of Table 1, p. 106, Sec. 3.3.4 of Knuth's
// _The Art of Computer Programming_, Volume 2, 3rd ed.
//
// Notes (SJC): I do not know what the initial state requirements
// of the Mersenne Twister are, but it seems this seeding generator
// could be better.  It achieves the maximum period for its modulus
// (2^30) iff x_initial is odd (p. 20-21, Sec. 3.2.1.2, Knuth); if
// x_initial can be even, you have sequences like 0, 0, 0, ...;
// 2^31, 2^31, 2^31, ...; 2^30, 2^30, 2^30, ...; 2^29, 2^29 + 2^31,
// 2^29, 2^29 + 2^31, ..., etc. so I force seed to be odd below.
//
// Even if x_initial is odd, if x_initial is 1 mod 4 then
//
//   the          lowest bit of x is always 1,
//   the  next-to-lowest bit of x is always 0,
//   the 2nd-from-lowest bit of x alternates  ... 0 1 0 1 0 1 0 1 ... ,
//   the 3rd-from-lowest bit of x 4-cycles    ... 0 1 1 0 0 1 1 0 ... ,
//   the 4th-from-lowest bit of x has 8-cycle ... 0 0 0 1 1 1 1 0 ... ,
//    ...
//
// and if x_initial is 3 mod 4 then
//
//   the          lowest bit of x is always 1,
//   the  next-to-lowest bit of x is always 1,
//   the 2nd-from-lowest bit of x alternates  ... 0 1 0 1 0 1 0 1 ... ,
//   the 3rd-from-lowest bit of x 4-cycles    ... 0 0 1 1 0 0 1 1 ... ,
//   the 4th-from-lowest bit of x has 8-cycle ... 0 0 1 1 1 1 0 0 ... ,
//    ...
//
// The generator's potency (min. s>=0 with (69069-1)^s = 0 mod 2^32) is
// 16, which seems to be alright by p. 25, Sec. 3.2.1.3 of Knuth.  It
// also does well in the dimension 2..5 spectral tests, but it could be
// better in dimension 6 (Line 15, Table 1, p. 106, Sec. 3.3.4, Knuth).
//
// Note that the random number user does not see the values generated
// here directly since reloadMT() will always munge them first, so maybe
// none of all of this matters.  In fact, the seed values made here could
// even be extra-special desirable if the Mersenne Twister theory says
// so-- that's why the only change I made is to restrict to odd seeds.
//
{

    register uint32 x = (seed | 1U) & 0xFFFFFFFFU, *s = state;
    register int    j;

    for(left=0, *s++=x, j=N; --j;
          *s++ = (x*=69069U) & 0xFFFFFFFFU);
 }


uint32 reloadMT(void)
{
   register uint32 *p0=state, *p2=state+2, *pM=state+M, s0, s1;
   register int    j;

   if(left < -1)
      seedMT(4357U);

   left=N-1, next=state+1;

   for(s0=state[0], s1=state[1], j=N-M+1; --j; s0=s1, s1=*p2++)
      *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

   for(pM=state, j=M; --j; s0=s1, s1=*p2++)
      *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);

   s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K : 0U);
   s1 ^= (s1 >> 11);
   s1 ^= (s1 <<  7) & 0x9D2C5680U;
   s1 ^= (s1 << 15) & 0xEFC60000U;
   return(s1 ^ (s1 >> 18));
}


uint32
randomMT
(void)
{
   uint32 y;

   if(--left < 0)
      return(reloadMT());

   y  = *next++;
   y ^= (y >> 11);
   y ^= (y <<  7) & 0x9D2C5680U;
   y ^= (y << 15) & 0xEFC60000U;
   return(y ^ (y >> 18));
}

// Tiny function to return a uniform pseudo random number between 0
// and 1 using the Mersenne twister algorithm above.

double runifMT (void) { return randomMT() / 4294967295.0; }




/*
 *  Mathlib : A C Library of Special Functions
 *  Copyright (C) 1998 Ross Ihaka
 *  Copyright (C) 2000-2014 The R Core Team
 *  Copyright (C) 2007 The R Foundation
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/
 *
 *  DESCRIPTION
 *
 *	Random variates from the binomial distribution.
 *
 *  REFERENCE
 *
 *	Kachitvichyanukul, V. and Schmeiser, B. W. (1988).
 *	Binomial random variate generation.
 *	Communications of the ACM 31, 216-222.
 *	(Algorithm BTPEC).
 */


#define repeat for(;;)

size_t
rbinom
(
   size_t n,
   double pp
)
{

   static double c, fm, npq, p1, p2, p3, p4, qn;
   static double xl, xll, xlr, xm, xr;

   static double psave = -1.0;
   static int nsave = -1;
   static int m;

   double f, f1, f2, u, v, w, w2, x, x1, x2, z, z2;
   double p, q, np, g, r, al, alv, amaxp, ffm, ynorm;
   int i, ix, k;

   if (n == 0 || pp == 0.) return 0;
   if (pp == 1.) return n;

   p = pp < .5 ? pp : 1. - pp;
   q = 1. - p;
   np = n * p;
   r = p / q;
   g = r * (n + 1);

// Setup, perform only when parameters change using static (globals).

   if (pp != psave || n != nsave) {
      psave = pp;
      nsave = n;
      if (np < 30.0) {
         /* inverse cdf logic for mean less than 30 */
         qn = pow(q, n);
         goto L_np_small;
      } else {
         ffm = np + p;
         m = (int) ffm;
         fm = m;
         npq = np * q;
         p1 = (int)(2.195 * sqrt(npq) - 4.6 * q) + 0.5;
         xm = fm + 0.5;
         xl = xm - p1;
         xr = xm + p1;
         c = 0.134 + 20.5 / (15.3 + fm);
         al = (ffm - xl) / (ffm - xl * p);
         xll = al * (1.0 + 0.5 * al);
         al = (xr - ffm) / (xr * q);
         xlr = al * (1.0 + 0.5 * al);
         p2 = p1 * (1.0 + c + c);
         p3 = p2 + c / xll;
         p4 = p3 + c / xlr;
      }
   } else if (n == nsave) {
      if (np < 30.0)
         goto L_np_small;
   }

   //-------------------------- np = n*p >= 30 : ------------------- //
   repeat {
      u = runifMT() * p4;
      v = runifMT();
      /* triangular region */
      if (u <= p1) {
         ix = (int)(xm - p1 * v + u);
         goto finis;
      }
      /* parallelogram region */
      if (u <= p2) {
         x = xl + (u - p1) / c;
         v = v * c + 1.0 - fabs(xm - x) / p1;
         if (v > 1.0 || v <= 0.)
            continue;
         ix = (int) x;
      } else {
         if (u > p3) {	/* right tail */
            ix = (int)(xr - log(v) / xlr);
            if (ix > n)
               continue;
            v = v * (u - p3) * xlr;
         } else {/* left tail */
            ix = (int)(xl + log(v) / xll);
            if (ix < 0)
               continue;
            v = v * (u - p2) * xll;
         }
      }
      // determine appropriate way to perform accept/reject test //
      k = abs(ix - m);
      if (k <= 20 || k >= npq / 2 - 1) {
         /* explicit evaluation */
         f = 1.0;
         if (m < ix) {
            for (i = m + 1; i <= ix; i++)
               f *= (g / i - r);
         } else if (m != ix) {
            for (i = ix + 1; i <= m; i++)
               f /= (g / i - r);
         }
         if (v <= f)
            goto finis;
      } else {
         // squeezing using upper and lower bounds on log(f(x)) //
         amaxp = (k / npq) *
            ((k * (k / 3. + 0.625) + 0.1666666666666) / npq + 0.5);
         ynorm = -k * k / (2.0 * npq);
         alv = log(v);
         if (alv < ynorm - amaxp)
            goto finis;
         if (alv <= ynorm + amaxp) {
            // stirling's formula to machine accuracy //
            // for the final acceptance/rejection test //
            x1 = ix + 1;
            f1 = fm + 1.0;
            z = n + 1 - fm;
            w = n - ix + 1.0;
            z2 = z * z;
            x2 = x1 * x1;
            f2 = f1 * f1;
            w2 = w * w;
            if (alv <= xm * log(f1 / x1) + (n - m + 0.5) * log(z / w) + (ix - m) * log(w * p / (x1 * q)) + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) / f1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) / z / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) / x1 / 166320.0 + (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) / w / 166320.)
               goto finis;
         }
      }
   }

L_np_small:
   //---------------------- np = n*p < 30 : ------------------------- //

   repeat {
      ix = 0;
      f = qn;
      u = runifMT();
      repeat {
         if (u < f)
            goto finis;
         if (ix > 110)
            break;
         u -= f;
         ix++;
         f *= (g / ix - r);
      }
   }
finis:
   if (psave > 0.5)
      ix = n - ix;
   return ix;
}
