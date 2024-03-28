#include <Rcpp.h> 
using namespace Rcpp;

#include <stdio.h>
#include <math.h>
#include <time.h>  /* for seeding the random number generator */

#define RAND01_PRECISION 100000 /* granularity for generating random numbers */
#define EPSILON 1.0e-5 /* used to test whether or not a number is zero */

typedef unsigned char byte;  /* just for convenience */

// [[Rcpp::export]]
int
  randFrac(double frac)
  {
    int i;
    double d;
    
    /* On the system I use, the random() function returns a random
     * 32-bit integer.  If you are using a computer which has a
     * function that returns a random floating-point number between
     * 0 and 1, then you can just use that function to set the value
     * of the variable "d" directly, which is basically what I achieve
     * below.
     */
    i = rand() % RAND01_PRECISION;
    d = ((double) i) / (double)RAND01_PRECISION;
    if (d < frac)
      return 1;
    else
      return 0;
  }

// [[Rcpp::export]]
double unirand()
{
  return (rand()+1.0)/(RAND_MAX+1.0);
}


/*
* Check to see whether a number is close enough to zero that we really
* should call it zero.
*/
// [[Rcpp::export]]
int
  fis_zero(double x)
  {
    if (fabs(x) > EPSILON)
      return 0;
    else
      return 1;
  }

void
  measureBlockCounts(int counts[4],  /* array for returning the counts */
byte *lattice,  /* pointer to the landscape array */
int xsize, int ysize)  /* size of the landscape array */
  {
    int i, x, y, x1, y1;
    
    /* First initialize the counts to zero */
    for (i=0; i < 4; i++){
      counts[i] = 0;
    }
    
    for (y=0; y < ysize; y++){
      for (x=0; x < xsize; x++) {
        /* north */
        x1 = x;
        y1 = (y-1+ysize)%ysize;
        counts[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]]++;
        
        /* south */
        x1 = x;
        y1 = (y+1)%ysize;
        counts[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]]++;
        
        /* east */
        x1 = (x+1)%xsize;
        y1 = y;
        counts[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]]++;
        
        /* west */
        x1 = (x-1+xsize)%xsize;
        y1 = y;
        counts[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]]++;
      }
    }
  }

/*
 * Compute how the block counts would be affected if we were
 * to flip the state of the cell at position (x,y)
 */
void
  computeCountChanges(int changes[4], /* array to return the vector of changes */
byte *lattice,  /* landscape array */
int xsize, int ysize, /* size of landscape */
int x, int y)  /* which site we want to flip */
  {
    int x1, y1, i;
    
    /* Clear out array just to be safe */
    for (i=0; i < 4; i++)
      changes[i] = 0;
    
    /* north */
    x1 = x;
    y1 = (y-1+ysize)%ysize;
    /* First decrease the appropriate count for the current configuration.
    * Note we always increase or decrease by 2, because every time we make
    * a flip, it will change the 2x1 block as seen from both ends, so to
    * speak, so we need to count the change twice. */
    changes[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]] -= 2;
    /* And then increase the count for the new configuration */
    changes[(1-lattice[y*xsize+x])*2 + lattice[y1*xsize+x1]] += 2;
    
    /* south */
    x1 = x;
    y1 = (y+1)%ysize;
    changes[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]] -= 2;
    changes[(1-lattice[y*xsize+x])*2 + lattice[y1*xsize+x1]] += 2;
    
    /* east */
    x1 = (x+1)%xsize;
    y1 = y;
    changes[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]] -= 2;
    changes[(1-lattice[y*xsize+x])*2 + lattice[y1*xsize+x1]] += 2;
    
    /* west */
    x1 = (x-1+xsize)%xsize;
    y1 = y;
    changes[lattice[y*xsize+x]*2 + lattice[y1*xsize+x1]] -= 2;
    changes[(1-lattice[y*xsize+x])*2 + lattice[y1*xsize+x1]] += 2;
    
    /* Now handle symmetry, lumping together the changes for [01] and [10]
    * blocks, since they are really the same by symmetry assumptions. */
    i = changes[1] + changes[2];
   
    changes[1] = changes[2] = i/2;
  }

/*
 * Given the parameters p0 and q00, compute the 2x1 block
 * probabilities, p00, p01, and p11.
 */
void
  condProbsTo2x1Probs(double p0, double q00,  /* input parameters */
double *p00, double *p01, double *p11)  /* return parms */
  {
    *p00 = p0 * q00;				/* p00 */
*p01 = p0 - *p00;				/* p01 */
*p11 = 1.0 - (*p00 + 2.0 * (*p01));		/* p11 */

/* Sometimes we end up with a parameter which is 1.0e-10, i.e. zero for
 * all practical purposes.  Here we set such values to actually BE zero,
 * for simplicity */
if (fis_zero(*p00))
  *p00 = 0.0;
if (fis_zero(*p01))
  *p01 = 0.0;
if (fis_zero(*p11))
  *p11 = 0.0;
return;
  }

/*
 * This is the heart of the program.  It generates the landscape with
 * spatially structured heterogeneities.
 */
void
  gen2x1(byte *lattice,  /* landscape array */
int xsize, int ysize,  /* size of landscape */
double p0, double q00,  /* the main parameter values */
int iterates, /*maximum number of times to loop through the algorithm*/
double maxErr)  /* the maximum difference between measured and desired
  * probabilities that we will accept when stopping */
  {
    int
    x, y, x1, y1, i, j,
    blockCounts[4],  /* used to count the numbers of occurances of the
               * various possible 2x1 blocks on the landscape */
tmpCounts[4],  /* a temporary array of counts */
targetCounts[4],  /* the counts we would see if we had the landscape
            * we are trying to get */
countChanges[4],  /* for computing changes to counts */
oldErr, newErr,  /* differences between desired and achieved COUNTS */
timesFlipped,  /* for counting how many times we flipped sites */
stop;  /* flag for stopping the algorithm */
double p00, p01, p11,  /* the 2x1 block probabilities */
p1,
probErr; /* differences between desired and
* achieved PROBABILITIES */

timesFlipped = 0;

/* Compute the 2x1 block probabilities */
condProbsTo2x1Probs(p0, q00, &p00, &p01, &p11);
p1 = 1.0 - p0;

/* First fill the lattice randomly, using the correct patch
* occupancy probability */
//bzero(lattice, xsize*ysize);  /* Clear out the lattice */
for (y=0; y < ysize; y++){
for (x=0; x < xsize; x++){
    if (randFrac(p1))
      lattice[y*xsize + x] = 1;
}
}
    
    /* Now measure the block counts from the lattice */
    measureBlockCounts(blockCounts, lattice, xsize, ysize);
    /* Compute target block counts from the block probabilities.
    * The reason for the factor of 4 is because when we count blocks,
    * we check 4 neighbors of each site.
    */
    targetCounts[0] = (int)(p00 * xsize*ysize*4);
    targetCounts[1] = targetCounts[2] = (int)(p01 * xsize*ysize*4);
    targetCounts[3] = (int)(p11 * xsize*ysize*4);
    
    oldErr = 0;
    
    /* Initialize the discrepancy between desired and achieved counts */
    for (i=0; i < 4; i++)
      oldErr += abs(targetCounts[i] - blockCounts[i]);
    
    /* Now comes the main loop.  We repeatedly pick a site at random, and
    * see if flipping its value (0 to 1, or 1 to 0) would move the block
    * counts (or equivalently, the 2x1 block probabilities) closer to the
    * desired values.  If so, flip the site.  Otherwise, try again. */
    stop = 0;
    for (i=0; (i < iterates) && (stop == 0); i++) {
      for (j=0; j < 4; j++)  /* copy the current counts into a temp. array */
    tmpCounts[j] = blockCounts[j];
      
      /* Pick a random site.
      * If your random number generator returns a number between 0 and 1,
      * you will need to do something like this [say your random number
      * generator is called rand01() ]:
      *   x = (int)(rand01() * xsize);
      *   y = (int)(rand01() * ysize);
      * That assumes rand01() returns a number greater than or equal
      * to zero, but strictly less than 1.
      */
      x = rand()%xsize;
      y = rand()%ysize;
      
      /* now see how the counts would be affected if we were to
      * change the state of this cell */
      computeCountChanges(countChanges, lattice, xsize, ysize, x, y);
      
      /* Merge the changes back into the count array */
      for (j=0; j < 4; j++)
        tmpCounts[j] += countChanges[j];
      
      /* And compute what the new discrepancy between desired and achieved
      * counts would be if we make this flip */
      newErr = 0;
      for (j=0; j < 4; j++)
        newErr += abs(targetCounts[j] - tmpCounts[j]);
      
      /* if the new error would be smaller, flip the site's value */
      if (newErr <= oldErr) {
        lattice[y*xsize+x] = 1-lattice[y*xsize+x];
        oldErr = newErr;
        /* Move the temporary copy of the counts into the regular array */
        for (j=0; j < 4; j++)
          blockCounts[j] = tmpCounts[j];
        timesFlipped++;  /* keep track of how many flips we make */
      }
      
      /* Convert the discrepancy in number of blocks into discrepancy
      * in probabilities, by normalizing */
      probErr = ((double)oldErr)/(double)(xsize*ysize*4.0);
      
      /* And if the error is small enough, stop. */
      if (probErr < maxErr)
        stop = 1;
    }
    
  }

// [[Rcpp::export]]
IntegerVector CreateLandscape(int xsize, int ysize,
                              double p0,double q00,int iterates, double maxErr) 
{
  IntegerVector paisaje(xsize*ysize);
  //int iterates=DEFAULT_MAXITERS,  /* maximum number of iterations */
  int x, y, i,
         blockCounts[4];
  double p00, p01, p11,  /* the 2x1 block probabilities */
      dsum, actualProbs[4], actualp0, actualq00;
  byte *lattice;  /* the main landscape array */
  
  /* Compute 2x1 block probabilities */
  condProbsTo2x1Probs(p0, q00, &p00, &p01, &p11);
  
  lattice=(byte *)malloc(xsize*ysize);

  /* Now generate the landscape. */
  gen2x1(lattice, xsize, ysize, p0, q00, iterates, maxErr);
    
    /* And output it. */
    /* First I output a standard header, which my simulation program
    * will recognize, to make sure I use the right kind of landscape files
    * with the simulations... */
    
    for (y=0; y < ysize; y++) {
      for (x=0; x < xsize; x++)
        paisaje[y*xsize+x]=1-lattice[y*xsize+x];
    }

    /* Now compute some diagnostics */
    measureBlockCounts(blockCounts, lattice, xsize, ysize);
    
    /* Compute the 2x1 block probabilities in the landscape we generated */
    for (i=0; i < 4; i++)
      actualProbs[i] = ((double)blockCounts[i])/(xsize*ysize*4.0);
    actualp0 = actualProbs[0] + actualProbs[1];
    if (actualp0 == 0.0)
      actualq00 = 0.0;
    else
      actualq00 = actualProbs[0]/actualp0;  /*  q00 = p[00]/p0  */
    
    dsum = fabs(actualProbs[0] - p00) + fabs(actualProbs[1] - p01) +
      fabs(actualProbs[2] - p01) + fabs(actualProbs[3] - p11);
    
    return(paisaje);
}


