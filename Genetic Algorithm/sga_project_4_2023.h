/*--------------------------------------------------------------------------------*
 * A Simple Genetic Algorithm (SGA), based on the algorithm presented in the book: *
 * Michalewicz, Z., "Genetic Algorithms + Data Structures = Evolution Programs",  *
 * Springer-Verlag, 3rd edition, 1996.                                            *
 *--------------------------------------------------------------------------------*/

/*--------------- SGA.H  Header File ---------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <time.h>

#define frand_ab() (2*((double) rand() / RAND_MAX)-1)
#define frand() ((double) rand()/RAND_MAX)

/* change any of these parameters to match your needs */

#define POPSIZE         10      /* Population size                  */
#define MAXGENS         50     /* Maximum number of generations    */
#define NVARS           32        /* Number of problem variables      */
#define PXOVER          0.9     /* Probability of crossover         */
#define PMUTATION       0.1      /* Probability of mutation          */
#define TRUE            1
#define FALSE           0
#define ANTIGENPOPSIZE  1000    /* Number of antigens               */

int Generation,                  /* Curent Generation number                */
    Best;                        /* The Best genotype in the population     */
FILE *galog;                     /* A log file to write results to          */

struct genotype                 /* Each genotype is a member of           */
{                               /* the population                         */
   double Gene[NVARS], 		/* A string of variables makes a genotype */
	  Fitness,              /* The genotype's fitness                 */
	  Upper[NVARS],	        /* The genotype's upper bound             */
	  Lower[NVARS],         /* The genotype's lower bound             */
	  RFitness,             /* The relative fitness                   */
	  CFitness,             /* The comulative fitness                 */
	  Selection_Prob,       /* The probability of a selection         */
	  Cumulative_Prob;      /* The cumulative probability             */
   int Survivor, 	      	/* Flag for selection procedure           */
       Mate,         	 	/* Flag for Crossover procedure           */
       Mutate;        		/* Flag for Mutation procedure            */
};

struct genotype Population[POPSIZE+1]; /* The population of genotypes */
struct genotype Best_Individual;

/*the struct for an antigen*/
struct antigen
{
   double Gene[NVARS], 	         /* A string of variables makes a genotype */
          Upper[NVARS],          /* The genotype's upper bound             */
          Lower[NVARS];          /* The genotype's lower bound             */
};
 
struct antigen Antigens[ANTIGENPOPSIZE]; /* The array of antigens */
    
/* the struct for the complementaries between antibodies and antigens*/
struct complementary
{
   double Gene[NVARS];
};

struct complementary complementaries[POPSIZE];

/* GA function prototypes */
void initialize(void);
double RandVal(double,double);
void evaluate(void);
void copy_genotypes(struct genotype*, struct genotype*);
void copy_population(struct genotype old_pop[POPSIZE+1], struct genotype new_pop[POPSIZE+1]);
void select(void);
void crossover(void);
void mutate(void);
void report(void);
