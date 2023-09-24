#include "sga_project_4_2023.h"	/* include the sga_project_4_2023 header file */

/*------------------------------------------------------*/
/*              RandVal                                 */
/* returns a pseudorandom value between the supplied    */
/* low and high arguments.                              */
/*                                                      */
double RandVal (double low, double high)
{
  return ((((double) rand())/RAND_MAX)*(high-low))+low;
}

/*-----------------------------------------------------------------*/
/*              INITIALIZE                                         */
/* The upper and lower bounds of each variable are set.            */
/* Then values between these bounds are randomly generated         */
/* and for each gene of each genotype in the population            */
/* a binary (0 or 1) value is assigned. Be aware that the          */
/* random number generator is not seeded using the system time.    */
/* The popoulation of the antigens is initialized in the same way. */
/*                                                                 */
void initialize(void)
{
  	int           i,j;

  	/* lower and upper bounds for each variable  */
  	/* we set the lower kai upper bounds to 0    */
  	/* and 1 for each variable-gene of each      */
  	/* antibody.                                 */
  	for(i=0;i<NVARS;i++)
    {
        for(j=0;j<=POPSIZE;j++)
        {
            Population[j].Lower[i]=0;
      		Population[j].Upper[i]=1;
    	}
  	}

  	/* Using the lower and upper bounds, randomly         */
  	/* generate a value using RandVal. If this value      */
  	/* is >= 0.5 then asign 1 to this gene, else 0.       */
  	/* Repeat this for each gene of each genotype.        */
  	for(i=0;i<=POPSIZE;i++)
  	{
    	for(j=0;j<NVARS;j++)
    	{
      		 if (RandVal(Population[i].Lower[j], Population[i].Upper[j])>=0.5)
             Population[i].Gene[j] = 1;
             else
             Population[i].Gene[j] = 0;
        }
    }

    /* lower and upper bounds for each variable  */
    /* we set the lower kai upper bounds to 0    */
  	/* and 1 for each variable of each antigen.  */
  	for(i=0;i<NVARS;i++)
    {
        for(j=0;j<ANTIGENPOPSIZE;j++)
        {
            Antigens[j].Lower[i]=0;
      		Antigens[j].Upper[i]=1;
    	}
    }
    
    /* Using the lower and upper bounds, randomly         */
  	/* generate a value using RandVal. If this value      */
  	/* is >= 0.5 then asign 1 to this variable, else 0.   */
  	/* Repeat this for each variable of each antigen.     */
  	for(i=0;i<ANTIGENPOPSIZE;i++)
  	{
    	for(j=0;j<NVARS;j++)
    	{
      		if (RandVal(Antigens[i].Lower[j], Antigens[i].Upper[j])>=0.5)
      		Antigens[i].Gene[j] = 1;
      		else
      		Antigens[i].Gene[j] = 0;
        }
    }
      		    
}//end of initialize

/*--------------------------------------------------------------*/
/*              EVALUATE                                        */
/* is a user defined function. In our problem the fitness       */
/* function is the number of 1 that come from complementary     */
/* genes between an antibody and an antigen for all the         */
/* existing antigens.                                           */
/*                                                              */
void evaluate(void)
{
	int mem=0,i=0,k=0,m=0,fitness;
  	Best=0;
  	/* We compute the fitness function of each genotype of the  */
  	/* population.*/
  	for(mem=0; mem<POPSIZE; mem++)
    {
               Population[mem].Fitness=0;
               /* For each genotype and for each antigen we compute the */
               /* varible named fitness. We then add it to the fitness  */
               /* of the genotype.*/
               for(k=0; k<ANTIGENPOPSIZE; k++)
               {
                        fitness=0;
                        /* Fitness is computed by adding the ones (1) that come */
                        /* from complementary genes between an antibody and     */
                        /* and an antigen.*/
                        for(m=0; m<NVARS; m++)
                        {
                                 complementaries[mem].Gene[m]=0;
                                 if (Population[mem].Gene[m]==1)
                                 {
                                     if (Antigens[k].Gene[m]==1)
                                     {
                                         complementaries[mem].Gene[m]=0;
                                     }
                                     if (Antigens[k].Gene[m]==0)
                                     {
                                         complementaries[mem].Gene[m]=1;
                                     }
                                 }//telos if
                                 if (Population[mem].Gene[m]==0)
                                 {
                                     if (Antigens[k].Gene[m]==1)
                                     {
                                         complementaries[mem].Gene[m]=1;
                                     }
                                     if (Antigens[k].Gene[m]==0)
                                     {
                                         complementaries[mem].Gene[m]=0;
                                     }
                                 }//telos if
                                 fitness+=complementaries[mem].Gene[m];
                        }//telos for m
                        Population[mem].Fitness+=fitness;
               }//telos for k
    

    	/* Keep track of the best member of the population  */
    	/* Note that the last member of the population holds*/
    	/* a copy of the best member.                       */

    	if(Population[mem].Fitness>Population[POPSIZE].Fitness){
      		Best=mem;
      		Population[POPSIZE].Fitness=Population[mem].Fitness;
      		for(i=0;i<NVARS;i++)
        		Population[POPSIZE].Gene[i]=Population[mem].Gene[i];
        }//telos if
  	}//telos for mem
}//end of evaluate

/*---------------------------------------------------------------*/
/*              COPY_GENOTYPES                                   */
/* 		Coppies a genotype to another                    */
/*                                                               */
void copy_genotypes(struct genotype* oldgenotype, struct genotype* newgenotype)
{
	int i=0; /* some counter variables */

  	for(i=0; i<NVARS; i++)
        	newgenotype->Gene[i]=oldgenotype->Gene[i];

  	newgenotype->Fitness=oldgenotype->Fitness;

  	for(i=0; i<NVARS; i++)
        	newgenotype->Upper[i]=oldgenotype->Upper[i];
  	for(i=0; i<NVARS; i++)
        	newgenotype->Lower[i]=oldgenotype->Lower[i];

  	newgenotype->RFitness=oldgenotype->RFitness;
  	newgenotype->CFitness=oldgenotype->CFitness;
  	newgenotype->Selection_Prob=oldgenotype->Selection_Prob;
  	newgenotype->Cumulative_Prob=oldgenotype->Cumulative_Prob;

  	newgenotype->Survivor=oldgenotype->Survivor;
  	newgenotype->Mate=oldgenotype->Mate;
  	newgenotype->Mutate=oldgenotype->Mutate;
}

/*--------------------------------------------------------------*/
/*              COPY_POPULATION                                 */
/* 		Copies a population to another population       */
/*                                                              */
void copy_population(struct genotype old_pop[POPSIZE+1],struct genotype new_pop[POPSIZE+1])
{
	int mem=0; /* some counter variables */

  	for(mem=0; mem<=POPSIZE; mem++)
    	copy_genotypes(&old_pop[mem],&new_pop[mem]);
}

/*--------------------------------------------------------------*/
/*              SELECT                                          */
/* This is an implementation of STANDARD PROPORTIONAL SELECTION */
/* (or ROULETTE WHEEL SELECTION) for MAXIMIZATION problems      */
/* It also checks to make sure that the best member survives    */
/* (i.e., elitest selection).                                   */
/*                                                              */
void select(void)
{
  	int    member, spin_num, mem; /* Some counter variables       */
  	double Total_Fitness;       /* The total population fitness */
  	double roulette=0;
  	/* a variable for temporary storing of the population */
  	struct genotype Buffered_Pop[POPSIZE+1];
  	int Found;                    /* A flag */

  	/* First, find the total fitness of the population    */
  	Total_Fitness=0;
  	for (member=0; member<POPSIZE; member++)
      	Total_Fitness += Population[member].Fitness;

  	/* Next, calculate the probability of a selection of each genotype      */
  	for(member=0; member<POPSIZE; member++)
    	Population[member].Selection_Prob = Population[member].Fitness/Total_Fitness;

	/* Now, calculate the cumulative probability of each genotype     */
  	Population[0].Cumulative_Prob=Population[0].Selection_Prob;

  	for(member=1; member<POPSIZE; member++)
    	Population[member].Cumulative_Prob = Population[member-1].Cumulative_Prob +
                                             Population[member].Selection_Prob;

  	/* Finally, select the survivors using the cumulative probability */
  	/* and create the new Population                                  */
  	for(spin_num=0; spin_num<POPSIZE; spin_num++) {
    	Found = FALSE;
    	roulette = (double) rand()/RAND_MAX;

    	if(roulette < Population[0].Cumulative_Prob) 			/* Select the first member of the Population */
      		copy_genotypes(&Population[0],&Buffered_Pop[spin_num]);
        else if(roulette >Population[POPSIZE-1].Cumulative_Prob) /* select the best member of the population */
             copy_genotypes(&Population[POPSIZE],&Buffered_Pop[spin_num]);
        else
        	for(mem=1; mem<POPSIZE && Found==FALSE; mem++)
           		if( (roulette > Population[mem-1].Cumulative_Prob) &&
              	    (roulette <=Population[mem].Cumulative_Prob) )
           		{
              		copy_genotypes(&Population[mem],&Buffered_Pop[spin_num]);
              		Found = TRUE;
           		}
	}

  	/* copy the best individual */
  	copy_genotypes(&Population[POPSIZE],&Buffered_Pop[POPSIZE]);

  	/* Copy the Buffered_Pop to the original Population */
  	copy_population(Buffered_Pop,Population);

  	/* Population , now is the new population           */
}

/*-----------------------------------------------------------------*/
/*              CROSSOVER                                          */
/* This is an implementation of STANDARD SINGLE POINT CROSSOVER.   */
/* Many other crossover operators developed specifically for       */
/* real-coded GA's may give better results. For simplicity only    */
/* the single point crossover is shown here.                       */
/*                                                                 */
void crossover(void)
{
  	int mem,i,            /* Counting variables   */
      	parent1,          /* Parent one           */
      	parent2,          /* Parent two           */
      	xover_point,      /* Crossover point      */
      	count=0,
      	lovers=0,         /* number of matting genotypes */
      	indiv=0;
  	struct genotype temp_parent;

  	for(mem=0; mem<=POPSIZE; mem++)
    	Population[mem].Mate=FALSE;

  	/* first find the individuals for the Crossover operation */
  	for (mem=0;mem<POPSIZE;mem++)
    	if (frand() < PXOVER){ /* frand returns a random number in the range [0,1] */
       		Population[mem].Mate = TRUE;
       		lovers++;
    	}

	/* We want an even number of "lovers"*/
  	if((lovers%2) != 0) {
    	do  	/* find an non "lover" */
    		indiv=rand()%POPSIZE;
    	while(Population[indiv].Mate==TRUE);
  		/* make it "lover"    */
  		Population[indiv].Mate=TRUE;
  		lovers++;
	}

  	/* second mate the "lovers" */
  	mem=0;
  	for(count=0; count<(lovers/2); count++) {
    	while(Population[mem].Mate==FALSE) mem++; /* find the first parent */
    	parent1=mem;
		mem++;
    	while(Population[mem].Mate==FALSE) mem++; /* find the second parent */
    	parent2=mem;
		mem++;

    	/* select the crossover point :1...NVARS-1 */
    	xover_point=(rand()%(NVARS-1))+1;

    	/* Perform the crossover */
    	/* copy parent1 to temp_parent */
    	copy_genotypes(&Population[parent1],&temp_parent);

    	for(i=xover_point; i<NVARS; i++)
        	Population[parent1].Gene[i]=Population[parent2].Gene[i];
    	for(i=xover_point;i<NVARS;i++)
        	Population[parent2].Gene[i]=temp_parent.Gene[i];
  	}
  	/* set Mate flag to False */
	for(mem=0;mem<=POPSIZE;mem++)
    	Population[mem].Mate=FALSE;
}

/*-----------------------------------------------------------------*/
/*              MUTATION                                           */
/* Because the genes of a genotype have binary values, a value     */
/* selected for mutation, based on the mutation probability,       */
/* is replaced by its compementary value.                          */
/*                                                                 */
void mutate(void)
{
  	double lbound,hbound;
  	int member,             /* The member to be mutated                 */
      	var;                /* The variable to be mutated               */


  	for(member=0; member<POPSIZE; member++) /* for every member */
    	for(var=0; var<NVARS; var++)		/* for every gene   */
      		if(frand()<PMUTATION)
            {
         		/* Change the value of the gene to its complementary value. */
         		if (Population[member].Gene[var]==0)
                Population[member].Gene[var]= 1;
                else
                Population[member].Gene[var]=0;
       		}
}

/*--------------------------------------------------------------*/
/*              REPORT                                          */
/* Report progress of the simulation. Data is dumped to a log   */
/* file in comma separeted value format which can be imported   */
/* and graphed using any commercial spreadsheet package.        */
/*                                                              */
void report(void)
{
  	double best_val,      /* Best population fitness                  */
           avg,           /* Average population fitness                   */
           stddev,        /* Std. deviation of population fitness         */
           sum_square,    /* Sum of the squares for std. dev calc.        */
           square_sum,    /* Square of the sums for std. dev. calc.       */
           sum;           /* The summed population fitness                */
	int i=0;	/* counter */

	sum=0.0;
  	sum_square=0.0;

  	/* Calculate the summed population fitness and the sum        */
  	/* of the squared individual fitnesses.                       */
  	for(i=0; i<POPSIZE; i++){
    	sum += (Population[i].Fitness);
    	sum_square += pow(Population[i].Fitness, 2);
	}

  	/* Calculate the average and standard deviations of the       */
  	/* population's fitness.                                      */
  	avg=sum/(double)POPSIZE;
  	square_sum=sum*sum/(double) POPSIZE;
  	stddev=sqrt((1.0/(double)(POPSIZE-1))*(sum_square-square_sum));
  	best_val=Population[POPSIZE].Fitness;

  	/* Print the generation, best, avg, and std. deviation to a   */
  	/* file in csv format.                                        */
  	fprintf(galog,"\n%10d  %15.4f  %15.4f  %15.4f", Generation, best_val, avg, stddev);
	/* Print the Best Genotype */
	fprintf(galog, "  (");
	for(i=0;i<NVARS;i++)
    	fprintf(galog," %5.3f ",Population[POPSIZE].Gene[i]);
	fprintf(galog, ") ");
}

/*--------------------------------------------------------------*/
/*                     MAIN                                     */
/* This is the main function. It loops for the specified number */
/* of generations. Each generation involves selecting survivors,*/
/* performing crossover and mutation, and then evaluating the   */
/* resulting population.                                        */
/*                                                              */
void main(void)
{
	int i,j;

  	srand(time(0));

  	if((galog=fopen("galog.txt","w"))==NULL)
    {
		printf("Can't open galog.txt \n");
    	exit(1);
  	}
	fprintf(galog,"%10s  %15s  %15s  %15s  %20s\n", "Generation", "Best Value", "Average", "StdDev", "Best Genotype");

  	Generation=0;
  	initialize();
  	/* We print the initial population of antibodies. */
  	printf("Initial Population\n");
	for (j=0; j<NVARS; j++)
		printf("  Gene%d   ", j);
	printf("\n");
	for(i=0;i<POPSIZE;i++)
    {
    	for(j=0;j<NVARS;j++)
        	printf("%10.5f", Population[i].Gene[j]);
    	printf("\n");
  	}
  	/* We print the antigens that we have created.    */
  	printf("Antigens\n");
	for (j=0; j<NVARS; j++)
		printf("  Gene%d   ", j);
	printf("\n");
	for(i=0;i<ANTIGENPOPSIZE;i++)
    {
    	for(j=0;j<NVARS;j++)
        	printf("%10.5f", Antigens[i].Gene[j]);
    	printf("\n");
  	}
  	printf("Press any key to start the evolution\n");
	getch();

  	evaluate();
  	while(Generation<MAXGENS)
    {
    	Generation++;
    	select();
    	crossover();
    	mutate();
    	evaluate();
    	printf("Generation : %d \n",Generation);
    	report();
  	}
	/* print to screen */
	printf("\n\nSimulation completed\n");
  	printf("\n   Best Member:");
	for(i=0; i<NVARS; i++)
    	printf(" %f ",Population[POPSIZE].Gene[i]);
	printf("  Fitness: %5f\n", Population[POPSIZE].Fitness);

    /* print to the log file*/
	fprintf(galog,"\n\nSimulation completed\n");
  	fprintf(galog,"\n   Best member :\n");

  	for(i=0;i<NVARS;i++)
    	fprintf(galog,"\n Var(%d) = %3.3f",i,Population[POPSIZE].Gene[i]);
	fprintf(galog, "  Fitness: %5f\n", Population[POPSIZE].Fitness);
	printf("Final Population\n");
	for (j=0; j<NVARS; j++)
		printf("  Gene%d   ", j);
	printf("\n");
	for(i=0;i<POPSIZE;i++)
    {
    	for(j=0;j<NVARS;j++)
        	printf("%10.5f", Population[i].Gene[j]);
    	printf("\n");
  	}

  	fclose(galog);
  	system("PAUSE");
}