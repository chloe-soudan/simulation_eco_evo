 
#### SIMULATION ECO EVO  ###

## Part 1 

**Aims**: comparing stochastic and deterministic models 

**Tasks**:
- Run multiple stochastic simulations in parallel using HPC
- Figure of time series of the (log scaled) total population size comparing initial condition of 100 adult individuals (red lines) and 10 adult individuals (blue lines) with a deterministic version (solid line) and several example stochastic versions (dahsed line)
- probability of population extinction against initial population size under the stochastic model

## Diary entry 1: 18/11/2026 - 16:30

**Aims**

Running a deterministic and a stochastic model on a test data set to familiarise myself with the funtions function from demographic.r file. Here, I will use a small test data set that I generated myself (Step 1) to test and organise the code. 


**Step 1: Creating a test data set** 
- growth_matrix <- 5 by 5 matrix where sum of values in each column add up to 1 or less 

```r
growth_matrix <- matrix(c(
  0.4, 0.0, 0.0, 0.0, 0.0,
  0.2, 0.2, 0.2, 0.1, 0.1,
  0.1, 0.2, 0.2, 0.2, 0.1,
  0.1, 0.1, 0.2, 0.2, 0.1,
  0.1, 0.1, 0.1, 0.2, 0.1
), nrow = 5, byrow = TRUE)
```
- state <- 5 by 1 matrix with random values (between 1 and 100) representing the population size at each age/stage

```r
state<- matrix(
  sample(1:100, 5, replace = TRUE),
  nrow = 5,
  ncol = 1)
```
- reproduction_matrix <- 5 by 5 matrix where top right corner value indicates the average amount of offspring an adult at the oldest age/stage has

```r
reproduction_matrix <- matrix(c(
  0, 0, 0, 0, 5,
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0,
  0, 0, 0, 0, 0
), nrow = 5, byrow = TRUE)
```
- clutch_distribution <- numeric vector with probabilities of having clutch sizes of 1,2,3,...
```r
clutch_distribution <- c(0.1, 0.3, 0.4, 0.2)
```


**Step 2: Deterministic modelling - test data set**
Running deterministic modelling to simulate the population distribution and size over the next n (here 10) generations of the test data set (taking into account survival, maturation, and offspring produced). To do so I need to compute:

a) The next generation's population distribution (based on survival and maturation)

> Finding population distribution of next generation based on survival and maturation rates obtained from the growth matrix. I will use the deterministic_step() function from demographic.r but will name it 'deterministic_survival_maturation' as it computes the population distribution of the next generation based on the growth matrix and that I want my deterministic_step() function to also take into account offspring produced (see below).

> Inputs: (state, growth_matrix)  
> Expected outputs: 5 by 1 matrix with the population distribution (new_state) of the next generation 

b) The probability of recruitment 

> Finding the probability that an individual at the oldest/latest age/stage produces offspring (recruitment). I will use the stochastic_recruitment() function from demographic.r but will rename it to 'deterministic_recruitment'. I will assign the obtained probability to be stored as 'recruitment_probability' (```recruitment_probability <- deterministic_recruitment(reproduction_matrix, clutch_distribution) ```) to be able to use it in future steps.

> Inputs: (reproduction_matrix, clutch_distribution)  
> Expected outputs: a single probability value

c) The amount of offspring produced in a generation 
    
> Calculating the total number of offspring produced during a generation. I will use the offspring_calc() function from demographic.r BUT will remove the lines that are adding the random variation (stochastic modelling) and will replace them by lines adapted for deterministic modelling:

> - Instead of randomly determining the number of recruiting adults (```number_clutch = rbinom(1,  recruitment_probability) ```), I will set x to be the number of recruiting adutls
> - Instead of randomly drawing a clutch size (`clutch = random_draw(clutch_distribution)`), I will multiply clutch size probability by clutch size (```expected_clutch_size <- sum(
    c   lutch_distribution * (1:length(clutch_distribution)))```)
> - I will remove the 'if' statement and the 'for' loop
> - I will calculate the amount of offspring to be the number of recruiting adults(x) * the recruitment probability (from previous substep) * expected_clutch_size
  
> Adapted code:

```r
deterministic_offspring_calc  <- function(state, clutch_distribution, recruitment_probability) {
  offspring = 0 
  x = state[length(state)] # number of adults
  
  expected_clutch_size <- sum(
    clutch_distribution * (1:length(clutch_distribution))
  )
  
  offspring <- x * recruitment_probability * expected_clutch_size

  return(offspring) # total number of offspring
}
```

> Inputs: (state, clutch_distribution, recruitment_probability)   
> Excpected outputs: one single integer representing the sum of offspring produced in the generation

d) The next generation's population distribution based not only on survival, maturation BUT ALSO on offspring production

> Finding the population distribution of the next generation, taking into account survival, maturation, and offspring. I will use the stochastic_step() function from demographic.r and rename it 'deterministic_step'. !! I also need to make sure that all the function calls within the function refer to the deterministic modelling step (and not stochastic modelling).

> Inputs: (state, growth_matrix, reproduction_matrix, clutch_distribution, recruitment_probability)  
> Expected outputs: 5 by 1 matrix showing the computed poulation distribution (state) of the next generation

I will then put these substeps together and compute the population sizes of the next n-1 (here n = 10) generations (based on survival, maturation, and offspring production)

> Putting everything together and finding the population sizes of the next n generations, taking into account survival, maturation, and offspring. I will use the stochastic_simulation() function from demographic.r and rename it 'deterministic_simulation'.  However, the last lines from this function, which, when the population reaches extinction, is supposed to set the next population sizes and distribtution with 0s is wrong. I will thus replace it by:

```r
if (sum(population_size) == 0) {  the following gen with 0s
    population_size[(i + 2):(simulation_length + 1)] <- 0 }
```
> !! I also need to make sure that all the function calls within the function refer to the deterministic modelling step (and not stochastic modelling).

> Inputs: (state, growth_matrix, reproduction_matrix, clutch_distribution, n)  
> Expected outputs: n values, representing the simulated population sizes of the current population and the next n-1 generations

**Step 3: Stochastic modelling - test data set** 

Running stochastic modelling (incorporating random variation) to simulate the population distribution and size over the next n (here 10) generations of the test data set, taking into account survival, maturation, and offspring produced. 

Here, I will use similar functions as in Step 2, HOWEVER, I will add random variation (stochastic modelling). To do so I need to run a function that (1) checks that the sum of probabilities I am using is within bounds (not greater than 1) and (2) randomly distributes individuals across all ages/stages according to their respective probabilities and accounting for probability of death. This function comes from demographic.r (multinomial()).

Then, following Step 2 workflow, I found need to compute:

a) The next generation's population distribution (based on survival and maturation)

> Here, I will use the survival_maturation() function from demographic.r and rename it 'stochastic_survival_maturation'. However this function calls for another function (sum_vect()), which hasn't been defined yet. I thus will first define sum_vect as a function that sums up two imput vectors (```sum_vect <- function(x, y) {
  return(x + y)
}```
)

> Inputs: (state, growth_matrix)  
> Expected outputs: 5 by 1 matrix with the population distribution (new_state) of the next generation 

b) The probability of recruitment. 

> Here I will use the stochastic_recruitment() function from demographic.r. As in Step 2, I will set the outcome from this function to be named 'recruitment_probability'.  However, I first need to define the random_draw() function (from demogrpahic.r) that gets called in the stochastic_recruitment() function. This is a function that adds random variation to the model by randomly drawing an outcome from the input. 
> Note: there was an error in random_draw(): size was 11 but should be 1 here 

> Inputs for random_draw(): a numeric vector  
> Expected outputs from random_draw(): one single value from the input vector

> Inputs for stochastic_recruitment(): (reproduction_matrix, clutch_distribution)  
> Expected outputs from stochastic_recruitment(): a single probability value


c) The amount of offspring produced in a generation

> Here I will use the offspring_calc() function from demographic.r and rename it 'stochastic_offspring_calc'. The stochastic elements of this function are the following:
-  Number of reproducing adults is determined by using binomial distribution based on recruitment_probability
- Clutch size is determined by random draw (random_draw())of a clutch size of clutch_distribution.

> Inputs: (state, clutch_distribution, recruitment_probability)  
> Excpected outputs: one single integer representing the sum of offspring produced in the generation

d) The next generation's population distribution based not only on survival, maturation BUT ALSO on offspring production

> Here, I will use the stochastic_step() from demographic.r (same formula as deterministic_step)

> Inputs: (state, growth_matrix, reproduction_matrix, clutch_distribution, recruitment_probability)  
> Expected outputs: 5 by 1 matrix showing the computed poulation distribution (state) of the next generation

I can then put these substeps together and compute the population sizes of the next n (here n = 10) generations (based on survival, maturation, and offspring production). For this, I will use the stochastic_simulation() function from demographic.r. However, the last lines from this function, which, when the population reaches extinction, are supposed to set the next population sizes and distribtution to 0s was wrong. I thus will replace it by:
```r
if (sum(population_size) == 0) {  the following gen with 0s
    population_size[(i + 2):(simulation_length + 1)] <- 0 }
```
> Inputs: (state, growth_matrix, reproduction_matrix, clutch_distribution, n)  
> Expected outputs: n values, representing the simulated population sizes of the current population and the next n-1 generations

**Step 4: Comparing results from determinisitc and stochastic modelling**

I will store the results from both simulation functions (here I will simulate 10 generations) as 'deterministic_results' and 'stochastic_results' by running:

```r

# Deterministic simulation
deterministic_results <- deterministic_simulation(
  state,growth_matrix,reproduction_matrix,clutch_distribution, 10)

# Stochastic simulation
stochastic_results <- stochastic_simulation(state, growth_matrix,
  reproduction_matrix,clutch_distribution,10)
  ```
  
**Step 5: plotting the results**

I will use basic R plotting commands to plot the deterministic results againt the stochastic results, and will thus obtain a graph showing the results of both modelling approaches.