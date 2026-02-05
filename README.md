# Neuroimmunological methods in stroke research course

## Instructions 

## Download and install R and Rstudio

Download and install the specific version of R for your OS [here:](https://ftp.fau.de/cran/)
(If using Windows 10/11, please install [Rtools](https://cran.r-project.org/bin/windows/Rtools/rtools45/files/rtools45-6608-6492.exe) as well) 

Download and install [Rstudio](https://posit.co/download/rstudio-desktop/)

### Before we begin 

**Note that everything below is subject to change before the course begins on 13/02/26, please install the necessary R packages before the day of the course if possible. 

Download the course content and extract the folder. 

![Alt text](img/download_instructions.gif)  [](img/download_instructions.gif)

Navigate to the course folder (Neuroimmunological-methods-in-stroke-research-course-2026) and open the "install.R" file in Rstudio. Run the code by selecting all of the code and pressing run. 

![Alt text](img/installation_instructions.gif)  [](img/installation_instructions.gif)

If everything installed successfully, you should receive the following message: 

"Installed packages successfully :)" 

Please check everything works before the course starts.

## Course content 

## Intro-to-R 

Since some of you may not have any experience using R, we will very quickly go through some key concepts to help you understand better what we will doing in the analysis section later. 

### Mathematical operations 

Type and test these operations in the console in Rstudio 

Addition
```1 + 1```

Subtraction 
```5 - 4``` 

Multiplication 
```2 * 10```

Division
```9 / 3```

Modulo 
```20 %/% 6```

Square
```4^2 or 4**2``` 

### Comparison

Greater than
```>``` 

Less than 
```<``` 

Equals 
```==```

Does not equal 
```!=```

### Basic types 

There are 6 basic data types in R: logical, numeric, integer, character, raw and complex. Today (and in general) only the first four are important. 

```logical``` Can only have two values ```TRUE``` or ```FALSE``` (Shorthand:```T```,```F```). 

```numeric``` All numbers with or without decimal ```100.34```

```integer``` Whole numbers ```123L``` A number appended with L suffix denotes an integer in R

```character``` Character or string values. Enclosed in either single 'microbe' or double "stroke" quotes. 

### Variable assignment 

Variables are assigned in R using the ```<-``` operator. A variable can be though of a a "box" you store a value in for later use.

They can be used to store a value which you want to be constant throughout your code e.g. a file name 
``` r
x <- "Bacteria_data.txt" 
``` 
Or to save time by storing long strings/numbers 

``` r
y <- 127575738398292929
```

### Data structures 

Four most commonly used data structures 

#### 1D
Vectors - ordered collection of data of a given length and all of the same type:
``` r
vector1 <- c(1,3,7,8,9)
vector2 <- c("A", "B", "C", "D")
print(vector2)
```

Lists - ordered collection of data/objects of a given length. Can be composed of heterogenous types. 
``` r
list1 <- list(vector1, vector2)
list2 <- list(1,2,3)
print(list1)
```

#### 2D 

Matrices - 2D data structure (rows x columns) of the same type
``` r
m <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
print(m)
```

Dataframes - 2D heterogenous data structure (rows x columns)
``` r
df <- data.frame(bacteria = c("E. coli", "L. reuteri", "C. difficile"),
                 sample1 = c(10, 27, 61),
                 sample2 = c(9, 200, 43) )
print(df) 
```

Objects such as vectors and dataframes can also be stored in variables as you saw above. ```df``` and ```vector2``` are both variables. 

### Control flow

As well as writing individual statements like we have done above we can also use logic to control the flow of our code, allowing us to repeat bits of code or only run if a given condition is met. 

for loops allow us to repeat code a specified number of times e.g.

``` r
for (i in vector2){
  print(i)
 }
```

This prints out each element in the ```vector2```variable we defined earlier. Not particularly interesting though.. 
We will discover a more interesting use case later. 

if/else statements allow us to control the flow of our code better: 

``` r
x <- 100
y <- 10

if (x > y) {
  print("x is higher")
} else if (x < y) {
  print("y is higher")
}
```
Change the values of x and y and see how the output changes. 

### Functions 

Functions are used to abstract repetitive code in something which is reusable and modular, making code easier to read and debug. 

Using the above if else example we could create a function called ```greater_than``` which tells us which value is highest in a more modular way

``` r
greater_than <- function(x, y) {
  if (x > y) {
    print("x is higher")
  } else if (x < y) {
    print("y is higher")
  }
}
```

``` r
greater_than(x=100, y=10)
greater_than(x=1, y=17)
```

More interesting functions might perform a calculation for us and return the value 

``` r
normalise <- function(x) {
  x_norm <- t(100 * t(x) / colSums(x))
  return(x_norm)
}
```
Let's run this function on the matrix we created earlier to see what it does. 

``` r
normalise(m)
```

So far we have covered built-in functions and custom functions, but R has a huge 
open source library of packages called CRAN, as well as a specific repository
for bioinformatics, called bioconductor. 

### Loading external functions

When working with R, at some point you may find yourself requiring a specific 
function to perform some kind of analysis and plot data etc. It's usually the 
case that someone has already written a function for whatever you want to do

You can install packages to access functions written by others. 
External packages can be downloaded using the ```install.packages``` function 
and loaded using the ```library``` function. Somewhat annoyingly, bioconductor
packages are installed a little differently, and require you to first install
the "BiocManager" package and then install using ```biocManager::install```. They
can however be loaded as normal. 

If you need help with a function you can also type ?functionname in the console e.g. ?log10 and the help for that function will show up, detailing what the function does, what inputs it expects and what value(s) it returns. 

Once you get into writing your own functions, it's good practice to store them
in a separate R file and import them using the ```source``` function. 

### Data wrangling 

A few key concepts on loading and manipulating data. 

Reading data 
``` r
# load tidyverse
library(tidyverse)
# read in the palmer penguins csv file 
penguins <- read_csv("data/penguins.csv")
```

Manipulating data 

Selecting columns 

``` r
# select certain columns 
bill_length <- select(penguins, bill_length_mm) 
print(bill_length)
# using pipe operator
bill_length <- penguins %>% 
    select(bill_length_mm)
print(bill_length)
```

Filtering rows 

``` r
gentoo  <- filter(penguins, species == "Gentoo")
print(gentoo)
# with pipe 
gentoo <- penguins %>%
    filter(species == "Gentoo")
print(gentoo)
```

Filtering and selecting in one with the pipe operator 
``` r
gentoo_bill_length <- penguins %>%
    filter(species == "Gentoo") %>%
    select(bill_length_mm)
print(gentoo_bill_length)
```
Now we understand a little bit of R we are prepared for some data analysis 
this afternoon. 

## Analysis of Microbiota data

Here we will perform a very streamlined analysis of shotgun metagenomic data from 
mouse stool. This dataset comes from our lab but is currently not published. 

## Data Analysis 

Source the "Microbiome_functions.R" script. This contains functions for analysing shotgun metagenomic data which we will use today. 

``` r
source("Microbiome_functions.R")
```

Firstly, we will read the species abundance table(s) and metadata into R, using
phyloseq. 

``` r
ps <- import_pseq_metag("data/mouse_stool_abundances.txt","data/mouse_stool_metadata.txt", level = "Species")
```

Phyloseq is used as it provides a nice way of storing all the associated
data in one object or class. As you can see, the metadata and species abundance
table are all combined here into `ps`. Importantly, each part can also
be accessed invidivually.

## Data Normalisation

Here will transform the data to relative abundance (percentages), and also 
perform a method called rarefaction for downstream analyses.

### Why do we need to normalise? 

``` r
barplot(sort(sample_sums(ps)), horiz = TRUE, las = 2, xlab=NULL, main="Library sizes")
```
![library_sizes](https://github.com/user-attachments/assets/73a51c10-84b2-4a4e-86bb-b76cabcbca9d)

### How does normalisation work?


Relative, transforms each sample into compositions, on a fixed scale of
0-1 or 0-100.

``` r
ps_rel <- transform(ps)
```

We can see how this works by looking at the column sums:

``` r
colSums(otu_table(ps_rel))
```
Rarefaction, on the other hand, works by choosing a fixed number of samples equal to or less than the sample with the lowest number of reads, then randomly discarding reads from larger samples until the remaining sample size equals this threshold.

``` r
ps_rar <- rarefy_even_depth(ps, rngseed = 42)
```
``` r
colSums(otu_table(ps_rar))
```

## Alpha Diversity

Here we will calculate two different measures of alpha-diversity:

-   Species richness, or the number of observed species
-   Shannon effective diversity, measuring richness and evenness 

The function `calc_alpha` wraps all of these calculations and only
requires the rarified phyloseq object as input, returning a
dataframe with a column for each above-mentioned dataframe.

``` r
alpha_div <- calc_alpha(ps_rar)
```

### How are alpha diversity metrics calculated?

Richness here is calculated as the the total number of observed species
greater than 0.5 mss normalised abundance. This threshold is used to
exclude noise (see Lagkouvardos et al 2017, PeerJ and Reitmeier et al
2021, ISME Communication for a more thorough explanation).

Shannon effective diversity is calculated as the exponent of the Shannon
index:

$$H = -\\sum\_{i=1}^{R} p\_i ln(p\_i)$$
where *R* = richness, *p<sub>i</sub>* is the relative abundance of
species *i* and *ln* is the natural logarithm.

This metric accounts for both the abundance and evenness of taxa.


### Plotting

To plot the alpha diversity metrics we will use a boxplot with jittered
points layered on top. The function `plot_boxplot`will do this for you,
we just need to set some parameters first.

Firstly, we will list the statistical comparisons we want to make, by
creating a list of vectors. In this case we only have two groups, stroke
and sham which can be written like: `list(c("Stroke", "Sham"))`. If we
had an extra group, for example, “Control”, we would then have three
comparisons to define like so:
`list(c("Control", "Sham"), c("Control", "Stroke"), c("Sham", "Stroke"))`


``` r
head(alpha_div)
```



Comparisons list

``` r
comparisons <- list(c("Stroke", "Sham"))
```

We can also specify the colours we want to use in our plots here by
creating a named vector of colours.

``` r
colour_pal <- c("sham" = "#98C1D9", "stroke"= "#EE6C4D")
```

To generate the plot we need to provide the dataframe, the name of the
column containing the grouping variable (in this case "condition”),
the name of the column containing the values to plot ("Richness"). To
colour by group we provide the column name of the grouping variable to
`fill_var`. We can then add the list of comparisons, some x and y-axis,
a title if desired, and the plot colours.

In instances where the alphabetical order of your group levels does not
match the order you would like to plot them in, you can specify the order
explicitly with the `group.order`parameter.

#### Richness

``` r
plot_boxplot(alpha_div, variable_col = "condition", value_col = "Richness", 
             comparisons_list = comparisons, fill_var = "condition", 
             group.order = c("sham", "stroke"), cols = colour_pal)
```
![richness](https://github.com/user-attachments/assets/6fbb9f26-d19f-48e1-81a5-b615e5770e93)


``` r
plot_boxplot(alpha_div, variable_col = "condition", value_col = "Shannon.Effective", 
             comparisons_list = comparisons, fill_var = "condition", 
             group.order = c("sham", "stroke"), cols = colour_pal)
```
![shannon_e](https://github.com/user-attachments/assets/0cb5fbaa-78a6-4909-82c9-ca0dfe34ed85)


Do stroke and sham mice show any differences in alpha-diversity? 

## Microbiome composition

In addition to examining within sample diversity (alpha-diversity), we also
usually want to know how different certain groups of samples are, in case, whether
sham and stroke differ in terms of microbiota composition. 

### Taxonomic composition - overview 

One very simple way to do this, which can give us a fairly rough picture of whether
our groups are different is to simply plot the taxonomic composition, looking
at the most abundant species. 

To plot composition we need to provide phyloseq object containing relative abundances, 
the taxonomic level we would like to plot, in this case species. In addition, 
we need to supply the group column (without quotes), 
the number of species we would like to plot and can additionally specify several additional 
parameters to customise the order and calculation of the most abundant taxa.
``` r
p <-
  plot_taxonomic_comp(
    ps_rel,
    "Species",
    condition,
    n_taxa = 12,
    ord = c("sham", "stroke"),
    per_group = T,
    groups = c("sham", "stroke")
  )
print(p)
```

![taxonomic_comp](https://github.com/user-attachments/assets/1c2e2b21-b381-401b-8059-431c53a41d95)

Is the taxonomic composition of stroke and sham mice different at the species level? 


### Beta Diversity

Here we will calculate beta-diversity based on Bray-Curtis
distance and plot an ordination of this using Non-metric
multidimensional scaling.

The `calc_betadiv` function calculates a distance matrix, and an
ordination of that matrix, returning both as a list.

Various dissimilarity indices are available:

-   Bray-Curtis - A count-based dissimilarity metric (beta-diversity),
    based on the fraction of overabundant counts.
-   Jaccard - Fraction of unique features. Does not consider abundance. 


Similarly, there are also various ordination options:

-   NMDS (default) - Non-Metric Multidimensional Scaling. An ordination
    method which attempts to represent the dissimilarity between
    samples, as closely as possible in a low-dimensional space.
-   MDS/PCoA - Principal Coordinate analysis (also known as Metric
    Multidimensional Scaling). An ordination method which attempts to
    preserve distance between samples in a low dimensional Euclidean
    space.

``` r
beta_div <- calc_betadiv(ps_rar, dist = "bray", ord_method = "NMDS")
```

To plot beta diversity a convenience function `plot-beta_div` is
provided. We just need to provide the phyloseq object and the betadiv 
object from above, a grouping variable and again the colour palette, if desired.

Within this function statistical testing of group separation is also
carried out using the adonis function of vegan. This function performs a
Permutational Multivariate Analysis of Variance or PERMANOVA test. The
resulting R<sup>2</sup> and p-value are added to the plot in the bottom
left.

The adonis R<sup>2</sup> represents the amount of variance in microbiome 
composition explained by the variable being tested, in this case condition. The
accompanying p-value denotes whether or not this effect is significant. 

``` r
plot_beta_div(ps_rar, beta_div, group_variable = "condition",cols = colour_pal, add_ellipse = T)
```
![betadiv](https://github.com/user-attachments/assets/b2533a80-6402-4a6d-968a-a893d4fdb138)

What do you interpret from this plot? 

## Differential abundance

The final step of this pipeline is to calculate differentially abundant
taxa between conditions.

This function performs the ancom-bc test, a compositionally aware
differential abundance method and returns significant results. As input,
only the phyloseq object and the column name of the grouping variable is
required.

``` r
da_taxa <- maaslin2_tax(
    ps_rel,
    out = "microbiome_DA_results",
    fixed = "condition",
    abun_thresh = 1e-4,
    prev_thresh = 0.3
)$results
```

To visualise differentially abundant taxa, we provide a function which
calculates fold change of significant taxa from above and plots
diverging dot plot coloured by group, providing a clear figure showing
which taxa change between conditions.

To this function, we need to provide the results of the ancom test
above, an ordered vector of the group levels e.g. ```c("Sham", "Stroke")```.
Additionally, we can provide the group colours to make interpretation
easier.

``` r
plot_da(da_taxa, groups = c("sham", "stroke"), cols = colour_pal)
```

![da_taxa](https://github.com/user-attachments/assets/7b481ab3-6355-4633-9d58-75a1f4d9fd7e)

What can we say about microbiome changes 3 days after experimental stroke?

