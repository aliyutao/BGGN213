BGGN213\_Class7
================
Yutao Wen
4/24/2019

``` r
rescale <- function(x, na.rm=TRUE, plot=FALSE, ...) {
   rng <-range(x, na.rm=na.rm)
   answer <- (x - rng[1]) / (rng[2] - rng[1])
   if(plot) {
     plot(answer, ...)
}
   return(answer)
}
```

``` r
rescale(1:10)
```

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000

``` r
#rescale(c(1,10,"string"))
```

``` r
rescale2 <- function(x, na.rm=TRUE, plot=FALSE, ...) {
   if( !is.numeric(x) ) {
      stop("Input x should be numeric", call.=FALSE)
   }
   rng <-range(x, na.rm=na.rm)
   answer <- (x - rng[1]) / (rng[2] - rng[1])
   if(plot) {
      plot(answer, ...)
}
   return(answer)
}
```

``` r
#rescale2(c(1,10,"string"))
```

Function Practice
-----------------

Write a function to identify NA elements in two vectors.

Start with a simple example input where I know what the answer is:

``` r
x <- c(1,2,NA,3,NA)
y <- c(NA,3,NA,3,4)
```

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

Take sum of two to find how many of TRUE I have:

``` r
sum(is.na(x) & is.na(y))
```

    ## [1] 1

Create a function based on previous work:

``` r
both_na <- function(x,y) {
  sum(is.na(x) & is.na(y))
}
```

``` r
both_na(c(NA,NA,NA), c(NA,NA,NA))
```

    ## [1] 3

Test abnormal cases of different lengths of vectors:

``` r
both_na(c(NA,NA,NA), c(1,NA,NA,NA,NA,NA))
```

    ## [1] 5

The shorter vector was recylcled to match with the longer one.

Check the lengths of the two vectors beforehand:

``` r
x <- c(NA,NA,NA)
y <- c(1,NA,NA,NA,NA,NA)
length(x) != length(y)
```

    ## [1] TRUE

Add that to the function:

``` r
both_na2 <- function(x,y) {
  if(length(x) != length(y)) {
    stop("Input x and y should be the same length")
  }
  sum(is.na(x) & is.na(y))
}
```

``` r
#both_na2(x,y)
```

Try the both\_na3 function with extra features:

``` r
both_na3 <- function(x, y) {
if(length(x) != length(y)) {
stop("Input x and y should be vectors of the same length")

}
  na.in.both <- ( is.na(x) & is.na(y) )
  na.number  <- sum(na.in.both)
  na.which   <- which(na.in.both)
  message("Found ", na.number, " NA's at position(s):",
          paste(na.which, collapse=", ") )
  return( list(number=na.number, which=na.which) )
}
```

``` r
#both_na3(x,y)
```

``` r
X <- c( 1, 2, NA, 3, NA) 
Y <- c(NA,3,NA,3, 4)
both_na3(X,Y)
```

    ## Found 1 NA's at position(s):3

    ## $number
    ## [1] 1
    ## 
    ## $which
    ## [1] 3

Write a function grade() to determine an overall grade from a vector of student homework assignment scores dropping the lowest single alignment score.

``` r
# Student 1
s1 <- c(100,100,100,100,100,100,100,90)
# Student 2
s2 <- c(100, NA, 90, 90, 90, 90, 97, 80)
```

``` r
grade <- function(x) {
  avg <- (sum(x, na.rm = TRUE)-min(x, na.rm = TRUE)) / (length(x)-1)
  return((avg))
}
```

``` r
grade(s1)
```

    ## [1] 100

``` r
grade(s2)
```

    ## [1] 79.57143

``` r
grade(c(s1,s2))
```

    ## [1] 89.8

``` r
url <- read.csv("https://tinyurl.com/gradeinput", row.names = 1)
ans <- apply(url, 1,grade)
sorted <- sort(ans, decreasing = TRUE)
sorted[1:5]
```

    ##  student-7  student-8 student-13  student-1 student-12 
    ##      94.00      93.75      92.25      91.75      91.75

``` r
# Start with a simple version of the problem
df1 <- data.frame(IDs=c("gene1", "gene2", "gene3"),
                  exp=c(2,1,1),
                  stringsAsFactors=FALSE)
df2 <- data.frame(IDs=c("gene2", "gene4", "gene3", "gene5"),
                  exp=c(-2, NA, 1, 2),
                  stringsAsFactors=FALSE)

# Simplify further to single vectors
x <- df1$IDs
y <- df2$IDs

intersect(x,y)
```

    ## [1] "gene2" "gene3"

``` r
x%in%y
```

    ## [1] FALSE  TRUE  TRUE

``` r
x[x%in%y]
```

    ## [1] "gene2" "gene3"

``` r
cbind( x[ x %in% y ], y[ y %in% x ] )
```

    ##      [,1]    [,2]   
    ## [1,] "gene2" "gene2"
    ## [2,] "gene3" "gene3"

``` r
gene_interest <- function(x,y) {
  cbind( x[ x %in% y ], y[ y %in% x ] )
}
```

``` r
merge(df1, df2, by="IDs")
```

    ##     IDs exp.x exp.y
    ## 1 gene2     1    -2
    ## 2 gene3     1     1

``` r
if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
```

    ## Loading required namespace: BiocManager

``` r
BiocManager::install()
```

    ## Bioconductor version 3.8 (BiocManager 1.30.4), R 3.5.3 (2019-03-11)

    ## Update old packages: 'boot', 'cluster'
