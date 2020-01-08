#######Notebook for R Tutorial (Wednesday, 25 December 2019)
######References: Grolemund
###Working Directory Query 
getwd()
###Matrix Filling from Vectors is autormatically by-column
#Compare 
hand2 <- c("ace", "spades", "king", "spades", "queen", "spades", "jack",   "spades", "ten", "spades")
matrix(hand2,nrow=5,byrow = TRUE)
#and 
matrix(hand2,nrow=5)
###Data Structures 
#Type of data: logical, double, characters(strings)
#Class: A special case of an atomic vector 

###Operations with Lists 
#Each data frame is a list with class data.frame 
typeof(df)
class(df)
str(df)

###Save Data
write.csv(deck, file = "cards.csv", row.names = FALSE)

###Different Ways of Indexing with deck[,] structure
#• Positive integers • Negative integers • Zero • Blank spaces • Logical values • Names
#Beside[,] indexing, data frames and lists admit '$' indexing 
#Subset a list with single brackets --> get back a list
#Subset a list with double brackets --> get back the value in the lists 

deal <- function(cards) {   
  cards[1, ] 
}
#Shuffling of rows of dataframe is possible 
shuffle <- function(cards) {   
  random <- sample(1:52, size = 52)   
  cards[random, ] 
}

###Logical Subsetting Examples pp.96 Grolemund
deck3$value[deck3$face == "ace"] <- 14 
# %in% is an operator. 

### Adding a column to a dataframe
df$columnname <- variable 

###Environments in R
#see the active environment 
environment()
#Scoping rules: R will look for the object in the active env, then the parent of the active env, until it reaches the empty env
#when a function is being called, we enter the runtime environment 

#Case1: Explain why the following function will not alter deck in the global environment:
deal <- function() {   
  card <- deck[1, ]   
  deck <- deck[-1, ]   
  card
}
#Answer: As the function is called, we enter the runtime environment. In the runtime environment, 
#the object 'deck' is created by deleting the first row of 'deck' accessed from the global environment
#the original object in the parent environment is not modified. 

#Case2: Explain what does the following function do?
shuffle <- function(){   
  random <- sample(1:52, size = 52)   
  assign("deck", DECK[random, ], envir = globalenv()) 
}
#The function creates an object random that is a permutation of the first 52 integers,
#Then the function accesses 'DECK' from the parent environment, shuffle it, and then assign this object
#to the original deck in the parent environment

### Vectorised Code (Chapter 10, Grolemund)
#Consider giving priority to element-wise operation. subsetting and logical tests

